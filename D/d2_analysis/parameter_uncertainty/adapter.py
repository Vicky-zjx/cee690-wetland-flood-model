"""Paired deterministic parameter evaluations with the frozen Stage D policy.

The caller supplies physical values; this module defines no distributions,
correlations, random draws, parameter multipliers or calibrated values. Each
case gets a fresh configuration and parameter object. Original modules, their
cached configurations and all historical files remain unchanged. Results are
returned in memory; only the caller chooses whether and where to save them.
"""
from __future__ import annotations

import copy
from collections.abc import Mapping
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys
from types import FunctionType, SimpleNamespace

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
from D.feedback_control import run_policy as frozen

UNITS = {'K_u': 'm/s', 'k_ud': '1/s (two branches)', 'k_ur': '1/s', 'R_r': 's/m2'}
CONDITIONING = ('Fixed September 2022 event weather and boundaries, physical initial storages, '
                'historical B-selected wetland design, Stage D hydraulics and 120 s feedback policy; '
                'only the explicitly recorded physical parameters vary. No uncertainty law is implied.')


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _normalise(physical_values):
    if not isinstance(physical_values, Mapping):
        raise TypeError('physical_values must be a mapping of explicit physical parameter values')
    unknown = set(physical_values) - set(UNITS)
    if unknown:
        raise ValueError('Only K_u, k_ud, k_ur and R_r may vary: ' + ', '.join(sorted(unknown)))
    values = {}
    for key, value in physical_values.items():
        if isinstance(value, (bool, np.bool_)):
            raise ValueError(key + ' requires a numerical physical value, not a Boolean')
        array = np.asarray(value, dtype=float)
        expected = (2,) if key == 'k_ud' else ()
        if array.shape != expected or not np.all(np.isfinite(array)):
            raise ValueError(key + ' requires finite values with shape ' + str(expected))
        if np.any(array < 0 if key == 'K_u' else array <= 0):
            raise ValueError(key + (' must be nonnegative' if key == 'K_u' else ' must be positive'))
        values[key] = array.tolist() if key == 'k_ud' else float(array)
    return values


def nominal_values():
    """Read the actual existing values without running a model or drawing data."""
    module = frozen.historical.adapter()
    fixed = frozen.historical.load_fixed(dt=15, gate_scale=1, design='best', constant_opening=None)
    params, _, _ = module._parameters(fixed)
    return {key:np.asarray(getattr(params, key)).tolist() for key in UNITS}


def _private_runner(physical_values, dt):
    values = _normalise(physical_values)
    if isinstance(dt, (bool, np.bool_)) or dt not in (15, 7.5):
        raise ValueError('This study supports the existing 15 s resolution and its verified 7.5 s refinement')
    base = frozen.run_case
    refinement_metadata = None
    if dt == 7.5:
        from D.d2_analysis.run_refinement import refined_runner
        base, refinement_metadata = refined_runner()
    records = []
    historical = frozen.historical
    module = historical.adapter()

    def load_private_fixed(*args, **kwargs):
        fixed = historical.load_fixed(*args, **kwargs)
        # The historical loader already deep-copies its cached configuration.
        # A second explicit copy makes this adapter's ownership unambiguous.
        fixed['model_config'] = copy.deepcopy(fixed['model_config'])
        config = fixed['model_config']
        for key, value in values.items():
            config['parameter_overrides'][key] = dict(
                value=copy.deepcopy(value), source_id='STAGE_D_EXPLICIT_PARAMETER_ARGUMENT',
                approved=True,
                evidence_category='Explicit caller-supplied evaluation point; not a measurement or fitted distribution',
                reason='Authorized software evaluation of an explicit physical parameter argument; no probability law, calibration or source inference is implied')
        params, initial8, gates = module._parameters(fixed)
        effective = {key:np.asarray(getattr(params, key)).tolist() for key in UNITS}
        forcing_path = module.wet.previous.resolve_path(config['forcing']['prepared_csv'])
        forcing_manifest = module.wet.previous.resolve_path(config['forcing']['manifest_json'])
        vector_id = hashlib.sha256(json.dumps(effective, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        record = dict(
            parameter_vector_sha256=vector_id, requested_physical_overrides=copy.deepcopy(values),
            effective_physical_values=effective, units=dict(UNITS),
            unspecified_values_rule='Retain existing nominal values; never sample or infer a value',
            soil_capacity_m3=float(params.C_m), soil_field_capacity_m3=float(params.S_fc),
            derived_initial_soil_fraction=float(params.soil_fraction_init),
            initial_wetland_m3=initial8.tolist(), initial_no_wetland_m3=initial8[module.dry.SHARED].tolist(),
            additional_initial_wetland_water_m3=float(initial8[2:4].sum()),
            gate_parameters=[asdict(gate) for gate in gates],
            initial_actual_opening_m=[0.0015, 0.0015], controller_period_s=float(frozen.CONTROL_PERIOD_S),
            policy_sha256=_sha(frozen.POLICY_PATH),
            design_categories=list(fixed['pair']), design_covers=fixed['covers'].tolist(),
            design_label='Historical B-selected wetland design; no new optimization',
            deterministic_weather_id='deterministic_september_2022',
            prepared_forcing_path=str(forcing_path.relative_to(ROOT)), prepared_forcing_sha256=_sha(forcing_path),
            prepared_forcing_manifest_sha256=_sha(forcing_manifest),
            boundaries=copy.deepcopy(config['forcing']['proposed_fixed_boundaries']),
            dt_s=float(dt), hydraulic_arm='D', Cw_US=2.8,
            conditioning=CONDITIONING,
            adapter_method='Private function globals with a per-call loader proxy; no original function/module/cache mutation',
            adapter_sha256=_sha(__file__), frozen_driver_sha256=_sha(frozen.__file__),
            refinement_metadata=copy.deepcopy(refinement_metadata))
        fixed['manifest']['parameter_evaluation'] = copy.deepcopy(record)
        records.append(record)
        return fixed

    proxy = SimpleNamespace(adapter=historical.adapter, load_fixed=load_private_fixed,
                            RHSObserver=historical.RHSObserver)
    namespace = dict(base.__globals__)
    namespace['historical'] = proxy
    runner = FunctionType(base.__code__, namespace, base.__name__, base.__defaults__, base.__closure__)
    runner.__kwdefaults__ = copy.copy(base.__kwdefaults__)
    return runner, records


def evaluate_pair(physical_values, dt=15):
    """Return both full results for one explicit, shared parameter vector.

    Allowed fields are K_u [m/s], k_ud [two 1/s values], k_ur [1/s], and
    R_r [s/m2, resistance]. Omitted fields retain their existing nominal value.
    The same fixed prepared storm, boundaries, initial physical storages and
    feedback policy are used in both models. Actual gate motions remain the
    policy's state-dependent outputs; no schedule is substituted or optimized.
    Completed unsafe simulations and any failures remain in their result objects.
    No files are written and no statistical distribution is chosen here.
    """
    values = _normalise(physical_values)
    outputs = {}
    records = {}
    for case in ('wetland', 'no_wetland'):
        runner, loaded = _private_runner(values, dt)
        result = runner(case=case, scenario_id=None, design='best', dt=float(dt), hydraulics='D', Cw_US=2.8)
        if len(loaded) != 1:
            raise RuntimeError('Expected one private physical configuration per model evaluation')
        record = loaded[0]
        diagnostics = result['diagnostics']
        diagnostics['parameter_evaluation'] = record
        diagnostics['fixed_configuration']['conditional_on'] = CONDITIONING
        if diagnostics['initial_storage_m3'] != record['initial_wetland_m3' if case == 'wetland' else 'initial_no_wetland_m3']:
            raise RuntimeError('A model did not use the recorded fixed physical initial state')
        records[case] = record
        outputs[case] = result
    for field in ('effective_physical_values', 'parameter_vector_sha256', 'policy_sha256',
                  'prepared_forcing_sha256', 'boundaries', 'gate_parameters', 'design_categories', 'design_covers'):
        if records['wetland'][field] != records['no_wetland'][field]:
            raise RuntimeError('Paired models differ in a required common setting: ' + field)
    return outputs
