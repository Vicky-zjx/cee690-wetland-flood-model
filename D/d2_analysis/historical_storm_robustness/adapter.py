"""Private observed-event input provider for the unchanged Stage D runner."""
from __future__ import annotations
import copy
from pathlib import Path
import sys
from types import FunctionType
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
from D.d2_analysis.parameter_uncertainty import adapter as parameter_adapter
from D.feedback_control import run_policy as frozen
import pandas as pd

CONDITIONING = ('Selected historical storm with unchanged shared physical initial states, fixed boundary conditions, '
                'wetland design, Stage D gate hydraulics and 120 s feedback policy; matched saved parameter vector. '
                'An experiment with historical forcing, not a reconstruction of event-specific antecedent storage.')


def resolve(path):
    p = Path(path)
    return p if p.is_absolute() else ROOT / p


def run_one(point, event, case, dt=15.):
    runner, records = parameter_adapter._private_runner(point['physical_values'], dt)
    source = resolve(event['forcing_file'])
    if frozen.sha(source) != event['forcing_sha256']:
        raise ValueError('Historical forcing changed after event selection')
    if not event['benchmark']:
        namespace = dict(runner.__globals__)
        manifest = dict(scenario_id=event['event_id'], forcing_file=str(source),
                        forcing_sha256=event['forcing_sha256'], seed_sequence_entropy=None)
        namespace['original_weather'] = lambda: [copy.deepcopy(manifest)]
        private = FunctionType(runner.__code__, namespace, runner.__name__, runner.__defaults__, runner.__closure__)
        private.__kwdefaults__ = copy.copy(runner.__kwdefaults__)
        runner = private
    result = runner(case=case, scenario_id=None if event['benchmark'] else event['event_id'],
                    design='best', dt=float(dt), hydraulics='D', Cw_US=2.8)
    if len(records) != 1: raise RuntimeError('Expected one private physical configuration')
    d = result['diagnostics']; record = copy.deepcopy(records[0])
    if event['benchmark'] and record['prepared_forcing_sha256'] != event['forcing_sha256']:
        raise ValueError('Benchmark manifest does not identify the internally used original forcing')
    if pd.Timestamp(d['event_start_utc']) != pd.Timestamp(event['start_utc']):
        raise ValueError('Applied weather start does not match the selected event')
    if (pd.Timestamp(event['end_utc'])-pd.Timestamp(event['start_utc'])).total_seconds() != d['fixed_configuration']['duration_s']:
        raise ValueError('Selected event does not cover the unchanged full horizon')
    record['nominal_2022_template_provenance'] = {k:record[k] for k in (
        'deterministic_weather_id', 'prepared_forcing_path', 'prepared_forcing_sha256',
        'prepared_forcing_manifest_sha256', 'adapter_sha256')}
    record.update(deterministic_weather_id=event['event_id'], prepared_forcing_path=str(source.relative_to(ROOT)),
        prepared_forcing_sha256=event['forcing_sha256'],
        prepared_forcing_manifest_sha256=frozen.sha(HERE/'events/manifest.json'),
        adapter_sha256=frozen.sha(__file__), parameter_template_adapter_sha256=frozen.sha(parameter_adapter.__file__),
        conditioning=CONDITIONING, event_start_utc=event['start_utc'], event_end_utc=event['end_utc'],
        forcing_method='unchanged original 2022 prepared knots' if event['benchmark'] else
                       'unchanged Stage C hourly rows; existing conservative hourly P/Rs and center-linear T/wind adapter',
        weather_seed=None, weather_source_kind='historical observations and existing grid-estimated solar radiation')
    d['parameter_evaluation'] = record
    d['fixed_configuration']['parameter_evaluation'] = copy.deepcopy(record)
    d['historical_event'] = copy.deepcopy(event)
    d['fixed_configuration']['conditional_on'] = CONDITIONING
    expected = record['initial_wetland_m3' if case == 'wetland' else 'initial_no_wetland_m3']
    if d['initial_storage_m3'] != expected: raise RuntimeError('Original physical initial state changed')
    if record['effective_physical_values'] != point['physical_values']:
        raise RuntimeError('Applied parameter vector differs from the saved draw')
    return result
