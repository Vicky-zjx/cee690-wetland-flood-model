"""Opt-in Stage D runs with frozen forcing, gates, design and initial water.

No run occurs on import. A revised RHS is supplied by the local model.make_rhs
factory; only its private recovered-adapter namespace is changed. A/B/C and
previous D modules are read-only dependencies. All weather IDs remain recorded,
including the four unchanged cold-domain exclusions. No fitting or optimization.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import sys
import time
import traceback
from types import SimpleNamespace

sys.dont_write_bytecode = True
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/cee690-stage-d-revision-matplotlib')
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import numpy as np
import pandas as pd
from D.g5_diagnostics import runner as historical
from D.gate_hydraulics_revision import preservation

RECOVERED = ROOT/'D/g5_diagnostics/recovered'
ADAPTER = RECOVERED/'model_adapters.py'
SCENARIOS = RECOVERED/'scenarios'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(historical.json_clean(value), indent=2, allow_nan=False)+'\n')


def weather_manifest():
    """Verify existing bytes and metadata; do not generate or select by outcome."""
    manifest = json.loads((RECOVERED/'scenario_manifest.json').read_text())
    records=[]
    for item in manifest['scenarios']:
        path=RECOVERED/item['forcing_file']
        side=json.loads(path.with_suffix('.json').read_text())
        digest=sha(path)
        if digest != item['forcing_sha256'] or digest != side['forcing_sha256']:
            raise ValueError('Original weather checksum mismatch: '+item['scenario_id'])
        if item['seed_sequence_entropy'] != side['seed_sequence_entropy']:
            raise ValueError('Original seed metadata mismatch: '+item['scenario_id'])
        frame,_=historical.adapter().read_forcing(path)
        cold=frame['temperature_C']<=0
        records.append(dict(scenario_id=item['scenario_id'],method=item['method'],
            seed_sequence_entropy=item['seed_sequence_entropy'],forcing_file=str(path),
            forcing_sha256=digest,rows=len(frame),start_utc=frame['timestamp_utc'].iloc[0].isoformat(),
            end_exclusive_utc=(frame['timestamp_utc'].iloc[-1]+pd.Timedelta(hours=1)).isoformat(),
            minimum_temperature_C=float(frame['temperature_C'].min()),
            nonpositive_temperature_hours=int(cold.sum()),original_liquid_domain_admissible=not bool(cold.any())))
    if len(records)!=20 or len({r['scenario_id'] for r in records})!=20:
        raise ValueError('Expected the unchanged twenty original weather IDs')
    if sum(r['original_liquid_domain_admissible'] for r in records)!=16:
        raise ValueError('The original 16 warm / 4 cold accounting changed')
    return records


def private_adapter(variant, hydraulic_diagnostics):
    """Independent function globals plus module proxies; no B monkey-patching."""
    from D.gate_hydraulics_revision import model
    spec=importlib.util.spec_from_file_location('_stage_d_revision_adapter', ADAPTER)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.ROOT=ROOT
    module.wet=SimpleNamespace(**vars(module.wet))
    module.dry=SimpleNamespace(**vars(module.dry))
    module.wet.make_b2_rhs=lambda gates, regimes=None: model.make_rhs(
        'wetland',gates,variant,diagnostics=hydraulic_diagnostics)
    module.dry.make_rhs=lambda gates, diagnostics=None: model.make_rhs(
        'no_wetland',gates,variant,diagnostics=hydraulic_diagnostics)
    return module



def _revised_failure_reporting(result, start, hydraulic_diagnostics):
    """Name the revised physical domain, not the adapter's legacy G5 wrapper."""
    d=result['diagnostics']
    violation=d.get('first_violation')
    physical=hydraulic_diagnostics.get('first_applicability_failure')
    if physical is not None:
        record=dict(violation or physical)
        reason=physical.get('reason') or physical.get('regime')
        record['reason']=reason
        record['regime']=reason
        record['legacy_G5_residual_m']=record.pop('residual_m',record.get('legacy_G5_residual_m'))
        record['legacy_G5_role']='Descriptive comparison with Stage B; not the revised failure criterion'
        record['detection']='First encountered official ode4u RHS/RK evaluation, not a continuous-time crossing estimate'
        if reason=='broad_crest_head_length_outside_supported_range':
            lo,hi=record['required_range'];value=record['head_length_ratio']
            record['applicability_residual']=max(lo-value,value-hi)
            record['applicability_residual_units']='1'
            record['applicability_rule']='0.05 <= upstream head / crest length <= 0.5 for this branch'
        elif reason=='tailwater_at_or_above_sill':
            sill=record['upstream_stage_m']-record['water_depth_m']
            record['applicability_residual']=record['tailwater_stage_m']-sill
            record['applicability_residual_units']='m'
            record['applicability_rule']='Tailwater must be strictly below gate sill in this bounded variant'
        elif reason=='transition_energy_inconsistent':
            record['applicability_rule']='Discharge must not exceed the independently calculated opening and ideal critical capacities'
        d.update(classification='revised_hydraulic_applicability_failure',
            failure_reason=reason,legacy_adapter_failure_reason=d.get('failure_reason'),
            failure_stage='ode_rhs_revised_hydraulic_applicability',
            failure_message='Stage D hydraulic applicability failure: '+str(reason),
            first_violation=record)
        hydraulic_diagnostics['first_applicability_failure']=dict(record)
        violation=record
    first=(d.get('all_rhs_observations') or {}).get('first_rhs_exception') or {}
    failure_time=(violation or {}).get('time_s',first.get('time_s'))
    if not result['valid'] and failure_time is not None:
        utc=(pd.Timestamp(start)+pd.to_timedelta(failure_time,unit='s')).isoformat()
        d.update(first_failure_time_s=float(failure_time),first_failure_timestamp_utc=utc)
        if violation is not None:
            violation['timestamp_utc']=utc
            violation['time_s_meaning']=('First excluded hourly input interval; no integration was started'
                if d.get('failure_reason')=='nonfreezing_liquid_water_domain' else
                'First evaluated failing RHS/RK stage; not an exact continuous-time crossing')
        if first:
            first['timestamp_utc']=utc
            first['time_s_meaning']='First evaluated failing RHS/RK stage'
            if physical is not None:
                first['physical_failure_reason']=physical.get('reason')
                first['diagnostic']=dict(violation)


def _saved_baseline_comparison(result, case, design):
    """Compare changed physics on saved 120 s nodes; never call it regression."""
    relative=('B/results/wetland_counterfactual/bypass_120s.npz' if case=='no_wetland' else
              'B/results/wetland_only/'+('fixed_reference_120s.npz' if design=='reference' else 'wetland_best_120s.npz'))
    path=ROOT/relative
    comparison=dict(kind='Structural hydraulic-model comparison, not unchanged-model regression',
        saved_trajectory=relative,saved_sha256=sha(path),expected_identical_trajectory=False,
        comparison_available=False,
        comparison_grid='Original saved 120 s nodes; new fine-grid peaks use the new native grid')
    with np.load(path,allow_pickle=False) as saved:
        oldtime=saved['time']; old_y=saved['y']
        names=list(historical.adapter().dry.OUTPUT_NAMES if case=='no_wetland' else historical.adapter().a2.OUTPUT_NAMES)
        river=old_y[names.index('H_r')]; outlet=old_y[names.index('Q_out')]
        comparison.update(saved_physical_river_peak_m=float(river.max()),
            saved_physical_peak_outlet_m3_s=float(outlet.max()),saved_output_nodes=len(oldtime))
        if result['valid']:
            indices=np.searchsorted(result['time'],oldtime)
            if np.any(indices>=len(result['time'])) or not np.array_equal(result['time'][indices],oldtime):
                comparison['reason_unavailable']='Original saved timestamps do not exactly align; no interpolation invented'
            else:
                differences={key:float(np.max(np.abs(result[key][...,indices]-saved[key])))
                             for key in ('x','dx','y','u') if key in saved}
                comparison.update(comparison_available=True,maximum_absolute_differences=differences,
                    identical_initial_states=bool(np.array_equal(result['x'][:,0],saved['x'][:,0])),
                    identical_inputs_at_original_nodes=bool(np.array_equal(result['u'][:,indices],saved['u'])),
                    river_peak_change_m=float(np.max(result['downstream_stage'])-np.max(river)),
                    peak_outlet_change_m3_s=float(np.max(result['downstream_discharge'])-np.max(outlet)),
                    max_absolute_stage_differences_m={name:float(np.max(abs(result['y'][names.index(name),indices]-old_y[names.index(name)])))
                        for name in names if name.startswith('H_')},
                    missing_saved_columns=[key for key in ('dx',) if key not in saved])
        else:
            comparison['reason_unavailable']='Revised trajectory is incomplete or invalid; no revised full-event peak is inferred'
    comparison['saved_file_unchanged']=sha(path)==comparison['saved_sha256']
    result['diagnostics']['saved_B_baseline_comparison']=comparison


def run_case(case, *, dt_s=120, scenario_id=None, design='best', variant=None, original=False, weather=None):
    """One deterministic or saved-weather run; fixed gate schedule is mandatory.

    model.make_rhs(case, gates, variant, diagnostics=None) must return the same
    eight/six-state RHS and named output schema. All timestamps use original UTC.
    """
    if case not in ('wetland','no_wetland') or dt_s not in (120,60,30):
        raise ValueError('Expected known physical case and 120/60/30 s resolution')
    variant={} if variant is None else dict(variant)
    sid=scenario_id or 'deterministic_september_2022'
    record=None
    if scenario_id is not None:
        candidates=weather_manifest() if weather is None else weather
        record=next((r for r in candidates if r['scenario_id']==scenario_id),None)
        if record is None:raise ValueError('Unknown original scenario ID: '+scenario_id)
        if sha(record['forcing_file'])!=record['forcing_sha256']:
            raise ValueError('Weather changed since execution planning: '+scenario_id)
    if original:
        if scenario_id is None:
            result=historical.run_deterministic(case,dt=dt_s,gate_scale=1,design=design,constant_opening=None)
        else:
            result=historical.run_saved(case,record['forcing_file'],dt=dt_s,gate_scale=1,design=design,
                scenario_id=sid,seed=record['seed_sequence_entropy'],constant_opening=None)
        result['diagnostics']['hydraulic_variant']='original_B'
        return result

    hydraulic_diagnostics={}
    module=private_adapter(variant,hydraulic_diagnostics)
    fixed=historical.load_fixed(dt=dt_s,gate_scale=1,design=design,constant_opening=None)
    if not np.array_equal(fixed['fixed_gates'],fixed['manifest']['original_fixed_gate_values']):
        raise ValueError('Revised hydraulic runs must retain every original gate knot')
    fixed['manifest'].update(hydraulic_variant='Stage D opt-in hydraulic revision',
        hydraulic_variant_settings=variant,
        stage_d_source_sha256={str(path.relative_to(ROOT)):sha(path) for path in
            (HERE/'run_study.py',HERE/'model.py',HERE/'hydraulics.py') if path.exists()},
        legacy_G5_note='The old H>=1.25a margin is diagnostic only for the revised law; new physical applicability is enforced inside the revised gate evaluator.')
    observer=historical.RHSObserver(case,fixed)
    fixed['official']['ode4u']=observer.wrap_solver(fixed['official']['ode4u'])
    start=fixed['model_config']['geography']['event_start_utc']
    try:
        if scenario_id is None:
            params=module.wet.base_parameters(fixed['model_config'])
            source_time,source_inputs=module.wet.fixed_forcing(fixed['model_config'],params)
            target=params.time
            inputs=module.dry.inputs(fixed['model_config'],fixed['fixed_gates'],target,(source_time,source_inputs))
            diagnostics=dict(forcing_mapping='Original deterministic linear weather knots; no hourly remapping',
                fixed_configuration=fixed['manifest'],event_start_utc=start,
                forcing_conservation=dict(source_precipitation_m=float(np.trapezoid(source_inputs[0],source_time)),
                    mapped_precipitation_m=float(np.trapezoid(inputs[module.a2.U['P']],target)),
                    source_solar_Wh_m2=float(np.trapezoid(source_inputs[2],source_time)/3600),
                    mapped_solar_Wh_m2=float(np.trapezoid(inputs[module.a2.U['R_s']],target)/3600)))
            result=module._integrate(case,[(target,inputs)],fixed,diagnostics)
        else:
            start=record['start_utc']
            result=module.run_case(case,record['forcing_file'],fixed)
        result=historical._finish(result,observer,sid,start,dt_s,1,design,
            seed=None if record is None else record['seed_sequence_entropy'])
    except Exception as exc:
        # Preserve one failure record and continue other predetermined cases.
        # No physical objective or full-horizon safety is invented here.
        diagnostics=dict(scenario_id=sid,case_name=case,design_source_label=design,dt_s=dt_s,
            original_seed=None if record is None else record['seed_sequence_entropy'],
            classification='unexpected_execution_error',failure_reason='unexpected_execution_error',
            failure_stage='input_or_execution_exception',exception_type=type(exc).__name__,
            failure_message=str(exc),exception_traceback=traceback.format_exc(),
            fixed_configuration=fixed['manifest'],all_rhs_observations=observer.summary(),
            ode_started=bool(observer.segment_calls),full_valid_trajectory_available=False,
            full_integration_horizon_completed=False,safety_violation=None,
            physical_river_peak_m=None,physical_river_peak_discharge_m3_s=None)
        if record is not None:
            diagnostics.update(forcing_file_sha256=record['forcing_sha256'],forcing_file=record['forcing_file'])
        result=module._empty_result(case,'execution_error',diagnostics)
    _revised_failure_reporting(result,start,hydraulic_diagnostics)
    d=result['diagnostics']
    d.update(hydraulic_variant='Stage D opt-in hydraulic revision',hydraulic_variant_settings=variant,
             hydraulic_diagnostics=hydraulic_diagnostics,
             original_G5_margin_role='Historical diagnostic, not the revised-law rejection rule')
    if result['valid']:
        k=int(np.argmax(result['downstream_discharge']))
        d.update(physical_river_peak_discharge_m3_s=float(result['downstream_discharge'][k]),
                 physical_river_peak_discharge_time_s=float(result['time'][k]))
    else:
        d.setdefault('physical_river_peak_m',None)
        d.setdefault('physical_river_peak_discharge_m3_s',None)
    if scenario_id is None:
        _saved_baseline_comparison(result,case,design)
    return result


def save_result(result, directory, name):
    """Preserve unsafe completed trajectories; invalid runs have no physical peak."""
    d=result['diagnostics'];directory.mkdir(parents=True,exist_ok=True)
    write_json(directory/(name+'.json'),d)
    if result['valid']:
        np.savez_compressed(directory/(name+'.npz'),**{k:result[k] for k in
            ('time','x','dx','y','u','balance','downstream_stage','downstream_discharge','state_names','output_names')})
    row=historical.failure_record(result)
    ob=d.get('all_rhs_observations',{});fc=d.get('forcing_conservation',{})
    row.update(hydraulic_variant=d.get('hydraulic_variant'),
        ode_started=d.get('ode_started',False),full_output_samples=len(result['time']),
        full_horizon_completed=d.get('full_integration_horizon_completed',False),
        safety_violation=d.get('safety_violation'),
        safety_violation_any_successful_rhs=d.get('safety_violation_any_successful_rhs'),
        physical_river_peak_discharge_m3_s=(float(np.max(result['downstream_discharge'])) if result['valid'] else None),
        balance_relative=d.get('balance_relative'),balance_max_m3=d.get('balance_max_m3'),
        derivative_mass_identity_m3_s=ob.get('derivative_physical_mass_identity_max_abs_m3_s'),
        forcing_sha256=d.get('forcing_file_sha256'),
        rain_mapping_error_m=fc.get('precipitation_error_m'),solar_mapping_error_Wh_m2=fc.get('solar_error_Wh_m2'),
        first_exception=(ob.get('first_rhs_exception') or {}).get('exception_type'))
    violation=d.get('first_violation') or {}
    row.update(first_failure_time_s=d.get('first_failure_time_s'),
        first_failure_timestamp_utc=d.get('first_failure_timestamp_utc'),
        legacy_G5_residual_m=violation.get('legacy_G5_residual_m'),
        applicability_residual=violation.get('applicability_residual'),
        applicability_residual_units=violation.get('applicability_residual_units'),
        applicability_rule=violation.get('applicability_rule'))
    for safety in ob.get('safety_successful_rhs',[]):
        row[safety['variable']+'_peak_m']=safety['peak_m']
        row[safety['variable']+'_margin_m']=safety['margin_m']
    return row


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['plan','run'])
    parser.add_argument('--hydraulics',choices=['original','revised'],default='revised')
    parser.add_argument('--suite',choices=['deterministic','weather','all'],default='all')
    parser.add_argument('--dt',nargs='+',type=int,choices=[120,60,30],default=[120])
    parser.add_argument('--variant-config',type=Path)
    parser.add_argument('--tag')
    args=parser.parse_args()
    variant=json.loads(args.variant_config.read_text()) if args.variant_config else {}
    if not isinstance(variant,dict):raise ValueError('Variant configuration must be a JSON object')
    check=preservation.verify()
    if not check['passed']:raise RuntimeError('Protected files changed: '+json.dumps(check))
    weather=weather_manifest()
    digest=hashlib.sha256(json.dumps(variant,sort_keys=True).encode()).hexdigest()[:10]
    tag=args.tag or f'{args.hydraulics}_{args.suite}_{digest}'
    if not re.fullmatch(r'[A-Za-z0-9_-]+',tag):raise ValueError('Use a simple output tag')
    output=HERE/'results'/tag
    if output.exists():raise FileExistsError('Use a new tag; existing results are preserved: '+str(output))
    plans=[]
    for dt in args.dt:
        if args.suite in ('deterministic','all'):
            plans.extend(dict(scenario_id=None,case=c,design=d,dt_s=dt) for c,d in
                [('wetland','reference'),('wetland','best'),('no_wetland','best')])
        if args.suite in ('weather','all'):
            plans.extend(dict(scenario_id=r['scenario_id'],case=c,design='best',dt_s=dt)
                         for r in weather for c in ('wetland','no_wetland'))
    output.mkdir(parents=True)
    write_json(output/'execution_plan.json',dict(action=args.action,hydraulics=args.hydraulics,
        variant=variant,variant_config_sha256=sha(args.variant_config) if args.variant_config else None,
        scenario_manifest_sha256=sha(RECOVERED/'scenario_manifest.json'),weather=weather,plans=plans,
        original_gate_schedule_required=True,weather_regenerated=False,
        fixed_configuration=historical.load_fixed()['manifest']))
    if args.action=='plan':
        print(json.dumps(dict(planned_cases=len(plans),simulations_run=0,output=str(output))))
        return 0
    rows=[]
    for plan in plans:
        sid=plan['scenario_id'] or 'deterministic_september_2022'
        name=f'{sid}__{plan["case"]}__{plan["design"]}__{plan["dt_s"]}s'
        print('START',name,flush=True);start=time.monotonic()
        try:
            result=run_case(**plan,variant=variant,original=args.hydraulics=='original',weather=weather)
            row=save_result(result,output,name)
        except Exception as exc:
            error=dict(**plan,classification='unexpected_execution_error',valid=False,
                exception_type=type(exc).__name__,failure_message=str(exc),exception_traceback=traceback.format_exc())
            write_json(output/(name+'.error.json'),error)
            row={k:v for k,v in error.items() if k!='exception_traceback'}
            row['case_name']=plan['case'];row['scenario_id']=sid
        row['elapsed_s']=time.monotonic()-start;rows.append(row)
        pd.DataFrame(rows).to_csv(output/'scenario_diagnostics.csv',index=False)
        print('DONE',name,row.get('classification'),flush=True)
    counts=pd.DataFrame(rows).groupby(['case_name','classification'],dropna=False).size().to_dict()
    write_json(output/'counts.json',{f'{k[0]} / {k[1]}':int(v) for k,v in counts.items()})
    check=preservation.verify();write_json(output/'preservation_check.json',check)
    if not check['passed']:raise RuntimeError('Protected files changed during execution')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
