"""Historical-storm pilot and matched reuse of the existing 64 parameter draws.

No new parameter sampling, fitting, optimization or physical equation changes.
"""
from __future__ import annotations
import argparse
import copy
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from multiprocessing import get_context
from pathlib import Path
import sys
import time
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
import numpy as np
import pandas as pd
from D.d2_analysis.parameter_uncertainty import run_study as prior
from D.d2_analysis.historical_storm_robustness import adapter, preservation
from D.feedback_control import run_policy as frozen

OLD = ROOT/'D/d2_analysis/parameter_uncertainty'
CASES = prior.CASES


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    frozen.write_json(path, value)


def prepare():
    source = OLD/'results/study_plan.json'
    old = json.loads(source.read_text())
    events_source = HERE/'events/manifest.json'
    events = json.loads(events_source.read_text())['events']
    if len(events) not in (3,4,5) or sum(bool(e['benchmark']) for e in events) != 1:
        raise ValueError('Require 3–5 selected events including exactly one original benchmark')
    scores = [(float(np.log([p[k] for k in ('m_K_u','m_routing','m_R_r')]).mean()), p['sample_id']) for p in old['mc']]
    median = float(np.median([s[0] for s in scores]))
    identifiers = [min(scores)[1], min(scores, key=lambda s:(abs(s[0]-median),s[1]))[1], max(scores)[1]]
    points = copy.deepcopy(old['mc'])
    plan = dict(source_parameter_plan=str(source.relative_to(ROOT)), source_parameter_plan_sha256=frozen.sha(source),
        events_manifest_sha256=frozen.sha(events_source), points=points, events=events,
        source_seed=old['seed'], source_distribution=old['config']['monte_carlo'], parameter_sampling_calls=0,
        pilot_sample_ids=identifiers, pilot_selection='Lowest, nearest median, highest mean log multiplier among original 64 vectors; not risk ranks',
        pilot_tolerances=old['config']['pilot'], numerics=old['config']['numerics'],
        primary_dt_s=15., refined_dt_s=7.5, original_initial_conditions_preserved=True,
        no_pooling_interpretation='Selected observed events have no occurrence probabilities; the same parameter rows are dependent across events',
        benchmark_reuse='Existing 64-pair 2022 outputs reused exactly after one complete trajectory adapter regression',
        runtime_source_sha256={p.name:frozen.sha(p) for p in (HERE/'run_study.py', HERE/'adapter.py')})
    dest = HERE/'study_plan.json'
    if dest.exists():
        saved = json.loads(dest.read_text())
        if saved != plan: raise ValueError('Study configuration/code differs from the saved immutable plan')
        return saved
    with dest.open('x') as stream: json.dump(plan, stream, indent=2); stream.write('\n')
    return plan


def evaluate(point, event, dt):
    results = {}
    for case in CASES:
        try: results[case] = adapter.run_one(point, event, case, dt)
        except Exception as exc: results[case] = prior._failure(case, exc, dt)
    return results


def extend_metrics(result):
    row = prior.metrics(result); d = result['diagnostics']
    row.update(completed=d.get('full_integration_horizon_completed'),
        balance_max_m3=d.get('balance_max_m3'), solver_termination=d.get('solver_termination'),
        initial_total_storage_m3=sum(d.get('initial_storage_m3',[])) or None,
        final_total_storage_m3=None, total_storage_change_m3=None, external_net_integrated_m3=None,
        end_storage_balance_residual_m3=None)
    if result['valid']:
        initial = float(result['x'][:,0].sum()); final = float(result['x'][:,-1].sum())
        residual = float(result['balance'][-1])
        row.update(initial_total_storage_m3=initial, final_total_storage_m3=final,
            total_storage_change_m3=final-initial, end_storage_balance_residual_m3=residual,
            external_net_integrated_m3=final-initial-residual,
            external_net_method='Recovered from unchanged per-control-block balance; not trapezoid of joined discontinuous rain inputs')
    return row


def save_case(result, point, event, dt, directory, full):
    case = result['case_name']; name = prior.stem(point['sample_id'], case, dt)
    d = result['diagnostics']
    d['historical_storm_study'] = dict(event_id=event['event_id'], sample_id=point['sample_id'],
        study_plan_sha256=frozen.sha(HERE/'study_plan.json'), sample_physical_values=point['physical_values'],
        reused_from=result.get('reused_from'), frozen_feedback_source_sha256=frozen.sha(frozen.__file__))
    metrics = extend_metrics(result); d['historical_storm_metrics'] = metrics
    write(directory/(name+'.json'), d)
    prior.history(result).to_csv(directory/(name+'_controller.csv.gz'), index=False, compression='gzip')
    if result['valid']:
        if full:
            arrays = {k:result[k] for k in prior.ARRAY_KEYS}
        else:
            record = d['parameter_evaluation']
            max_open = np.array([g['gate_max_opening'] for g in record['gate_parameters']])[:,None]
            names = list(result['output_names']); indices = [i for i,n in enumerate(names) if str(n).startswith('H_')]
            arrays = {k:result[k] for k in ('time','x','balance','downstream_stage','downstream_discharge','state_names')}
            arrays.update(actual_opening_m=result['u'][:2]*max_open,
                stages=result['y'][indices], stage_names=np.array(names)[indices],
                external_net_integrated_m3=result['x'].sum(axis=0)-result['x'][:,0].sum()-result['balance'])
        np.savez_compressed(directory/(name+'.npz'), **arrays)
    return metrics


def run_point(point, event, phase, plan):
    directory = HERE/'results'/phase/event['event_id']; rows=[]; models=[]; by_dt={}
    for dt in ([15.,7.5] if phase=='pilot' else [15.]):
        print('START',phase,event['event_id'],point['sample_id'],dt,flush=True); start=time.monotonic()
        olddir = HERE/'results/pilot'/event['event_id']
        if phase=='ensemble' and point['sample_id'] in plan['pilot_sample_ids']:
            results={c:prior._read_saved(olddir,prior.stem(point['sample_id'],c,dt),c) for c in CASES}
        else: results=evaluate(point,event,dt)
        row=prior.paired_row(point,dt,results); row['event_id']=event['event_id']
        for case,result in results.items():
            metrics=save_case(result,point,event,dt,directory,full=phase=='pilot')
            models.append(dict(event_id=event['event_id'],sample_id=point['sample_id'],case_name=case,dt_s=dt,**metrics))
            row.update({case+'_'+key:value for key,value in metrics.items()})
        row['elapsed_s']=time.monotonic()-start;rows.append(row);by_dt[dt]=results
        print('DONE',event['event_id'],point['sample_id'],dt,row['paired_valid'],
              row['wetland_classification'],row['no_wetland_classification'],round(row['elapsed_s'],2),flush=True)
    checks=prior.convergence(point,by_dt[15.],by_dt[7.5],
        dict(pilot=plan['pilot_tolerances'],numerics=plan['numerics'])) if phase=='pilot' else []
    for check in checks:check['event_id']=event['event_id']
    return rows,models,checks


def regression(plan):
    output=HERE/'results/regression'
    if output.exists():raise FileExistsError(output)
    output.mkdir(parents=True)
    point=plan['points'][0];event=next(e for e in plan['events'] if e['benchmark'])
    results=evaluate(point,event,15.);records=[]
    for case,result in results.items():
        path=OLD/f'results/pilot/pilot_mc_00000__{case}__15s.npz'
        same={}
        if result['valid']:
            with np.load(path,allow_pickle=False) as previous:
                for key in ('time','x','dx','y','u','balance','downstream_stage','downstream_discharge'):
                    same[key]=bool(np.array_equal(previous[key],result[key]))
        records.append(dict(case_name=case,valid=result['valid'],exact_arrays=same,
            passed=bool(result['valid'] and all(same.values())),reference_path=str(path.relative_to(ROOT)),reference_sha256=frozen.sha(path)))
        save_case(result,point,event,15.,output,full=True)
    write(output/'regression_summary.json',dict(passed=all(r['passed'] for r in records),checks=records,
        purpose='Exact original September 2022 forcing/parameter/controller regression before new events'))


def reuse_benchmark(event, plan):
    directory=HERE/'results/ensemble'/event['event_id'];directory.mkdir(parents=True)
    source=OLD/'results/monte_carlo'
    pair=pd.read_csv(source/'paired_results.csv',float_precision='round_trip');pair['event_id']=event['event_id']
    if pair.sample_id.tolist()!=[p['sample_id'] for p in plan['points']]:raise ValueError('Benchmark sample order changed')
    pair['historical_benchmark_reused']=True
    pair.to_csv(directory/'paired_results.csv',index=False)
    model=pd.read_csv(source/'model_results.csv',float_precision='round_trip');model['event_id']=event['event_id']
    for point in plan['points']:
        for case in CASES:
            name=prior.stem(point['sample_id'],case,15.);path=source/(name+'.json');d=json.loads(path.read_text())
            d['historical_event']=copy.deepcopy(event)
            d['historical_benchmark_reuse']=dict(source=str(path.relative_to(ROOT)),sha256=frozen.sha(path),
                controller_source=str((source/(name+'_controller.csv.gz')).relative_to(ROOT)),
                primary_metrics_source=str((source/'paired_results.csv').relative_to(ROOT)),
                note='Exact previously saved numerical outputs; no new trajectory fabricated or repeated parameter draw')
            write(directory/(name+'.json'),d)
    model.to_csv(directory/'model_results.csv',index=False)
    write(directory/'counts.json',json.loads((source/'counts.json').read_text())|
        dict(event_id=event['event_id'],reused_existing_benchmark=True,source_csv_sha256=frozen.sha(source/'paired_results.csv')))


def execute(phase, plan, workers):
    output=HERE/'results'/phase
    if output.exists():raise FileExistsError(output)
    if not json.loads((HERE/'results/regression/regression_summary.json').read_text())['passed']:
        raise RuntimeError('Original deterministic regression did not pass')
    if phase=='ensemble' and not json.loads((HERE/'results/pilot/pilot_summary.json').read_text())['passed']:
        raise RuntimeError('Historical-event numerical/hydraulic pilot did not pass; review failures without relaxing criteria')
    output.mkdir(parents=True)
    events=[e for e in plan['events'] if not e['benchmark']]
    points=[p for p in plan['points'] if phase!='pilot' or p['sample_id'] in plan['pilot_sample_ids']]
    for event in events:(output/event['event_id']).mkdir()
    if phase=='ensemble':reuse_benchmark(next(e for e in plan['events'] if e['benchmark']),plan)
    write(output/'execution_plan.json',dict(study_plan_sha256=frozen.sha(HERE/'study_plan.json'),workers=workers,
        event_ids=[e['event_id'] for e in events],sample_ids=[p['sample_id'] for p in points],
        all_planned_cases_retained=True,safety_is_reported_separately_from_model_validity=True))
    rows=[];models=[];checks=[]
    def collect(value):
        a,b,c=value;rows.extend(a);models.extend(b);checks.extend(c)
        eid=a[0]['event_id'];directory=output/eid
        pd.DataFrame([r for r in rows if r['event_id']==eid]).sort_values(['sample_id','dt_s']).to_csv(directory/'paired_results.csv',index=False)
        pd.DataFrame([r for r in models if r['event_id']==eid]).sort_values(['sample_id','dt_s','case_name']).to_csv(directory/'model_results.csv',index=False)
        if c:pd.DataFrame([r for r in checks if r['event_id']==eid]).to_csv(directory/'numerical_comparison.csv',index=False)
    with ProcessPoolExecutor(max_workers=workers,mp_context=get_context('spawn')) as pool:
        pending=[pool.submit(run_point,p,e,phase,plan) for e in events for p in points]
        for future in as_completed(pending):collect(future.result())
    for event in events:
        selected=[r for r in rows if r['event_id']==event['event_id']]
        write(output/event['event_id']/'counts.json',dict(event_id=event['event_id'],planned_parameter_vectors=len(points),
            paired_evaluations=len(selected),valid_pairs=sum(r['paired_valid'] for r in selected),
            classifications={c:pd.Series([r[c+'_classification'] for r in selected]).value_counts().to_dict() for c in CASES}))
    if phase=='pilot':
        write(output/'pilot_summary.json',dict(passed=bool(checks and all(r['passed'] for r in checks)),
            selected_sample_ids=plan['pilot_sample_ids'],new_historical_events=len(events),checks=checks,
            safety_exceeded_cases=sum(bool(r.get('unsafe')) for r in models),
            interpretation='Numerical/hydraulic applicability pilot; valid unsafe cases are reported and never discarded'))
    check=preservation.verify();write(output/'preservation_check.json',check)
    if not check['passed']:raise RuntimeError('Protected files changed during execution')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phase',required=True,choices=['prepare','regression','pilot','ensemble'])
    parser.add_argument('--workers',type=int,default=4,choices=[1,2,3,4]);args=parser.parse_args()
    check=preservation.verify()
    if not check['passed']:raise RuntimeError('Protected previous files changed: '+json.dumps(check))
    plan=prepare()
    if args.phase=='prepare':print('Prepared fixed historical-event and unchanged64-vector plan; no random draws')
    elif args.phase=='regression':
        regression(plan)
        check=preservation.verify();write(HERE/'results/regression/preservation_check.json',check)
        if not check['passed']:raise RuntimeError('Protected files changed during regression')
    else:execute(args.phase,plan,args.workers)


if __name__=='__main__':main()
