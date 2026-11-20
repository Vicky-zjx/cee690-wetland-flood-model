"""Run only the predeclared antecedent experiment, preserving prior results.

Twelve physical pairs: 2 historical events x 2 prefixes x 3 saved vectors.
A single longest-horizon pair is additionally refined at 7.5 s as a numerical
check, not a new parameter draw or a thirteenth physical experiment.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import copy
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
from D.feedback_control import run_policy as frozen
from D.d2_analysis.parameter_uncertainty import run_study as prior
from D.d2_analysis.antecedent_conditions import adapter, preservation

EVENT_IDS = ['historical_20190825T0200Z', 'historical_20240804T1500Z']
SAMPLE_IDS = ['mc_00016', 'mc_00024', 'mc_00018']
PILOT = ('historical_20240804T1500Z', 72, 'mc_00018')


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    frozen.write_json(path, value)


def key(event, hours, point):
    return f"{event['event_id']}__{hours}h__{point['sample_id']}"


def prepare():
    source = ROOT/'D/d2_analysis/historical_storm_robustness/study_plan.json'
    previous = json.loads(source.read_text())
    if previous['pilot_sample_ids'] != SAMPLE_IDS:
        raise ValueError('Predeclared parameter IDs changed')
    points = [next(p for p in previous['points'] if p['sample_id'] == sid) for sid in SAMPLE_IDS]
    events = [next(e for e in previous['events'] if e['event_id'] == eid) for eid in EVENT_IDS]
    input_paths = [source, ROOT/'D/d2_analysis/robustness_checks/FINDINGS.md',
                   HERE/'forcing/manifest.json', HERE/'adapter.py', HERE/'forcing/provider.py',
                   Path(frozen.__file__), frozen.POLICY_PATH, HERE/'run_study.py']
    plan = dict(
        experiment='Predeclared 24 h / 72 h continuous antecedent conditioning',
        events=events, points=points, antecedent_hours=[24,72],
        design_selection='Existing low/central/high mean-log-multiplier pilot IDs, selected before these simulations; not benefit or severity ranks',
        primary_dt_s=15., refined_dt_s=7.5, control_period_s=120.,
        unique_physical_pairs=12, primary_model_runs=24, numerical_refinement_model_runs=2,
        pilot=dict(event_id=PILOT[0], antecedent_hours=PILOT[1], sample_id=PILOT[2],
            selection='Longest 72 h prefix, latest event, highest predeclared mean-log point; selected before observing new outcomes'),
        pilot_tolerances={k:v for k,v in previous['pilot_tolerances'].items() if k!='selection'}, numerics=previous['numerics'],
        existing_no_antecedent_results='Reuse original saved pilot arrays at 15 s; no original trajectory rerun',
        initial_conditions='Accepted physical states, initial openings 0.0015 m and original enabled latches at prefix start; no reset at event boundary',
        event_weather='Original 144 h solver weather inputs unchanged; separately mapped prefix with explicit boundary endpoint holds',
        parameter_sampling_calls=0, weather_generation_calls=0,
        source_sha256={str(p.relative_to(ROOT)): frozen.sha(p) for p in input_paths},
        source_parameter_plan_sha256=previous['source_parameter_plan_sha256'])
    dest = HERE/'study_plan.json'
    if dest.exists():
        if json.loads(dest.read_text()) != plan: raise ValueError('Existing frozen experimental plan differs')
    else: write(dest, plan)
    return plan


def save(result, point, event, hours, dt):
    directory = HERE/'results'/key(event,hours,point)
    directory.mkdir(parents=True, exist_ok=True)
    name = prior.stem(point['sample_id'], result['case_name'], dt)
    dest = directory/(name+'.json')
    if dest.exists(): raise FileExistsError(dest)
    result['diagnostics']['antecedent_study_plan_sha256'] = frozen.sha(HERE/'study_plan.json')
    write(dest, result['diagnostics'])
    prior.history(result).to_csv(directory/(name+'_controller.csv.gz'), index=False, compression='gzip')
    if all(k in result for k in prior.ARRAY_KEYS):
        np.savez_compressed(directory/(name+'.npz'), **{k:result[k] for k in prior.ARRAY_KEYS})
    if 'block_audit' in result:
        pd.DataFrame(result['block_audit']).to_csv(directory/(name+'_blocks.csv.gz'),index=False,compression='gzip')
    return directory


def evaluate(point, event, hours, dt):
    results = {}
    for case in prior.CASES:
        try: results[case] = adapter.run_one(point,event,hours,case,dt)
        except Exception as exc: results[case] = prior._failure(case,exc,dt)
        save(results[case],point,event,hours,dt)
    row = prior.paired_row(point,dt,results)
    row.update(event_id=event['event_id'],antecedent_hours=hours,metric_scope='Entire continuous antecedent plus 144 h event; event-window metrics in separate analysis tables')
    return results,row


def convergence(coarse, fine, point, event, hours, plan):
    records=[]
    module=frozen.historical.adapter();fixed=frozen.historical.load_fixed(dt=15,design='best')
    params,_,_=module._parameters(fixed)
    den8=np.r_[params.A_u,params.A_u,params.A_w,params.A_d,params.S_y*params.A_g,params.A_r]
    tol=plan['pilot_tolerances']; duration=(144+hours)*3600.; boundary=hours*3600.
    for case in prior.CASES:
        a,b=coarse[case],fine[case]
        row=dict(event_id=event['event_id'],sample_id=point['sample_id'],antecedent_hours=hours,
            case_name=case,coarse_dt_s=15,refined_dt_s=7.5,both_valid=a['valid'] and b['valid'],passed=False)
        if not row['both_valid']:
            row['reason']='At least one trajectory invalid; inspect retained diagnostics';records.append(row);continue
        ix=np.searchsorted(b['time'],a['time']); clocks=np.array_equal(a['time'],b['time'][ix])
        den=den8 if case=='wetland' else den8[module.dry.SHARED]
        hi=[i for i,n in enumerate(a['output_names']) if str(n).startswith('H_')]
        states=float(np.max(np.abs(a['x']-b['x'][:,ix])/den[:,None]))
        stages=float(np.max(np.abs(a['y'][hi]-b['y'][hi][:,ix])))
        peaks=float(np.max(np.abs(a['y'][hi].max(axis=1)-b['y'][hi].max(axis=1))))
        gate=float(np.max(np.abs((a['u'][:2]-b['u'][:2,ix])*.3)))
        ha,hb=prior.history(a),prior.history(b);hc=prior._history_compare(ha,hb)
        count=int(duration/120);clock_ok=bool(len(ha)==count*2 and np.array_equal(ha.time_s.to_numpy()[::2],np.arange(count)*120))
        weather=bool(np.array_equal(a['u'][2:],b['u'][2:,ix]))
        mass=all(r['diagnostics']['balance_relative']<plan['numerics']['mass_balance_relative_tolerance'] for r in (a,b))
        eventpeak=abs(a['downstream_stage'][a['time']>=boundary].max()-b['downstream_stage'][b['time']>=boundary].max())
        row.update(common_grid_exact=bool(clocks),controller_clock_exact=clock_ok,weather_inputs_exact=weather,
            max_state_equivalent_depth_difference_m=states,max_stage_difference_m=stages,max_all_stage_native_peak_difference_m=peaks,
            max_actual_opening_difference_m=gate,latch_mismatch_count=hc['latch_mismatch_count'],
            event_peak_stage_difference_m=float(eventpeak),mass_pass=mass,
            coarse_unsafe=a['diagnostics'].get('safety_violation'),fine_unsafe=b['diagnostics'].get('safety_violation'))
        row['passed']=bool(clocks and clock_ok and weather and mass and hc['clock_equal'] and hc['latch_mismatch_count']==0
            and states<=tol['max_state_equivalent_depth_difference_m'] and stages<=tol['max_stage_difference_m']
            and peaks<=tol['max_physical_peak_stage_difference_m'] and eventpeak<=tol['max_physical_peak_stage_difference_m']
            and gate<=tol['max_actual_gate_opening_difference_m'])
        records.append(row)
    return records


def run_job(point,event,hours,refine,plan):
    start=time.monotonic(); print('START',key(event,hours,point),flush=True)
    coarse,row=evaluate(point,event,hours,15.)
    rows=[row];checks=[]
    if refine:
        fine,other=evaluate(point,event,hours,7.5);rows.append(other)
        checks=convergence(coarse,fine,point,event,hours,plan)
    directory=HERE/'results'/key(event,hours,point)
    pd.DataFrame(rows).to_csv(directory/'paired_results.csv',index=False)
    if checks:
        pd.DataFrame(checks).to_csv(directory/'numerical_comparison.csv',index=False)
        write(HERE/'results/pilot_summary.json',dict(passed=all(r['passed'] for r in checks),checks=checks,
            extra_parameter_draws=0,physical_pair_count=1,numerical_model_executions=4))
    print('DONE',key(event,hours,point),row['paired_valid'],row['wetland_classification'],row['no_wetland_classification'],round(time.monotonic()-start,2),flush=True)
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--phase',required=True,choices=['prepare','pilot','remaining'])
    p.add_argument('--workers',type=int,choices=[1,2,3,4],default=3)
    args=p.parse_args();guard=preservation.verify()
    if not guard['passed']: raise RuntimeError('Protected prior files changed')
    plan=prepare()
    if args.phase=='prepare': print('Prepared immutable 12-pair plan; no ODE runs');return
    if args.phase=='pilot':
        e=next(e for e in plan['events'] if e['event_id']==PILOT[0]);q=next(q for q in plan['points'] if q['sample_id']==PILOT[2])
        run_job(q,e,PILOT[1],True,plan)
    else:
        if not json.loads((HERE/'results/pilot_summary.json').read_text())['passed']:
            raise RuntimeError('Numerical pilot failed; remaining experiments not started')
        jobs=[(q,e,h) for e in plan['events'] for h in plan['antecedent_hours'] for q in plan['points']
              if (e['event_id'],h,q['sample_id'])!=PILOT]
        with ProcessPoolExecutor(max_workers=args.workers,mp_context=get_context('spawn')) as pool:
            pending=[pool.submit(run_job,q,e,h,False,plan) for q,e,h in jobs]
            for f in as_completed(pending): f.result()
    write(HERE/'preservation_check.json',preservation.verify())


if __name__=='__main__':main()
