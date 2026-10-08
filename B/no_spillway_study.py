"""Bounded revised B.2 study; all historical scientific files/results are read-only.

Run --screen, --run, --compare, then --report. Official NMS and ode4u only.
"""
from pathlib import Path
from collections import Counter, OrderedDict
import argparse, itertools, json, time
import numpy as np
from B import optimize_wetland as historical
from B import no_spillway_model as model
from B.official_opt import run_opt, load_official
from B.multistart_check import save_trajectory

ROOT=historical.ROOT
OUT=ROOT/'results/no_spillway_revision'
CONFIG=ROOT/'no_spillway_config.json'
CALLBACKS=100
BATCH_SECONDS=1500
TOTAL_ODE_LIMIT=1300
MAP='odds_cube'


def to_design(q):
    q=np.asarray(q,dtype=float); x=q.copy()
    a=q[:-2]**3; b=(1-q[:-2])**3
    x[:-2]=a/(a+b)
    return x


def to_coordinates(x):
    x=np.asarray(x,dtype=float); q=x.copy()
    a=np.cbrt(x[:-2]); b=np.cbrt(1-x[:-2]); q[:-2]=a/(a+b)
    return q


def load_trajectory(path):
    with np.load(historical.resolve_path(path)) as d:
        result={k:d[k].copy() for k in d.files if k!='metadata_json'}
        result.update(json.loads(str(d['metadata_json'])))
    return result


def save(name,data):historical.save_json(OUT/name,data)

def assert_frozen(config):
    if model.preflight(config):raise RuntimeError('Unapproved revised configuration')
    lock=json.loads((OUT/'configuration_lock.json').read_text())['configuration']
    if config != lock:raise RuntimeError('Configuration changed since revised-study lock; no mixed-case execution.')


def candidates(reference):
    """Deterministic finite pool; retain all existing 24 gate knot times."""
    ref=np.asarray(reference); pool=[('revised_reference',ref.copy())]
    patterns=[
        ('timing_cover_30_80',[.8,.8,1.25,1.25,1.05,.8,.95,1.08,1.08,1.05,1.0,1.0],
         [.75,.75,1.1,1.1,.9,.8,1.04,1.08,1.08,1.03,1.0,1.0],[.30,.80]),
        ('asymmetric_cover_80_35',[1.1,1.1,.85,.95,1.05,.9,1.06,1.05,1.04,1.0,1.02,1.0],
         [.9,.9,1.1,1.05,1.0,1.1,.96,1.06,1.04,1.03,1.0,1.0],[.80,.35]),
        ('gentle_timing_cover_30_80',[.95,.95,1.03,1.03,1.01,.98,1.01,1.02,1.02,1.01,1.,1.],
         [.96,.96,1.02,1.02,.98,.98,1.01,1.02,1.02,1.01,1.,1.],[.30,.80]),
        ('gentle_asymmetric_cover_80_35',[1.02,1.02,.98,.99,1.01,.99,1.02,1.02,1.01,1.,1.,1.],
         [.98,.98,1.02,1.02,1.,1.02,1.01,1.02,1.02,1.01,1.,1.],[.80,.35]),
    ]
    for label,f1,f2,covers in patterns:
        x=ref.copy(); x[:12]*=f1;x[12:24]*=f2;x[-2:]=covers;pool.append((label,x))
    for j,amp in enumerate([.01,-.01,.003,-.003,.001]):
        x=ref.copy()
        x[:12]*=1+amp*np.array([1,1,-1,-1,1,0,1,1,-1,0,0,0])
        x[12:24]*=1+amp*np.array([-1,-1,1,1,0,1,-1,1,1,0,0,0])
        x[-2:]=[.30,.80] if j%2==0 else [.80,.35]
        pool.append((f'small_timing_{j+1}',x))
    return pool


def screen(config):
    if (OUT/'start_screening.json').exists():raise RuntimeError('Screening exists; not overwriting.')
    if model.preflight(config):raise RuntimeError('Unapproved configuration')
    ref=load_trajectory(OUT/'reference.npz')
    forcing=model.fixed_forcing(config,model.base_parameters(config)); records=[]; accepted=[]
    pool=candidates(ref['theta'])
    save('screening_plan.json',dict(maximum_candidates_per_pair=10,candidates=[dict(label=l,theta=x) for l,x in pool],
        coordinate_map=MAP,formula='g=q^3/(q^3+(1-q)^3); cover=q',full_physical_bounds=[0,1],
        rationale='One fixed bijection for every gate and start; accommodates both dry-period small openings and larger storm releases without narrowing physical bounds.'))
    for pair in itertools.product(*config['categories']['admissible']):
        found=0
        for label,theta in pool:
            row=dict(category_pair=list(pair),start_id=label,requested_design=theta)
            if np.any(theta<0) or np.any(theta>1):
                row.update(status='native_bound_rejection',feasible=False,ode_attempted=False)
            else:
                q=to_coordinates(theta); effective=to_design(np.clip(q,.1,.9))
                row.update(solver_coordinates=q,expected_first_design=effective)
                if not np.allclose(effective,theta,rtol=0,atol=3e-15):
                    row.update(status='initializer_would_shift',feasible=False,ode_attempted=False)
                else:
                    row['ode_attempted']=True
                    try:
                        result=model.simulate(theta,pair,config,forcing=forcing)
                        ok=model.is_feasible(result,config)
                        row.update(status='feasible' if ok else 'valid_constraint_violation',feasible=ok,
                                   metrics=result['metrics'],audit=result['audit'])
                        if ok:
                            save_trajectory(OUT/('_'.join(pair)+'_'+label+'_screen.npz'),result)
                            accepted.append(dict(category_pair=list(pair),start_id=label,theta=theta));found+=1
                    except (historical.TrajectoryFailure,FloatingPointError) as exc:
                        row.update(status='invalid',feasible=False,failure=str(exc),diagnostic=getattr(exc,'diagnostic',None))
            records.append(row);print('Screen',pair,label,row['status'],flush=True)
            if found==3:break
    save('start_screening.json',dict(records=records,accepted=accepted,maximum_candidates_per_pair=10,
        actual_ode_attempts=sum(r['ode_attempted'] for r in records)))


def run_batch(config):
    assert_frozen(config)
    if (OUT/'batch_status.json').exists():raise RuntimeError('Completed batch exists; not overwriting.')
    record=json.loads((OUT/'start_screening.json').read_text());starts=record['accepted']
    reference_record=json.loads((OUT/'reference_search.json').read_text())
    previous_attempts=reference_record.get('actual_ode_attempts',len(reference_record['records']))+record['actual_ode_attempts']
    forcing=model.fixed_forcing(config,model.base_parameters(config))
    settings={**config['optimizer'],'max_evals':CALLBACKS}
    official=load_official()
    save('batch_plan.json',dict(configuration=config,configuration_sha256=historical.config_hash(config),settings=settings,
        requested_runs=len(starts),maximum_callbacks=CALLBACKS*len(starts),maximum_total_ode_attempts=TOTAL_ODE_LIMIT,
        maximum_seconds=BATCH_SECONDS,coordinate_map=MAP,full_physical_bounds=[0,1],full_wrapper_coordinate_bounds=[0,1],
        native_official_bounds=[-1,1],constraint_and_objective_scaling_changed=False,
        reference='Recomputed no-dam-spillway A/A reference, fixed before joint optimization',
        official_source={k:v for k,v in official.items() if k not in ['nms','ode4u','opt_hyp']}))
    began=time.monotonic();total_odes=previous_attempts;reports=[];global_best=None;abort=False
    for start in starts:
        theta=np.asarray(start['theta']);pair=start['category_pair'];label=start['start_id'];slug='_'.join(pair)+'_'+label
        records=[];cache=OrderedDict();best=None;first=[];odes=0;before=time.monotonic();solver_result={}
        def callback(q,constants=None):
            nonlocal best,odes
            if time.monotonic()-began>BATCH_SECONDS or total_odes+odes>=TOTAL_ODE_LIMIT-40:
                raise historical.RunBudgetExceeded('Declared batch budget; reserved fixed-design checks.')
            design=to_design(q)
            if not first:
                first.extend(design.tolist())
                if not np.allclose(design,theta,rtol=0,atol=3e-15):raise RuntimeError('First physical design moved')
            key=design.tobytes();cached=key in cache
            if cached: objective,residuals,result,error,diagnostic=cache[key]
            else:
                odes+=1
                try:
                    result=model.simulate(design,pair,config,forcing=forcing)
                    objective=(result['metrics']['peak_stage_m']-config['numerics']['objective_origin_m'])/config['numerics']['height_scale_m']
                    residuals=result['residuals'];error=None;diagnostic=None
                except (historical.TrajectoryFailure,FloatingPointError) as exc:
                    result=None;objective=config['numerics']['failure_objective']
                    residuals=np.full(historical.CONSTRAINT_COUNT,config['numerics']['failure_residual'])
                    error=type(exc).__name__+': '+str(exc);diagnostic=getattr(exc,'diagnostic',None)
                cache[key]=(objective,residuals,result,error,diagnostic)
                if len(cache)>8:cache.popitem(last=False)
            feasible=result is not None and model.is_feasible(result,config)
            if feasible and (best is None or result['metrics']['peak_stage_m']<best['metrics']['peak_stage_m']):best=result
            records.append(dict(evaluation=len(records)+1,category_pair='/'.join(pair),cached=cached,
                raw_peak_stage_m=None if result is None else result['metrics']['peak_stage_m'],conditioned_objective=float(objective),
                max_scaled_violation=float(max(0,max(residuals))),feasible=feasible,
                best_feasible_peak_m=None if best is None else best['metrics']['peak_stage_m'],failure=error or '',
                diagnostic=diagnostic,hydraulically_invalid=bool(error and 'HydraulicApplicabilityFailure' in error),
                decision_json=json.dumps(design.tolist())))
            return float(objective),residuals.copy()
        try:
            solver_result=run_opt(callback,to_coordinates(theta),np.zeros(26),np.ones(26),settings)
            termination=solver_result['status']
        except historical.RunBudgetExceeded as exc:termination='batch_budget: '+str(exc);abort=True
        total_odes+=odes
        if best is not None:
            save_trajectory(OUT/(slug+'_best.npz'),best)
            if global_best is None or best['metrics']['peak_stage_m']<global_best['metrics']['peak_stage_m']:global_best=best
        report=dict(category_pair=pair,start_id=label,requested_initial_design=theta,first_optimizer_design=first,
            first_design_feasible=bool(records and records[0]['feasible']),optimizer_coordinate_settings=solver_result,
            physical_callback_log=records,termination=termination,callbacks=len(records),cached_callbacks=sum(r['cached'] for r in records),
            actual_ode_attempts=odes,seconds=time.monotonic()-before,
            rejection_counts=dict(Counter(r['failure'].split(':',1)[0] for r in records if r['failure'])),
            best_design=None if best is None else best['theta'],best_metrics=None if best is None else best['metrics'],
            best_audit=None if best is None else best['audit'])
        save(slug+'_run.json',report)
        reports.append({k:v for k,v in report.items() if k not in ['physical_callback_log','optimizer_coordinate_settings','best_audit']})
        print('OPT',pair,label,termination,'peak',None if best is None else best['metrics']['peak_stage_m'],'ODE',odes,flush=True)
        if abort:break
    if global_best is not None:save_trajectory(OUT/'batch_best.npz',global_best)
    save('batch_status.json',dict(status='bounded_batch_stopped' if abort else 'bounded_batch_finished',runs=reports,
        seconds=time.monotonic()-began,ode_attempts_including_reference_and_screening=total_odes,
        global_best_pair=None if global_best is None else global_best['pair'],
        global_best_design=None if global_best is None else global_best['theta']))


def compare(config):
    assert_frozen(config)
    plan=json.loads((OUT/'batch_plan.json').read_text())
    if historical.config_hash(config)!=plan['configuration_sha256']:raise RuntimeError('Configuration differs from optimizer batch.')
    if (OUT/'controlled_comparisons.json').exists():raise RuntimeError('Comparison exists; not overwriting.')
    forcing=model.fixed_forcing(config,model.base_parameters(config));ref=load_trajectory(OUT/'reference.npz')
    designs=[('revised_reference',ref)]
    if (OUT/'batch_best.npz').exists():designs.append(('joint_best',load_trajectory(OUT/'batch_best.npz')))
    # Preserve a former solution only as a fixed design rerun under the NEW model, not an old outcome.
    prior=load_trajectory(ROOT/'results/multistart/batch_best.npz')
    designs.append(('historical_design_recomputed',prior))
    rows=[]; fixed_cache={}
    def evaluate(label,theta,pair,dt):
        try:
            key=(dt,tuple(pair),np.asarray(theta,dtype=float).tobytes())
            for old_key, (old_run, _) in fixed_cache.items():
                if old_key[:2]==key[:2] and np.allclose(old_run['theta'],theta,rtol=0,atol=3e-15):key=old_key;break
            reused=key in fixed_cache
            if reused:r,first_label=fixed_cache[key]
            else:
                r=model.simulate(theta,pair,config,refine=int(120/dt),forcing=forcing)
                first_label=label;fixed_cache[key]=(r,label)
            filename=f'{label}_{dt}s.npz';save_trajectory(OUT/filename,r)
            row=dict(label=label,dt_s=dt,design=theta,pair=pair,status='valid',feasible=model.is_feasible(r,config),
                metrics=r['metrics'],audit=r['audit'],trajectory_file=filename,ode_attempted=not reused,simulation_reused_from=first_label if reused else None,
                rain_depth_m=float(np.trapezoid(r['u'][historical.a2.U['P']],r['time'])),
                solar_energy_Wh_m2=float(np.trapezoid(r['u'][historical.a2.U['R_s']],r['time'])/3600))
        except (historical.TrajectoryFailure,FloatingPointError) as exc:
            row=dict(label=label,dt_s=dt,design=theta,pair=pair,status='invalid',feasible=False,ode_attempted=True,
                     failure=str(exc),diagnostic=getattr(exc,'diagnostic',None))
        rows.append(row);print('Resolution',label,dt,row['status'],row['feasible'],flush=True)
        return row
    for label,r in designs:
        for dt in [120,60]:evaluate(label,r['theta'],r['pair'],dt)
    comparisons=[];need30=False
    for label,_ in designs[1:]:
        group={r['dt_s']:r for r in rows if r['label']==label};reference={r['dt_s']:r for r in rows if r['label']=='revised_reference'}
        if not all(r['status']=='valid' for r in [*group.values(),*reference.values()]):continue
        d120=group[120]['metrics']['peak_stage_m']-reference[120]['metrics']['peak_stage_m']
        d60=group[60]['metrics']['peak_stage_m']-reference[60]['metrics']['peak_stage_m']
        delta=sum(abs(g[120]['metrics']['peak_stage_m']-g[60]['metrics']['peak_stage_m']) for g in [group,reference])
        feasible_stable=all(g[120]['feasible']==g[60]['feasible'] for g in [group,reference])
        identical=group[60]['pair']==reference[60]['pair'] and np.allclose(group[60]['design'],reference[60]['design'],rtol=0,atol=3e-15)
        resolved=identical or (np.sign(d120)==np.sign(d60) and delta<max(.1*abs(d60),1e-10) and feasible_stable)
        if label=='joint_best' and (not resolved or not all(r['feasible'] for r in group.values())):need30=True
        comparisons.append(dict(first='revised_reference',second=label,difference_120s_m=d120,difference_60s_m=d60,sum_grid_changes_m=delta,resolved=bool(resolved),identical_design=identical))
    if any(r['status']=='invalid' for r in rows if r['label'] in ('revised_reference','joint_best')):need30=True
    if need30:
        for label,r in designs[:2]:evaluate(label,r['theta'],r['pair'],30)
    if need30:
        final={r['label']:r for r in rows if r['dt_s']==30}
        coarse={r['label']:r for r in rows if r['dt_s']==60}
        if all(k in final and final[k]['status']=='valid' and coarse[k]['status']=='valid' for k in ['revised_reference','joint_best']):
            d30=final['joint_best']['metrics']['peak_stage_m']-final['revised_reference']['metrics']['peak_stage_m']
            d60=coarse['joint_best']['metrics']['peak_stage_m']-coarse['revised_reference']['metrics']['peak_stage_m']
            delta=sum(abs(final[k]['metrics']['peak_stage_m']-coarse[k]['metrics']['peak_stage_m']) for k in ['revised_reference','joint_best'])
            comparisons.append(dict(first='revised_reference',second='joint_best',refinement='60_vs_30',difference_60s_m=d60,difference_30s_m=d30,sum_grid_changes_m=delta,resolved=bool(np.sign(d30)==np.sign(d60) and delta<max(.1*abs(d30),1e-10) and all(final[k]['feasible'] and coarse[k]['feasible'] for k in ['revised_reference','joint_best']))))
    save('fixed_design_timesteps.json',dict(records=rows,comparisons=comparisons,thirty_second_required=need30))
    blocks=[];controlled=[];controlled_gate_keys=[];skipped_duplicate_blocks=[]
    for label,r in designs[:2]:
        gate_key=r['theta']
        if any(np.allclose(gate_key,old,rtol=0,atol=3e-15) for old in controlled_gate_keys):
            skipped_duplicate_blocks.append(dict(block=label,identical_to='revised_reference',reason='All gates and covers match within 3e-15 coordinate-roundtrip tolerance; same four class comparisons reused.'));continue
        controlled_gate_keys.append(gate_key.copy())
        blocks.append(dict(id=label,description='Identical prescribed gate schedules and cover fractions within this block; only response-class pair varies.',gate_design=r['theta'][:-2],covers=r['theta'][-2:]))
        for pair in itertools.product(*config['categories']['admissible']):
            slug='controlled_'+label+'_'+'_'.join(pair)
            row=dict(block=label,label=slug,pair=list(pair),covers=r['theta'][-2:],theta=r['theta'],dt_s=60)
            try:
                case=model.simulate(r['theta'],pair,config,refine=2,forcing=forcing)
                filename=slug+'.npz';save_trajectory(OUT/filename,case)
                row.update(status='valid',feasible=model.is_feasible(case,config),metrics=case['metrics'],audit=case['audit'],trajectory_file=filename)
            except (historical.TrajectoryFailure,FloatingPointError) as exc:
                row.update(status='invalid',feasible=False,failure=str(exc),diagnostic=getattr(exc,'diagnostic',None))
            controlled.append(row);print('Controlled',label,pair,row['status'],row['feasible'],flush=True)
    save('controlled_comparisons.json',dict(records=controlled,blocks=blocks,dt_s=60,skipped_duplicate_blocks=skipped_duplicate_blocks))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--screen',action='store_true');p.add_argument('--run',action='store_true');p.add_argument('--compare',action='store_true');p.add_argument('--report',action='store_true');a=p.parse_args()
    config=historical.load_config(CONFIG);OUT.mkdir(parents=True,exist_ok=True)
    if a.screen:screen(config)
    if a.run:run_batch(config)
    if a.compare:compare(config)
    if a.report:
        from B.no_spillway_reporting import export_study
        export_study(config,OUT)

if __name__=='__main__':main()
