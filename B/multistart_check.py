"""Small fixed-problem initial-guess check; separate artifacts, unchanged official OPT.

python3 -B -m B.multistart_check --screen
python3 -B -m B.multistart_check --run
"""
from pathlib import Path
import argparse
from collections import Counter
import copy,csv,itertools,json,time
import numpy as np
from B import optimize_wetland as w
from B.official_opt import run_opt,load_official

OUT=w.ROOT/'results/multistart'
P=.04  # One fixed numerical coordinate map for every gate, category and start.
CALLBACKS=100
BATCH_SECONDS=1500.
BATCH_ODE_LIMIT=1250

def to_design(q):
    q=np.asarray(q,dtype=float);theta=q.copy()
    theta[:-2]=P*q[:-2]/((1-P)*(1-q[:-2])+P*q[:-2])
    return theta

def to_coordinates(theta):
    theta=np.asarray(theta,dtype=float);q=theta.copy()
    q[:-2]=(1-P)*theta[:-2]/(P*(1-theta[:-2])+(1-P)*theta[:-2])
    return q

def load_saved(name):
    data=np.load(w.ROOT/'results'/f'{name}.npz')
    meta=json.loads(str(data['metadata_json']))
    return np.array(data['theta']),meta['pair'],meta

def save_trajectory(path,run):
    np.savez_compressed(path,**{k:run[k] for k in ['time','x','u','y','theta','residuals','balance']},
        metadata_json=np.asarray(json.dumps(dict(pair=run['pair'],metrics=run['metrics'],audit=run['audit']),default=w.serial)))

def candidate_pool():
    ref,_,_=load_saved('reference');old,_,_=load_saved('best_feasible')
    first=ref.copy();first[:12]*=np.array([.8,.8,1.25,1.25,1.05,.8,.85,1.2,1.2,1.15,1.05,1.1])
    first[12:24]*=np.array([.75,.75,1.1,1.1,.9,.7,.8,1.25,1.25,1.1,1.05,1.1]);first[-2:]=[.30,.80]
    second=np.r_[[.03,.03,.015,.02,.045,.03,.04,.03,.025,.055,.045,.04],
                 [.04,.04,.025,.025,.05,.035,.042,.035,.035,.065,.055,.05],[.80,.35]]
    fallback=ref.copy();fallback[:12]=[.025,.026,.021,.019,.052,.028,.035,.027,.03,.058,.042,.04]
    fallback[12:24]=[.034,.033,.028,.026,.058,.032,.043,.036,.04,.065,.05,.048];fallback[-2:]=[.45,.70]
    return [('reference',ref),('early_retention_cover_30_80',first),('asymmetric_cover_80_35',second),
            ('previous_AA_best_design',old),('timing_variant_cover_45_70',fallback)]

def screen(config):
    forcing=w.fixed_forcing(config,w.base_parameters(config));records=[];accepted=[]
    for pair in itertools.product(*config['categories']['admissible']):
        found=[]
        for label,theta in candidate_pool():
            q=to_coordinates(theta);effective=to_design(np.clip(q,.1,.9))
            if not np.allclose(effective,theta,rtol=0,atol=2e-15):
                raise RuntimeError('Proposed fixed mapping still shifts requested start: '+label)
            row=dict(category_pair=list(pair),start_id=label,requested_design=theta,
                solver_coordinates=q,expected_first_design=effective)
            try:
                result=w.simulate(theta,pair,config,forcing=forcing)
                ok=w.is_feasible(result,config)
                row.update(feasible=ok,status='feasible' if ok else 'valid_constraint_violation',metrics=result['metrics'],audit=result['audit'])
                if ok:
                    save_trajectory(OUT/('_'.join(pair)+'_'+label+'_screen.npz'),result)
                    found.append(dict(category_pair=list(pair),start_id=label,theta=theta))
            except (w.TrajectoryFailure,FloatingPointError) as exc:
                row.update(feasible=False,status='invalid',failure=str(exc),diagnostic=getattr(exc,'diagnostic',None))
            records.append(row);print('Screen:',pair,label,row['status'],flush=True)
            if len(found)==3:break
        accepted.extend(found)
    w.save_json(OUT/'start_screening.json',dict(maximum_candidates_per_pair=5,records=records,accepted=accepted,
        selection='First three fully checked feasible candidates from declared deterministic pool; no assumed convex feasibility.',
        map=dict(p=P,formula='g=p*q/((1-p)*(1-q)+p*q); covers unchanged; q and all physical controls retain full [0,1] bounds')))
    return accepted

def run_batch(config):
    if (OUT/'batch_status.json').exists():raise RuntimeError('Batch already exists; preserving results. No automatic overwrite.')
    record=json.loads((OUT/'start_screening.json').read_text());starts=record['accepted']
    settings={**config['optimizer'],'max_evals':CALLBACKS}
    w.save_json(OUT/'batch_plan.json',dict(configuration_sha256=w.config_hash(config),configuration=config,
        settings=settings,requested_runs=len(starts),maximum_callbacks=CALLBACKS*len(starts),
        maximum_total_ode_attempts=BATCH_ODE_LIMIT,maximum_seconds=BATCH_SECONDS,coordinate_map_p=P,
        full_physical_bounds=[0,1],full_wrapper_coordinate_bounds=[0,1],native_official_bounds=[-1,1],
        constraint_and_objective_scaling_changed=False,reference='Unchanged previous B.2 replacement reference A/A',
        official_source={k:v for k,v in load_official().items() if k not in ['nms','ode4u','opt_hyp']}))
    forcing=w.fixed_forcing(config,w.base_parameters(config));began=time.monotonic();total_odes=len(record['records'])+len(json.loads((OUT/'affine_initialization_probe.json').read_text())['probes'])
    summaries=[];global_best=None;abort=False
    for start in starts:
        theta=np.asarray(start['theta']);pair=start['category_pair'];label=start['start_id'];slug='_'.join(pair)+'_'+label
        budget=dict(evaluations=0,start=time.monotonic());ev=w.Evaluator(config,pair,forcing,budget)
        first=[];result={};before=time.monotonic()
        def callback(q,constants=None):
            if time.monotonic()-began>BATCH_SECONDS or total_odes+budget['evaluations']>=BATCH_ODE_LIMIT-9:
                raise w.RunBudgetExceeded('Declared small-batch budget exhausted.')
            design=to_design(q)
            if not first:
                first.extend(design.tolist())
                if not np.allclose(design,theta,rtol=0,atol=2e-15):raise RuntimeError('Actual OPT first design shifted from requested feasible start.')
            return ev(design,constants)
        try:
            result=run_opt(callback,to_coordinates(theta),np.zeros(26),np.ones(26),settings)
            termination=result['status']
        except w.RunBudgetExceeded as exc:termination='batch_budget: '+str(exc);abort=True
        total_odes+=budget['evaluations']
        best=ev.best
        if best is not None:
            save_trajectory(OUT/(slug+'_best.npz'),best)
            if global_best is None or best['metrics']['peak_stage_m']<global_best['metrics']['peak_stage_m']:global_best=best
        reasons=Counter(r['failure'].split(':',1)[0] for r in ev.records if r['failure'])
        report=dict(category_pair=pair,start_id=label,requested_initial_design=theta,first_optimizer_design=first,
            first_design_feasible=bool(ev.records and ev.records[0]['feasible']),
            optimizer_coordinate_settings=result,physical_callback_log=ev.records,
            termination=termination,callbacks=len(ev.records),cached_callbacks=sum(r['cached'] for r in ev.records),
            actual_ode_attempts=budget['evaluations'],seconds=time.monotonic()-before,rejection_counts=dict(reasons),
            best_design=None if best is None else best['theta'],best_metrics=None if best is None else best['metrics'],
            best_audit=None if best is None else best['audit'])
        w.save_json(OUT/(slug+'_run.json'),report)
        summaries.append({k:v for k,v in report.items() if k not in ['optimizer_coordinate_settings','physical_callback_log','best_audit']})
        print('OPT:',pair,label,termination,'peak',None if best is None else best['metrics']['peak_stage_m'],
              'callbacks',len(ev.records),'rejected',sum(reasons.values()),flush=True)
        if abort:break
    if global_best is not None:save_trajectory(OUT/'batch_best.npz',global_best)
    w.save_json(OUT/'batch_status.json',dict(status='bounded_batch_finished' if not abort else 'bounded_batch_stopped',
        runs=summaries,seconds=time.monotonic()-began,ode_attempts_including_screening_and_probe=total_odes,
        global_best_pair=None if global_best is None else global_best['pair'],
        global_best_design=None if global_best is None else global_best['theta']))

def timestep_check(config):
    forcing=w.fixed_forcing(config,w.base_parameters(config));designs=[]
    for label,name in [('fixed_reference','reference'),('previous_best','best_feasible')]:
        theta,pair,_=load_saved(name);designs.append((label,theta,pair))
    p=OUT/'batch_best.npz'
    if p.exists():
        d=np.load(p);m=json.loads(str(d['metadata_json']));theta=d['theta'];pair=m['pair']
        if not any(pair==p2 and np.array_equal(theta,t2) for _,t2,p2 in designs):designs.append(('new_batch_best',theta,pair))
    rows=[]
    def evaluate(label,theta,pair,factor):
        dt=120//factor
        try:
            r=w.simulate(theta,pair,config,refine=factor,forcing=forcing)
            save_trajectory(OUT/f'{label}_{dt}s.npz',r)
            row=dict(label=label,dt_s=dt,design=theta,pair=pair,status='valid',feasible=w.is_feasible(r,config),metrics=r['metrics'],audit=r['audit'],
                rain_depth_m=float(np.trapezoid(r['u'][w.a2.U['P']],r['time'])),solar_energy_Wh_m2=float(np.trapezoid(r['u'][w.a2.U['R_s']],r['time'])/3600))
        except (w.TrajectoryFailure,FloatingPointError) as exc:row=dict(label=label,dt_s=dt,design=theta,pair=pair,status='invalid',feasible=False,failure=str(exc))
        rows.append(row);print('Resolution:',label,dt,row['status'],flush=True)
    for label,theta,pair in designs:
        for factor in [1,2]:evaluate(label,theta,pair,factor)
    by={(r['label'],r['dt_s']):r for r in rows};unresolved=not all(r['feasible'] for r in rows)
    comparisons=[]
    for i,(label,_,_) in enumerate(designs):
        if by[label,120]['status']=='valid' and by[label,60]['status']=='valid':
            change=abs(by[label,120]['metrics']['peak_stage_m']-by[label,60]['metrics']['peak_stage_m'])
            if change>config['numerics']['peak_grid_tolerance_m']:unresolved=True
        for label2,_,_ in designs[:i]:
            if not all(by[l,dt]['status']=='valid' for l in [label,label2] for dt in [120,60]):continue
            diff120=by[label,120]['metrics']['peak_stage_m']-by[label2,120]['metrics']['peak_stage_m']
            diff60=by[label,60]['metrics']['peak_stage_m']-by[label2,60]['metrics']['peak_stage_m']
            uncertainty=sum(abs(by[l,120]['metrics']['peak_stage_m']-by[l,60]['metrics']['peak_stage_m']) for l in [label,label2])
            resolved=np.sign(diff120)==np.sign(diff60) and uncertainty<.1*abs(diff60)
            comparisons.append(dict(first=label,second=label2,peak_difference_120_s_m=diff120,peak_difference_60_s_m=diff60,
                sum_grid_changes_m=uncertainty,resolved_by_declared_screen=bool(resolved)))
            unresolved|=not resolved
    if unresolved:
        for label,theta,pair in designs:evaluate(label,theta,pair,4)
    w.save_json(OUT/'fixed_design_timesteps.json',dict(rows=rows,comparisons_120_60=comparisons,ran_30s=bool(unresolved),
        criterion='30 s only if feasibility fails, any peak change exceeds existing .001m check, ranking changes, or sum of grid changes is >=10% of a compared peak difference. This is a resolution screen, not a rigorous error bound.',
        forcing_reconstruction='Same original 120s knots and linear gate histories at every refinement; no weather reinterpretation or reoptimization.'))

def export_report(config):
    """Export the completed bounded batch; no simulations or optimizer calls."""
    from collections import Counter
    import csv
    import itertools
    import json
    from pathlib import Path
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from B import optimize_wetland as w

    root=w.ROOT; out=root/'results/multistart'
    read=lambda name: json.loads((out/name).read_text())
    status=read('batch_status.json'); plan=read('batch_plan.json')
    screening=read('start_screening.json'); scaling=read('scaling_audit.json')
    resolution=read('fixed_design_timesteps.json'); probe=read('affine_initialization_probe.json')
    reports=[read('_'.join(s['category_pair'])+'_'+s['start_id']+'_run.json') for s in status['runs']]
    if not reports: raise RuntimeError('No completed multistart run records to report.')
    grid={(r['label'],r['dt_s']):r for r in resolution['rows']}
    reference=grid.get(('fixed_reference',120))
    if reference is None or reference.get('status')!='valid':
        raise RuntimeError('A valid fixed-reference 120-second reevaluation is required for reporting.')
    reference_peak=reference['metrics']['peak_stage_m']
    categories=lambda pair: '/'.join(s.removeprefix('type_').upper() for s in pair)
    short={'reference':'Reference start','early_retention_cover_30_80':'Timing / covers 0.30, 0.80',
           'asymmetric_cover_80_35':'Asymmetric / covers 0.80, 0.35',
           'previous_AA_best_design':'Previous A/A solution','timing_variant_cover_45_70':'Timing / covers 0.45, 0.70'}
    screen_by={(tuple(r['category_pair']),r['start_id']):r for r in screening['records']}
    good=[r for r in reports if r['best_metrics'] is not None]
    best=min(good,key=lambda r:r['best_metrics']['peak_stage_m']) if good else None
    rows=[]
    for r in reports:
        initial=r['requested_initial_design']; actual=r['first_optimizer_design']; final=r['best_design']
        audit={a['label']:a for a in (r['best_audit'] or [])}
        callback=r['physical_callback_log']; metric=r['best_metrics']
        screened=screen_by[tuple(r['category_pair']),r['start_id']]
        row=dict(category_pair=categories(r['category_pair']),start_id=r['start_id'],
            requested_initial_design_json=json.dumps(initial),first_optimizer_design_json=json.dumps(actual),
            final_best_feasible_design_json=json.dumps(final),
            initial_cover1=initial[-2],initial_cover2=initial[-1],
            first_cover1=actual[-2] if actual else None,first_cover2=actual[-1] if actual else None,
            final_cover1=final[-2] if final else None,final_cover2=final[-1] if final else None,
            first_design_max_absolute_shift=float(np.max(np.abs(np.asarray(initial)-np.asarray(actual)))) if actual else None,
            first_design_feasible=r['first_design_feasible'],
            screened_initial_peak_120s_m=screened.get('metrics',{}).get('peak_stage_m'),
            first_optimizer_peak_120s_m=callback[0]['raw_peak_stage_m'] if callback else None,
            physical_peak_120s_m=metric['peak_stage_m'] if metric else None,
            fixed_reference_peak_120s_m=reference_peak,
            reduction_from_fixed_reference_m=reference_peak-metric['peak_stage_m'] if metric else None,
            feasible=bool(metric and all(a['pass_at_tolerance'] for a in audit.values())),
            max_scaled_violation=max(0.,max(a['scaled_residual'] for a in audit.values())) if audit else None,
            callbacks=r['callbacks'],cached_callbacks=r['cached_callbacks'],actual_ode_attempts=r['actual_ode_attempts'],
            completed_valid_ode_simulations=sum(not x['cached'] and x['raw_peak_stage_m'] is not None for x in callback),
            hydraulically_invalid_callbacks=sum(x['hydraulically_invalid'] for x in callback),
            hydraulically_invalid_ode_attempts=sum(not x['cached'] and x['hydraulically_invalid'] for x in callback),
            valid_but_constraint_infeasible_callbacks=sum(not x['failure'] and not x['feasible'] for x in callback),
            rejection_counts_json=json.dumps(r['rejection_counts']),termination=r['termination'],
            converged=bool(r['optimizer_coordinate_settings'].get('converged',False)),elapsed_seconds=r['seconds'],
            balance_max_m3=metric['balance_max_m3'] if metric else None,
            balance_relative=metric['balance_relative'] if metric else None)
        for label in ['H_d1_safe','H_d2_safe','H_w1_safe','H_w2_safe','gate1_head_opening','gate2_head_opening','soil_capacity']:
            unit='m3' if label=='soil_capacity' else 'm'
            row[label+'_margin_'+unit]=audit.get(label,{}).get('margin')
        rows.append(row)
    table=root/'tables/multistart_summary.csv'; table.parent.mkdir(parents=True,exist_ok=True)
    with table.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

    # These two figures contain physical peaks only; invalid objective sentinels never enter either plot.
    pairs=list(itertools.product(*config['categories']['admissible']))
    colors=['#20639b','#d76a21','#39805a','#8e5b9e','#666666']
    fig,axs=plt.subplots(2,2,figsize=(10.5,7.3),sharey=True,constrained_layout=True)
    for ax,pair in zip(axs.flat,pairs):
        group=[r for r in reports if tuple(r['category_pair'])==pair]
        for i,r in enumerate(group):
            peak=None if r['best_metrics'] is None else r['best_metrics']['peak_stage_m']
            if peak is not None: ax.scatter(i,peak,s=72,color=colors[i],zorder=3)
        ax.axhline(reference_peak,color='black',linestyle='--',linewidth=1,label='Fixed A/A reference')
        ax.set_title('Synthetic classes '+categories(pair));ax.set_xticks(range(len(group)))
        ax.set_xticklabels([short.get(r['start_id'],r['start_id']).replace(' / ','\n') for r in group],fontsize=8)
        ax.ticklabel_format(axis='y',style='plain',useOffset=False);ax.grid(axis='y',alpha=.25)
        ax.set_ylabel('Best feasible physical river peak (m)')
    if good:
        low=min(r['best_metrics']['peak_stage_m'] for r in good);high=max(reference_peak,max(r['best_metrics']['peak_stage_m'] for r in good))
        pad=max((high-low)*.18,1e-5)
        for ax in axs.flat: ax.set_ylim(low-pad,high+pad)
    axs.flat[0].legend(fontsize=8);fig.suptitle('Bounded initial-guess comparison · 120-second integration',fontsize=13)
    plots=root/'plots';plots.mkdir(parents=True,exist_ok=True)
    fig.savefig(plots/'B2_multistart_objectives.png',dpi=170);plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(10.5,7.3),sharex=True,sharey=True,constrained_layout=True)
    for ax,pair in zip(axs.flat,pairs):
        group=[r for r in reports if tuple(r['category_pair'])==pair]
        for i,r in enumerate(group):
            count=0;xs=[];ys=[]
            for item in r['physical_callback_log']:
                if item['cached']:continue
                count+=1
                if item['best_feasible_peak_m'] is not None:xs.append(count);ys.append(item['best_feasible_peak_m'])
            ax.step(xs,ys,where='post',color=colors[i],label=short.get(r['start_id'],r['start_id']),linewidth=1.6)
        ax.axhline(reference_peak,color='black',linestyle='--',linewidth=.8)
        ax.set_title('Synthetic classes '+categories(pair));ax.grid(alpha=.2)
        ax.set_xlabel('Actual ODE attempts (cache hits excluded)');ax.set_ylabel('Best feasible physical river peak (m)')
        ax.ticklabel_format(axis='y',style='plain',useOffset=False);ax.legend(fontsize=7,loc='best')
    fig.suptitle('Feasible-best progress · rejected ODE attempts count on the horizontal axis',fontsize=12)
    fig.savefig(plots/'B2_multistart_convergence.png',dpi=170);plt.close(fig)

    fmt=lambda value,digits=6: '—' if value is None else f'{value:.{digits}f}'
    total_callbacks=sum(r['callbacks'] for r in reports);total_cache=sum(r['cached_callbacks'] for r in reports)
    total_odes=sum(r['actual_ode_attempts'] for r in reports)
    rejection=Counter();terminations=Counter(r['termination'] for r in reports)
    for r in reports:rejection.update(r['rejection_counts'])
    invalid_hyd=sum(x['hydraulically_invalid'] for r in reports for x in r['physical_callback_log'])
    valid_infeasible=sum(not x['failure'] and not x['feasible'] for r in reports for x in r['physical_callback_log'])
    max_shift=max(x['first_design_max_absolute_shift'] or 0 for x in rows)
    lines=['# Controlled B.2 initial-guess and scaling check','',
        'This bounded study retains the approved hypothetical engineered system, eight-state equations, gate hydraulics, G1–G6 assumptions, B.1 objective and constraints, forcing, initial storages, physical bounds and control-knot times. A/B remain synthetic hydraulic/ET response classes. The reference remains the previously accepted A/A replacement schedule at 0.1 times the original gate openings; it was not retuned.','',
        f"Completed {len(reports)} OPT runs after {len(screening['records'])} full-event starting-candidate checks; {sum(not r['feasible'] for r in screening['records'])} screened candidates were rejected. Each category received the feasible starts recorded in the [summary CSV](../tables/multistart_summary.csv).",'']
    if best:
        lines += [f"The lowest feasible 120-second physical peak in this batch is **{best['best_metrics']['peak_stage_m']:.9f} m**, from {categories(best['category_pair'])}, `{best['start_id']}`. Its reduction from the unchanged fixed reference ({reference_peak:.9f} m) is **{(reference_peak-best['best_metrics']['peak_stage_m'])*1000:.4f} mm**. These are physical river levels, separate from the penalized solver objective.",'']
    lines += ['| Synthetic classes | Start | Final covers | Peak (m, 120 s) | Reduction (mm) | ODE attempts | Termination |',
              '|---|---|---|---:|---:|---:|---|']
    for r in rows:
        lines.append(f"| {r['category_pair']} | {short.get(r['start_id'],r['start_id'])} | {fmt(r['final_cover1'],4)}, {fmt(r['final_cover2'],4)} | {fmt(r['physical_peak_120s_m'],9)} | {fmt(None if r['reduction_from_fixed_reference_m'] is None else 1000*r['reduction_from_fixed_reference_m'],4)} | {r['actual_ode_attempts']} | {r['termination']} |")
    lines += ['',f"The optimizer made {total_callbacks} callbacks: {total_cache} cache hits and {total_odes} actual ODE attempts. Separate work comprised {len(screening['records'])} screening attempts, {len(probe['probes'])} initialization-probe attempts and {len(resolution['rows'])} fixed-design timestep reevaluations. An ODE attempt can stop at its first hydraulic violation. Terminations: {dict(terminations)}. Rejection counts by callback: {dict(rejection)}; hydraulically invalid callbacks: {invalid_hyd}; valid trajectories that failed a constraint: {valid_infeasible}.",'',
        'The cap was declared before execution: 100 callbacks and 180 seconds per run, at most 1,200 optimization callbacks, plus the recorded batch caps of 1,250 ODE attempts and 25 minutes. All runs use the same tolerances (`tol_x=1e-4`, `tol_f=1e-6`, feasibility tolerance `1e-6`). Callback budgets are enforced before the next evaluation; an in-flight integration is allowed to finish. Official stagnation means no strict improvement of the augmented best for more than 20% of the evaluation budget. First-feasible stopping is disabled. Formal convergence requires more than 26²=676 iterations, so a 100-callback run cannot satisfy it. No run is evidence of a global optimum.','',
        '## Starting-point dependence','',
        '| Classes | Final peak range across starts (mm) | Largest gate-knot difference between final designs | Final cover-1 range | Final cover-2 range |',
        '|---|---:|---:|---|---|']
    for pair in pairs:
        group=[r for r in good if tuple(r['category_pair'])==pair]
        if not group:continue
        theta=np.asarray([r['best_design'] for r in group]);peaks=np.array([r['best_metrics']['peak_stage_m'] for r in group])
        distance=max((float(np.max(np.abs(a[:-2]-b[:-2]))) for a,b in itertools.combinations(theta,2)),default=0.)
        lines.append(f"| {categories(pair)} | {np.ptp(peaks)*1000:.4f} | {distance:.6f} | {theta[:,-2].min():.4f}–{theta[:,-2].max():.4f} | {theta[:,-1].min():.4f}–{theta[:,-1].max():.4f} |")
    near_pairs=[(a,b) for a,b in itertools.combinations(good,2) if a['category_pair']==b['category_pair']]
    if near_pairs:
        a,b=min(near_pairs,key=lambda pair:abs(pair[0]['best_metrics']['peak_stage_m']-pair[1]['best_metrics']['peak_stage_m']))
        ta=np.asarray(a['best_design']);tb=np.asarray(b['best_design'])
        lines += ['',f"For example, the closest same-category pair ({categories(a['category_pair'])}, {short[a['start_id']]} versus {short[b['start_id']]}) differs by only {1000*abs(a['best_metrics']['peak_stage_m']-b['best_metrics']['peak_stage_m']):.6f} mm in its 120-second peak, while its largest gate-knot difference is {np.max(abs(ta[:-2]-tb[:-2])):.6f} and its covers are ({ta[-2]:.6f}, {ta[-1]:.6f}) versus ({tb[-2]:.6f}, {tb[-1]:.6f}). Thus near-equal objectives conceal different designs. This near-tied pair is not itself a complete two-design refinement comparison; no precise ranking of all near-tied starts is claimed.",'']
    lines += ['', 'These ranges describe sensitivity of this short search to its starting design. Gate differences are normalized openings at the existing knots; covers are fractions. Similar peak values do not establish a unique operating schedule. Full requested, first-evaluated and final physical vectors, both covers and constraint margins are retained in the CSV. The final designs are the best feasible candidates observed, not necessarily the official solver’s penalized return point.','',
        '## Verified OPT interface and numerical remedy','',
        'The decision vector contains 12 independent gate-1 knot values, 12 independent gate-2 knot values and two independent covers. All 26 physical bounds remain [0,1]. The unchanged official adapter passes native bounds [−1,1] to the pinned `multivarious` NMS implementation. Its initial normalized point is clipped internally to [−0.8,0.8]. With the old affine full-box map, requested gate openings below 0.1 are therefore moved to 0.1; the two explicit probes recorded the resulting collapse. This is an initialization effect, not evidence that the requested starts were infeasible.','',
        'This batch uses one fixed gate map, `g = 0.04 q / (0.96 (1−q) + 0.04 q)`, with `q∈[0,1]`; covers remain affine. The map is a monotone bijection onto the full physical [0,1] interval and is identical for every start and category. All selected starts map inside the official initialization clip. It changes numerical search geometry, and its effect must not be attributed to a pure comparison with the previous affine search. It does not narrow physical bounds or modify the official solver. The official 27-vertex simplex perturbs all native coordinates by 0.03423609 and one coordinate by 0.24636812 for each added vertex. Individual gates and covers can all vary.','',
        f"Each actual first physical callback was checked against its requested design (maximum component difference {max_shift:.3g}); all were required to be feasible. The run JSON’s `optimizer_coordinate_settings` vectors/history are q-space; `physical_callback_log` and `best_design` contain physical openings/covers. The pinned official revision is `{plan['official_source']['revision']}`.",'',
        '## Constraint scaling and stopping resolution','',
        'Scaling is unchanged; no before/after scaling run was needed. Residual ≤0 means satisfied. Water levels use `(H−H_limit)/1 m`: a fixed depth scale comparable with the allowable 1.20 m pool and 0.90 m wetland depths, not an absolute-elevation divisor. Nonnegative storages use the eight fixed positive initial-storage scales. Soil capacity retains `(S_m−120000 m³)/90000 m³`; this is a valid fixed scale. Using capacity instead would reduce that residual by 0.75, change its physical tolerance from 0.09 to 0.12 m³, and alter its isolated squared penalty by a factor of 1/1.7778. That optional rescaling was deliberately not mixed with initial-guess changes. G5 uses the fixed 1 m length scale. Residuals are never clipped to [−1,1].','',
        'The unchanged official penalty is `F = f + 100 [sum(g_i for g_i > 1e-6)]²`, with `f=(physical H_peak−97 m)/1 m`. This is the square of a sum, and uses the full violating residual. Invalid hydraulic candidates receive finite numerical sentinels and diagnostics; no physical river peak is assigned to an incomplete trajectory. The relative objective-spread tolerance corresponds to about 0.0006 mm near these feasible peaks, which can distinguish 1.44 mm nominally. It is not an integration-error bound or a convergence certificate.','',
        'All returned reporting nodes, including the final time, enter the audit. Existing G5 and hard hydraulic-geometry checks also execute at every official RK stage. Ordinary safety/storage constraints and the additional 0.001 m geometry margin are audited at reporting nodes; they are not claimed to be checked at every intermediate RK state. Finite-time checks do not prove continuous-time feasibility. G5 is a strict runtime condition even though its reported row includes the nominal 1e-6 m audit tolerance.','']
    if best:
        expressions={r['label']:r['physical_residual_expression'] for r in scaling['constraint_records_from_previous_refined_runs']}
        lines += ['<details>','<summary>All 28 residual rows for the new batch best at 120 seconds</summary>','',
            'Physical margin = −physical residual. The first 20 rows are supplied to the penalty interface; the remaining eight rows independently audit native/reconstructed bounds.','',
            '| Constraint / physical expression | Units | Physical residual | Scale | Scaled residual | Physical tolerance | Physical margin |',
            '|---|---|---:|---:|---:|---:|---:|']
        for a in best['best_audit']:
            expression=expressions.get(a['label'],a['label']).replace('|','/').replace('<','&lt;')
            lines.append(f"| `{a['label']}`: {expression} | {a['units']} | {-a['margin']:.7g} | {a['scale']:.7g} | {a['scaled_residual']:.7g} | {a['physical_tolerance']:.7g} | {a['margin']:.7g} |")
        lines += ['','</details>','']
    lines += ['## Fixed-design timestep check','',
        'The fixed reference, previous best and new batch best were reevaluated without reoptimization, removing duplicate designs. Refinement interpolates the same original 120-second forcing knots and unchanged gate histories. Rainfall and solar energy are not reinterpreted at finer steps.','',
        '| Fixed design | Step (s) | Peak (m) | Reduction from same-step fixed reference (mm) | Max balance residual (m³) | Relative balance residual | Feasible |',
        '|---|---:|---:|---:|---:|---:|---|']
    for r in resolution['rows']:
        m=r.get('metrics',{});ref=grid.get(('fixed_reference',r['dt_s']),{}).get('metrics',{}).get('peak_stage_m')
        reduction=None if ref is None or 'peak_stage_m' not in m else (ref-m['peak_stage_m'])*1000
        lines.append(f"| {r['label']} | {r['dt_s']} | {fmt(m.get('peak_stage_m'),9)} | {fmt(reduction,4)} | {fmt(m.get('balance_max_m3'),6)} | {m.get('balance_relative','—')} | {r['feasible']} |")
    valid_grid=[r for r in resolution['rows'] if r['status']=='valid']
    if valid_grid:
        rain=np.array([r['rain_depth_m'] for r in valid_grid]);sun=np.array([r['solar_energy_Wh_m2'] for r in valid_grid])
        lines += ['',f"Across these runs, integrated rainfall is {rain[0]*1000:.6f} mm (maximum spread {np.ptp(rain):.3g} m) and solar energy is {sun[0]:.6f} Wh/m² (maximum spread {np.ptp(sun):.3g} Wh/m²). These checks preserve the prepared observation/grid/transfer assumptions.",'']
    labels=list(dict.fromkeys(r['label'] for r in resolution['rows']))
    finest=30 if resolution['ran_30s'] else 60;coarser=60 if finest==30 else 120
    lines += [f"The comparison below uses {coarser} versus {finest} seconds. Thirty seconds was {'triggered by the recorded 120/60-second screen' if resolution['ran_30s'] else 'not required by the recorded 120/60-second screen'}.",'',
        '| Difference (second design minus first) | Coarser difference (mm) | Finer difference (mm) | Sum of individual grid changes (mm) | Resolution screen |',
        '|---|---:|---:|---:|---|']
    for a,b in itertools.combinations(labels,2):
        rr=[grid.get((name,dt)) for name in (a,b) for dt in (coarser,finest)]
        if not all(r and r['status']=='valid' for r in rr):
            lines.append(f'| {b} − {a} | — | — | — | Unresolved: invalid/missing run |');continue
        ac,af,bc,bf=[r['metrics']['peak_stage_m'] for r in rr]
        dc=bc-ac;df=bf-af;change=abs(ac-af)+abs(bc-bf)
        resolved=all(r['feasible'] for r in rr) and np.sign(dc)==np.sign(df) and change<.1*abs(df)
        result='Stable at tested steps' if resolved else 'Magnitude/ranking not numerically resolved'
        lines.append(f'| {b} − {a} | {1000*dc:.6f} | {1000*df:.6f} | {1000*change:.6f} | {result} |')
    lines += ['', '“Stable at tested steps” requires unchanged ranking, feasible fixed-design runs and summed grid changes below 10% of the compared peak difference. This is an empirical resolution screen, not a rigorous error bound. Unresolved comparisons must not support precise millimeter-scale improvement claims. Multiple starts assess this bounded search only; they neither prove global optimality nor calibrate the regional hypothetical model.','',
        '![Physical objective by start](../plots/B2_multistart_objectives.png)','',
        '![Best-feasible objective by actual ODE attempts](../plots/B2_multistart_convergence.png)','',
        'Reproduce from the project root using the existing environment: `python3 -B -m B.multistart_check --screen`, then `python3 -B -m B.multistart_check --run`. The driver refuses to overwrite a completed batch. For the existing completed batch, `python3 -B -m B.multistart_check --report` regenerates only these reports/figures. Existing results and raw callback logs are local under `B/results/multistart/`; approved configuration, A.2 and previous B.2 outputs are preserved. No Git publication is part of this check.','']
    import re
    content='\n'.join(lines)
    # Keep the main explanation compact; complete per-run designs/results stay in the one CSV.
    content=re.sub(r'\| Synthetic classes \| Start \|.*?(?=The optimizer made)', '', content, flags=re.S)
    content=content.replace('These ranges describe sensitivity of this short search to its starting design. Gate differences are normalized openings at the existing knots; covers are fractions. Similar peak values do not establish a unique operating schedule. Full requested, first-evaluated and final physical vectors, both covers and constraint margins are retained in the CSV. The final designs are the best feasible candidates observed, not necessarily the official solver’s penalized return point.',
        'The CSV contains every requested, first-evaluated and final physical vector, both covers, physical margins and best feasible results. These results are separate from the solver’s penalized return point.')
    ref_runs=[r for r in reports if r['start_id']=='reference']
    reference_rejections=sum(sum(x['hydraulically_invalid'] for x in r['physical_callback_log']) for r in ref_runs)
    content=content.replace('## Starting-point dependence',
        f'The {len(ref_runs)} reference-start runs rejected {reference_rejections} of their {sum(r["callbacks"] for r in ref_runs)} callbacks for gate applicability; their only valid callback was the feasible first point. This explains their early stagnation. Other starts explored more valid candidates without changing the constraints.\n\n## Starting-point dependence')
    content=content.replace('The cap was declared before execution: 100 callbacks and 180 seconds per run, at most 1,200 optimization callbacks, plus the recorded batch caps of 1,250 ODE attempts and 25 minutes. All runs use the same tolerances (`tol_x=1e-4`, `tol_f=1e-6`, feasibility tolerance `1e-6`). Callback budgets are enforced before the next evaluation; an in-flight integration is allowed to finish. Official stagnation means no strict improvement of the augmented best for more than 20% of the evaluation budget. First-feasible stopping is disabled. Formal convergence requires more than 26²=676 iterations, so a 100-callback run cannot satisfy it. No run is evidence of a global optimum.',
        'Declared caps: 100 callbacks/180 s per run, 1,200 callbacks overall, 1,250 total ODE attempts and 25 minutes. Settings are identical: `tol_x=1e-4`, `tol_f=1e-6`, feasibility tolerance `1e-6`. Stagnation means no strict augmented-best improvement for over 20% of the evaluation budget; first-feasible stopping is disabled. Official convergence requires more than 676 iterations and is impossible under this small callback budget.')
    document=root/'docs/10_multistart_check.md';document.write_text(content)
    return dict(summary_csv=str(table),report=str(document),figures=[str(plots/'B2_multistart_objectives.png'),str(plots/'B2_multistart_convergence.png')])


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--screen',action='store_true');ap.add_argument('--run',action='store_true');ap.add_argument('--resolution',action='store_true');ap.add_argument('--report',action='store_true');args=ap.parse_args()
    OUT.mkdir(parents=True,exist_ok=True);config=w.load_config();assert not w.preflight(config)
    if args.screen:screen(config)
    if args.run:run_batch(config);timestep_check(config);export_report(config)
    elif args.resolution:timestep_check(config)
    if args.report:print(export_report(config))

if __name__=='__main__':main()
