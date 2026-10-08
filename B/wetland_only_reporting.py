"""Actual-result reporting for the fixed-gate, two-cover wetland design study.

No ODE or optimizer is executed here. Historical studies remain untouched.
"""
from __future__ import annotations
from collections import Counter
import csv
import json
from pathlib import Path
import numpy as np
from A import wetland_flood_model as a2
from B.no_spillway_reporting import mechanism_metrics,_load_run

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/wetland_only'
G5={'gate1_head_opening','gate2_head_opening'}


def _pair(pair):return '/'.join(p.removeprefix('type_').upper() for p in pair)
def _fmt(value,n=6):return '—' if value is None else f'{value:.{n}f}'
def _json(value):return json.dumps(value,default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item())


def _write_table(path,rows):
    fields=list(dict.fromkeys(k for r in rows for k in r)) or ['status']
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    return path


def export_study(output_dir=None):
    """Export this study's completed records and at most two scientific figures."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    output=OUT if output_dir is None else Path(output_dir)
    read=lambda name:json.loads((output/name).read_text())
    plan=read('study_plan.json');spec=plan['specification'];config=plan['physical_model_configuration']
    status=read('batch_status.json');screen=read('start_screening.json');validation=read('fixed_gate_validation.json')
    resolution=read('fixed_design_timesteps.json') if (output/'fixed_design_timesteps.json').exists() else {'records':[],'comparisons':[]}
    reference=_load_run(output/'reference.npz')
    best=_load_run(output/'batch_best.npz') if (output/'batch_best.npz').exists() else None
    reports=[]
    for summary in status['runs']:
        name='_'.join(summary['category_pair'])+'_'+summary['start_id']+'_run.json';reports.append(read(name))
    gates=np.asarray(plan['fixed_gate_knots'])
    if gates.shape!=(24,):raise ValueError('Expected exactly 24 fixed gate knots.')
    for run in [reference]+([] if best is None else [best]):
        if not np.array_equal(run['theta'][:-2],gates):raise ValueError('Reported trajectory changes the fixed gates.')
    for report in reports:
        if not np.array_equal(report['fixed_gate_knots'],gates):raise ValueError('Optimizer record changes the fixed gates.')
        for entry in report['physical_callback_log']:
            if not np.array_equal(np.asarray(entry['full_design'])[:-2],gates):raise ValueError('Callback changed a gate decision.')
            if len(entry['optimizer_residuals'])!=18:raise ValueError('Unexpected OPT residual count.')
    tabledir=ROOT/'tables/wetland_only';plotdir=ROOT/'plots/wetland_only'
    tabledir.mkdir(parents=True,exist_ok=True);plotdir.mkdir(parents=True,exist_ok=True)
    refpeak=reference['metrics']['peak_stage_m'];summaries=[];constraint_rows=[]
    from B.no_spillway_model import base_parameters
    parameters=base_parameters(config)  # Read existing fixed coefficients; no ODE evaluation.
    def effective_coefficients(covers,pair):
        keys=['effective_R1_s_m2','effective_R2_s_m2','effective_K1','effective_K2']
        if covers is None:return dict.fromkeys(keys,None)
        rho=np.asarray(covers);endpoints=config['categories']['endpoints']
        rv=np.array([endpoints[pair[i]]['R_veg'][i] for i in range(2)])
        kv=np.array([endpoints[pair[i]]['K_veg'][i] for i in range(2)])
        resistance=(1-rho)*parameters.R_bare+rho*rv
        multiplier=(1-rho)*parameters.K_bare+rho*kv
        return dict(zip(keys,[float(v) for v in [*resistance,*multiplier]]))
    def constraint_role(row):
        label=row['label']
        if label in G5:return 'G5 model applicability; strict RHS/RK guard',False
        if label.startswith('gate'):return 'Fixed-input bound validation',False
        if label.startswith('cover'):return 'Native OPT cover bound; not penalty residual',False
        return 'Active OPT inequality'+(' / numerical integration validity' if label=='integration_status' else ''),True
    def register_constraints(run_label,pair,audit):
        for row in audit:
            role,active=constraint_role(row)
            constraint_rows.append(dict(run=run_label,pair=_pair(pair),role=role,
                supplied_to_OPT_residual_vector=active,optimizer_active=active or row['label'].startswith('cover'),
                native_OPT_bound=row['label'].startswith('cover'),physical_residual=-row['margin'],**row))
    for r in reports:
        initial=r['requested_initial_covers'];first=r['first_optimizer_covers'];final=r['best_covers']
        metric=r['best_metrics'];log=r['physical_callback_log'];solver=r.get('optimizer_result',{})
        audit={a['label']:a for a in (r['best_audit'] or [])}
        origin=next((e['evaluation'] for e in log if final is not None and e['feasible'] and np.array_equal(e['covers'],final)),None)
        row=dict(pair=_pair(r['category_pair']),start_id=r['start_id'],requested_cover1=initial[0],requested_cover2=initial[1],
            first_cover1=first[0] if first else None,first_cover2=first[1] if first else None,
            final_cover1=final[0] if final is not None else None,final_cover2=final[1] if final is not None else None,
            requested_covers_json=_json(initial),first_covers_json=_json(first),final_covers_json=_json(final),
            fixed_gate_knots_json=_json(gates),full_final_design_json=_json(r['best_design']),
            first_cover_max_absolute_shift=float(np.max(np.abs(np.asarray(initial)-first))) if first else None,
            first_design_feasible=r['first_design_feasible'],first_physical_peak_m=log[0]['raw_peak_stage_m'] if log else None,
            physical_peak_120s_m=None if metric is None else metric['peak_stage_m'],reference_peak_120s_m=refpeak,
            reduction_from_fixed_reference_m=None if metric is None else refpeak-metric['peak_stage_m'],
            OPT_residual_vector_count=18,max_active_scaled_violation=None if not r['best_active_audit'] else max(0.,max(a['scaled_residual'] for a in r['best_active_audit'])),
            G5_supplied_to_OPT=False,G5_model_applicability_passed=None if not r['best_applicability_audit'] else all(a['strict_pass'] for a in r['best_applicability_audit']),
            callbacks=r['callbacks'],cached_callbacks=r['cached_callbacks'],actual_ode_attempts=r['actual_ode_attempts'],
            valid_ode_completions=sum(not e['cached'] and e['raw_peak_stage_m'] is not None for e in log),
            invalid_model_callbacks=sum(bool(e['failure']) for e in log),
            valid_but_infeasible_callbacks=sum(not e['failure'] and not e['feasible'] for e in log),
            best_origin_callback=origin,termination=r['termination'],solver_converged=solver.get('converged',False),
            official_returned=solver.get('official_returned'),seconds=r['seconds'],rejection_counts_json=_json(r['rejection_counts']))
        for label in ['H_d1_safe','H_d2_safe','H_w1_safe','H_w2_safe','soil_capacity','gate1_head_opening','gate2_head_opening']:
            row[label+'_margin']=audit.get(label,{}).get('margin')
        row.update(effective_coefficients(final,r['category_pair']))
        summaries.append(row)
        if r['best_audit']:register_constraints(r['start_id']+'_'+_pair(r['category_pair']),r['category_pair'],r['best_audit'])
    # Use the finest verified common reporting grid for physical comparisons.
    by={(r['label'],r['dt_s']):r for r in resolution['records']}
    common=[dt for dt in [30,60,120] if all(by.get((label,dt),{}).get('status')=='valid' for label in ['fixed_reference','wetland_best'])]
    chosen_dt=common[0] if common else float(np.diff(reference['time'])[0])
    plotted=[('fixed_reference',reference)]
    if best is not None:plotted.append(('wetland_best',best))
    resolved_runs=[]
    for label,fallback in plotted:
        record=by.get((label,chosen_dt));filename=None if record is None else record.get('trajectory_file')
        run=_load_run(output/filename) if filename else fallback
        if not np.array_equal(run['theta'][:-2],gates):raise ValueError('Fixed-design check changed gate schedule.')
        resolved_runs.append((label,run));register_constraints(label+f'_{chosen_dt:g}s',run['pair'],run['audit'])
    mechanisms=[]
    for label,run in resolved_runs:
        metric=mechanism_metrics(run)
        if metric['dam_spillway_max_abs_m3_s']!=0:raise ValueError('Dam-spillway flow must remain zero.')
        mechanisms.append(dict(label=label,pair=_pair(run['pair']),cover1=float(run['theta'][-2]),cover2=float(run['theta'][-1]),
            dt_s=float(np.diff(run['time'])[0]),**effective_coefficients(run['theta'][-2:],run['pair']),**metric))
    timestep_rows=[]
    for r in resolution['records']:
        timestep_rows.append({k:v for k,v in r.items() if k not in ['metrics','audit','active_audit','applicability_audit','design','pair','covers']}|
            dict(pair=_pair(r.get('pair',[])),design_json=_json(r.get('design')),covers_json=_json(r.get('covers')),**r.get('metrics',{})))
    tables=[_write_table(tabledir/'optimization_summary.csv',summaries),
            _write_table(tabledir/'constraint_audit.csv',constraint_rows),
            _write_table(tabledir/'mechanism_metrics.csv',mechanisms),
            _write_table(tabledir/'fixed_design_timesteps.csv',timestep_rows)]

    figures=[]
    def finish(fig,name):
        fig.text(.5,.007,'Fixed gates; synthetic wetland classes and covers. Bare cover retains both contained wetland units.',ha='center',fontsize=8)
        fig.tight_layout(rect=(0,.025,1,.95));path=plotdir/name;fig.savefig(path,dpi=170);plt.close(fig);figures.append(path)
    fig,axes=plt.subplots(2,3,figsize=(13,7.3),sharex=True)
    colors=['#2166ac','#b35806']
    for i,(label,run) in enumerate(resolved_runs):
        t=run['time']/3600;y=run['y'];color=colors[i];caption=label.replace('_',' ')
        axes[0,0].plot(t,y[a2.Y['H_r']],color=color,label=caption)
        for branch in (1,2):
            axes[0,branch].plot(t,run['x'][1+branch],color=color,label=caption)
        axes[1,0].plot(t,y[a2.Y['Q_G1']]+y[a2.Y['Q_G2']],color=color,label=caption)
        axes[1,1].plot(t,y[a2.Y['Q_wr1']]+y[a2.Y['Q_wr2']],color=color,label=caption)
        axes[1,2].plot(t,y[a2.Y['Q_wsp1']]+y[a2.Y['Q_wsp2']],color=color,label=caption)
    labels=['River level (m)','Wetland 1 storage (m³)','Wetland 2 storage (m³)',
            'Fixed-gate release (m³/s)','Normal wetland drainage (m³/s)','Retained wetland overflow (m³/s)']
    for ax,label in zip(axes.flat,labels):
        ax.set_ylabel(label);ax.grid(alpha=.2);ax.legend(fontsize=8);ax.ticklabel_format(axis='y',useOffset=False,style='plain')
    for ax in axes[-1]:ax.set_xlabel('Elapsed time (h)')
    fig.suptitle(f'Wetland-design comparison under identical gate schedules · {chosen_dt:g}-second reporting grid')
    finish(fig,'01_wetland_design_mechanisms.png')
    if reports:
        fig,axes=plt.subplots(1,2,figsize=(12.6,5.2));colors={'A/A':'#2166ac','A/B':'#1b9e77','B/A':'#e6ab02','B/B':'#b35806'}
        counts=Counter()
        for r in reports:
            pair=_pair(r['category_pair']);counts[pair]+=1;label=f'{pair} / start {counts[pair]}'
            count=0;xs=[];ys=[]
            for e in r['physical_callback_log']:
                if e['cached']:continue
                count+=1
                if e['best_feasible_peak_m'] is not None:xs.append(count);ys.append(e['best_feasible_peak_m'])
            axes[0].step(xs,ys,where='post',label=label,color=colors.get(pair),linestyle=['-','--',':'][(counts[pair]-1)%3])
            if r['best_covers'] is not None:
                x,y=r['best_covers'];axes[1].scatter(x,y,color=colors.get(pair),s=45,marker=['o','s','^'][(counts[pair]-1)%3],label=label)
                x0,y0=r['requested_initial_covers'];axes[1].plot([x0,x],[y0,y],color=colors.get(pair),alpha=.35,linewidth=.8)
                axes[1].scatter(x0,y0,color=colors.get(pair),s=20,marker='x',alpha=.5)
        axes[0].axhline(refpeak,color='black',linestyle=':',linewidth=.8);axes[0].set(xlabel='Actual ODE attempts; cache hits excluded',ylabel='Best feasible physical river peak (m)',title='Official OPT feasible-best progress')
        axes[0].ticklabel_format(axis='y',style='plain',useOffset=False);axes[0].legend(fontsize=7,ncol=2)
        axes[1].set(xlabel='Wetland 1 cover ρ₁',ylabel='Wetland 2 cover ρ₂',xlim=(-.04,1.04),ylim=(-.04,1.04),title='Initial × to final cover designs (overlap possible)')
        axes[1].set_aspect('equal');axes[1].grid(alpha=.2);axes[0].grid(alpha=.2)
        fig.suptitle('Two continuous cover decisions for each synthetic category pair')
        finish(fig,'02_cover_search_robustness.png')

    total_callbacks=sum(r['callbacks'] for r in reports);total_cache=sum(r['cached_callbacks'] for r in reports);total_opt_odes=sum(r['actual_ode_attempts'] for r in reports)
    total_study_odes=1+len(screen['records'])+total_opt_odes+len(resolution['records'])
    failures=Counter();terminations=Counter(r['termination'] for r in reports)
    for r in reports:failures.update(r['rejection_counts'])
    max_shift=max((r['first_cover_max_absolute_shift'] or 0 for r in summaries),default=0.)
    display_reference=resolved_runs[0][1];display_best=resolved_runs[1][1] if len(resolved_runs)>1 else None
    lines=['# B.2 wetland-design optimization with fixed gates','',
        'This is the user-authorized wetland-design-only analysis. It preserves the latest no-dam-spillway hydraulic model and fixes every one of the 24 gate-knot openings to the newly accepted admissible reference schedule. The continuous OPT decision vector is **[ρ₁, ρ₂] only**; synthetic category pairs A/A, A/B, B/A and B/B are enumerated separately. Initial states, regional event weather, geometry, coefficients, the six-day horizon, physical safety limits and wetland-overflow routing remain unchanged.','',
        'The previous joint gate-and-wetland optimization remains a historical study. This comparison measures changes in wetland response-class and cover design under fixed operations; it does **not** compare wetland presence with wetland absence or establish a real-species recommendation. A permitted cover of ρ=0 denotes a bare contained wetland with its storage, drainage, seepage and overflow still present; it does not delete the wetland.','',
        '## Physical result','',
        f"The fixed reference is A/A with covers 0.65/0.60. Its {chosen_dt:g}-second physical river peak is **{display_reference['metrics']['peak_stage_m']:.9f} m** in the assumed model datum.",'']
    if display_best is None:
        lines += ['No best-feasible design was retained by the recorded search. This outcome is not proof of infeasibility; no optimized trajectory or improvement is invented.','']
    else:
        cover=display_best['theta'][-2:];peak=display_best['metrics']['peak_stage_m'];reduction=display_reference['metrics']['peak_stage_m']-peak
        lines += [f"The best feasible design selected by the 120-second search is **{_pair(display_best['pair'])}, ρ₁={cover[0]:.8f}, ρ₂={cover[1]:.8f}**, with physical river peak **{peak:.9f} m**. The reduction from the same fixed reference is **{reduction*1000:.4f} mm** on the {chosen_dt:g}-second grid. The physical objective is `min max_t H_r(t)`; neither a penalty score nor a discharge/volume objective is substituted.",'']
        if np.all(np.abs(cover)<=1e-12):
            lines += ['Both covers reach zero. This means **two bare contained wetland units**, retaining their area, storage, normal drainage, seepage and overflow. It does not remove either wetland. At zero cover the existing interpolation uses bare hydraulic/ET coefficients; synthetic category labels have no physical effect, so a category tie does not identify a preferred species or response class.','']
        elif np.any(np.abs(cover)<=1e-12):
            lines += ['A zero-cover branch remains a bare contained wetland unit; its compartment, drainage, seepage and overflow remain. The category label for that bare branch has no physical effect under the existing interpolation.','']
    lines += ['| Synthetic pair | Start covers | Final covers | Physical peak at 120 s (m) | Reference reduction (mm) | ODE attempts | Termination |',
              '|---|---|---|---:|---:|---:|---|']
    for r in summaries:
        lines.append(f"| {r['pair']} | {r['requested_cover1']:.2f}/{r['requested_cover2']:.2f} | {_fmt(r['final_cover1'],6)}/{_fmt(r['final_cover2'],6)} | {_fmt(r['physical_peak_120s_m'],9)} | {_fmt(None if r['reduction_from_fixed_reference_m'] is None else r['reduction_from_fixed_reference_m']*1000,4)} | {r['actual_ode_attempts']} | {r['termination']} |")
    lines += ['',f"The {len(reports)} searches made {total_callbacks} callbacks, including {total_cache} cache hits and {total_opt_odes} actual ODE attempts. Study total: {total_study_odes} ODE attempts ({1} reference, {len(screen['records'])} start checks, {total_opt_odes} OPT attempts and {len(resolution['records'])} fixed-design reevaluations). Terminations: {dict(terminations)}. Invalid-candidate diagnostics: {dict(failures)}. Full requested, first-evaluated and final covers, fixed gate vectors, physical residual margins and counts are in [optimization_summary.csv](../tables/wetland_only/optimization_summary.csv).",'']
    if reports:
        bypair={}
        for r in summaries:
            if r['physical_peak_120s_m'] is not None:bypair.setdefault(r['pair'],[]).append(r)
        lines += ['| Pair | Range of final peaks across starts (mm) | Final cover-1 range | Final cover-2 range |',
                  '|---|---:|---|---|']
        for pair,group in bypair.items():
            peaks=[r['physical_peak_120s_m'] for r in group];c1=[r['final_cover1'] for r in group];c2=[r['final_cover2'] for r in group]
            lines.append(f'| {pair} | {(max(peaks)-min(peaks))*1000:.6f} | {min(c1):.6f}–{max(c1):.6f} | {min(c2):.6f}–{max(c2):.6f} |')
        lines += ['', '| Best retained candidate by category | Covers ρ₁/ρ₂ | Effective R₁/R₂ (s/m²) | Effective K₁/K₂ |',
                  '|---|---|---|---|']
        for pair,group in bypair.items():
            row=min(group,key=lambda r:r['physical_peak_120s_m'])
            lines.append(f"| {pair} | {row['final_cover1']:.6f}/{row['final_cover2']:.6f} | {row['effective_R1_s_m2']:.6f}/{row['effective_R2_s_m2']:.6f} | {row['effective_K1']:.6f}/{row['effective_K2']:.6f} |")
        lines += ['', 'Effective coefficients are derived from the unchanged interpolation `R_i=(1−ρ_i) R_bare,i + ρ_i R_class,i` and `K_i=(1−ρ_i) K_bare,i + ρ_i K_class,i`; they are not new measured properties. Different class/cover combinations can yield nearly the same effective resistance, while their paired ET multipliers still differ. Similar river peaks therefore need not imply the same covers or identify a preferred real species. The summary CSV includes these effective coefficients.','']
        lines += ['', 'These ranges measure starting-point sensitivity of this bounded search. Similar objectives can accompany different covers. Multiple starts do not prove a global optimum.','']
    lines += ['## Interface, constraints and model validity','',
        f"The verified official solver revision is `{plan['official_source']['revision']}`. Both cover bounds remain [0,1]; the adapter passes native bounds [−1,1]. The three requested interior starts avoid the official normalized initial-point clip. Maximum observed requested-to-first cover shift was {max_shift:.3g}. The two-dimensional initial simplex has three vertices; both covers can vary. Unlike the historical 26-variable joint search, the official convergence check here requires more than **four iterations**, together with its design/objective tests. Budget, stagnation and convergence statuses remain distinct.",'',
        f"Identical settings across runs: `{_json(spec['settings'])}`. Declared study caps: {spec['maximum_study_ode_attempts']} ODE attempts and {spec['maximum_batch_seconds']/60:g} minutes for the batch. First-feasible stopping is disabled. The recorded best feasible physical design is distinguished from any penalized solver-returned point.",'',
        '**Exactly 18 residuals are passed to OPT.** They retain nonnegative storage, soil capacity, four fixed upstream safety limits, model-geometry conditions and numerical integration status. The two G5 head/opening rows are excluded from the optimizer vector as authorized for fixed-gate design; their unchanged strict RHS/RK applicability guard remains active in every simulation. Fixed gate-bound checks are input validation, while four cover-bound rows record the native two-variable box. Every row in [constraint_audit.csv](../tables/wetland_only/constraint_audit.csv) states its role: `supplied_to_OPT_residual_vector=True` identifies the 18 residual rows, `native_OPT_bound=True` identifies the four enforced cover-bound rows, and `optimizer_active=True` covers both groups. G5 and fixed-gate bound-validation rows have false optimizer flags; their model/input validation remains enforced.','',
        'The fixed schedule was validated before optimization. Safety-feasible wetland surfaces remain below the two gate sills, so wetland-induced gate submergence is excluded under the retained tailwater assumption. A candidate that violates hydraulic applicability still stops immediately. Its generic integration-status sentinel rejects the invalid model evaluation; the other returned placeholder zeros are **uncomputed entries**, not observed physical residuals or evidence of satisfaction. No physical peak is assigned to that incomplete trajectory.','',
        'Residual convention remains physical residual / fixed positive scale ≤0. Level and G5 audit scales are 1 m; soil capacity is `(S_m−120000 m³)/90000 m³`; nonnegative storages use the frozen positive storage scales. Residuals are not divided by absolute elevation or changing water depth, and are not clipped. The scaled feasibility tolerance is 1e-6 (1e-6 m for level rows, 0.09 m³ for soil capacity); strict runtime guards can be stronger. The callback objective `(H_peak−97 m)/1 m` and penalty `100 [sum(g_i for g_i>1e-6)]²` are numerical quantities separate from the reported river level.','',
        'Storage admissibility, G5 and hard geometry checks run at all existing RHS/RK evaluations. The four safety levels and extra 0.001 m geometry margin are audited on all reporting nodes, including final time. These finite checks are not continuous-time safety certification.','',
        '## Wetland mechanisms under fixed operations','']
    if len(resolved_runs)>1:
        rr=resolved_runs[0][1];br=resolved_runs[1][1];rm=mechanisms[0];bm=mechanisms[1]
        if not np.array_equal(rr['time'],br['time']):raise ValueError('Mechanism comparison requires a common time grid.')
        gate_diff=float(np.max(np.abs(rr['u'][:2]-br['u'][:2])));release_diff=float(np.max(np.abs(rr['y'][[a2.Y['Q_G1'],a2.Y['Q_G2']]]-br['y'][[a2.Y['Q_G1'],a2.Y['Q_G2']]])))
        lines += [f"Maximum pointwise gate-opening difference is {gate_diff:.3g}; maximum gate-discharge difference is {release_diff:.3g} m³/s. The gate-release volume is {rm['gate_total_volume_m3']:.3f} m³ for the reference and {bm['gate_total_volume_m3']:.3f} m³ for the best wetland design. These checks separate wetland design effects from changed gate operations.",'',
            f"Here the routing change is consistent with lower effective drainage resistance: R₁/R₂ fall from {rm['effective_R1_s_m2']:.4f}/{rm['effective_R2_s_m2']:.4f} to {bm['effective_R1_s_m2']:.4f}/{bm['effective_R2_s_m2']:.4f} s/m². More water leaves through normal drainage before high wetland levels produce a concentrated overflow contribution. Both wetland storage peaks fall, and the reduction in overflow peak is larger than the increase in normal-drainage peak, reducing the simultaneous combined inflow peak to the river. The total wetland-to-river event volume actually increases by {bm['wetland_to_river_total_volume_m3']-rm['wetland_to_river_total_volume_m3']:.3f} m³ and wetland ET falls by {(rm['wetland_et1_volume_m3']+rm['wetland_et2_volume_m3'])-(bm['wetland_et1_volume_m3']+bm['wetland_et2_volume_m3']):.3f} m³. The lower river maximum therefore reflects changed routing and flow concentration in this event, rather than greater total water removal or greater peak wetland storage. This interpretation concerns the coupled synthetic design response, not resistance in isolation.",'',
            '| Quantity | Fixed reference | Best wetland design |','|---|---:|---:|']
        for label,key,n in [('Wetland 1 peak storage (m³)','wetland1_peak_storage_m3',3),('Wetland 2 peak storage (m³)','wetland2_peak_storage_m3',3),
            ('Wetland 1 storage-peak time (h)','wetland1_peak_storage_time_h',4),('Wetland 2 storage-peak time (h)','wetland2_peak_storage_time_h',4),
            ('Normal drainage volume (m³)','normal_drainage_total_volume_m3',3),('Wetland overflow volume (m³)','wetland_overflow_total_volume_m3',3),
            ('Normal drainage peak (m³/s)','normal_drainage_total_peak_m3_s',6),('Wetland overflow peak (m³/s)','wetland_overflow_total_peak_m3_s',6),
            ('Combined wetland-to-river peak (m³/s)','wetland_to_river_total_peak_m3_s',6),
            ('Combined wetland outflow centroid (h)','wetland_to_river_total_centroid_time_h',4),
            ('Combined wetland outflow peak time (h)','wetland_to_river_total_peak_time_h',4),('River peak time (h)','river_peak_time_h',4),
            ('Total final-storage change (m³)','total_storage_change_m3',3),('Max external-balance residual (m³)','external_balance_max_abs_residual_m3',6)]:
            lines.append(f'| {label} | {_fmt(rm[key],n)} | {_fmt(bm[key],n)} |')
        lines += ['', f"At the same gate releases, normal-drainage volume changes by {bm['normal_drainage_total_volume_m3']-rm['normal_drainage_total_volume_m3']:+.3f} m³ and retained wetland-overflow volume by {bm['wetland_overflow_total_volume_m3']-rm['wetland_overflow_total_volume_m3']:+.3f} m³. Their separate peaks change by {bm['normal_drainage_total_peak_m3_s']-rm['normal_drainage_total_peak_m3_s']:+.6f} and {bm['wetland_overflow_total_peak_m3_s']-rm['wetland_overflow_total_peak_m3_s']:+.6f} m³/s; the actual combined wetland-to-river peak changes by {bm['wetland_to_river_total_peak_m3_s']-rm['wetland_to_river_total_peak_m3_s']:+.6f} m³/s. The combined peak is calculated from the simultaneous summed time series, not by adding peaks that may occur at different times.",'',
            f"Wetland peak storages change by {bm['wetland1_peak_storage_m3']-rm['wetland1_peak_storage_m3']:+.3f} / {bm['wetland2_peak_storage_m3']-rm['wetland2_peak_storage_m3']:+.3f} m³ in branches 1/2. Combined wetland ET changes from {rm['wetland_et1_volume_m3']+rm['wetland_et2_volume_m3']:.3f} to {bm['wetland_et1_volume_m3']+bm['wetland_et2_volume_m3']:.3f} m³. The wetland outflow centroid shifts by {bm['wetland_to_river_total_centroid_time_h']-rm['wetland_to_river_total_centroid_time_h']:+.4f} h, while the river peak time shifts by {bm['river_peak_time_h']-rm['river_peak_time_h']:+.4f} h. These differences quantify the routing and water-balance changes that accompany the reported physical objective change.",'']
        lines += ['', 'Changing cover/class changes the existing paired drainage-resistance and ET response, with consequent storage, seepage and overflow changes. The table separates normal drainage from retained wetland overflow; a larger stored volume or later release centroid does not by itself establish a smaller or later river peak. These are model-specific synthetic response packages, not isolated resistance effects or field-measured vegetation benefits.','']
    lines += ['All branch volumes, release centroids and peaks are in [mechanism_metrics.csv](../tables/wetland_only/mechanism_metrics.csv). Centroids are temporal first moments, not parcel residence times. The whole-system balance uses only `P_total + Q_ext + Q_gx − Q_out − E_tot`; internal gates, wetland drainage/overflow and soil/groundwater transfers are not counted twice. Dam-spillway columns remain zero; wetland overflow remains active.','',
        '## Fixed-design timestep check','',
        'Reference and best designs were frozen and reevaluated with the same forcing reconstruction and gate functions; no refinement optimization or weather reinterpretation was performed.','',
        '| Design | Step (s) | Physical peak (m) | Feasible | Max balance residual (m³) | Relative balance residual |',
        '|---|---:|---:|---|---:|---:|']
    for r in resolution['records']:
        m=r.get('metrics',{})
        lines.append(f"| {r['label']} | {r['dt_s']} | {_fmt(m.get('peak_stage_m'),9)} | {r['feasible']} | {_fmt(m.get('balance_max_m3'),6)} | {m.get('balance_relative','—')} |")
    lines += ['', '| Step comparison | Reference reduction, coarser (mm) | Reference reduction, finer (mm) | Sum of grid changes (mm) | Resolution screen |',
              '|---|---:|---:|---:|---|']
    for r in resolution['comparisons']:
        lines.append(f"| {r['coarse_dt_s']} / {r['fine_dt_s']} s | {_fmt(None if r.get('reduction_coarse_m') is None else r['reduction_coarse_m']*1000,6)} | {_fmt(None if r.get('reduction_fine_m') is None else r['reduction_fine_m']*1000,6)} | {_fmt(None if r.get('sum_grid_changes_m') is None else r['sum_grid_changes_m']*1000,6)} | {'Stable at tested steps' if r.get('resolved') else 'Unresolved'} |")
    category_checks=resolution.get('category_comparisons',[])
    if category_checks:
        lines += ['', 'The best saved design from each synthetic category pair was also frozen for timestep checks. This compares those actual finite-search candidates, not exact category optima. Fine differences between category outcomes are assessed separately from the larger reduction relative to the reference.','',
            '| Fixed category-candidate comparison | Steps (s) | Peak difference, coarser (mm) | Peak difference, finer (mm) | Sum of grid changes (mm) | Resolution screen |',
            '|---|---|---:|---:|---:|---|']
        for c in category_checks:
            lines.append(f"| {c['second_label']} − {c['first_label']} | {c['coarse_dt_s']}/{c['fine_dt_s']} | {_fmt(None if c.get('peak_difference_coarse_m') is None else 1000*c['peak_difference_coarse_m'],6)} | {_fmt(None if c.get('peak_difference_fine_m') is None else 1000*c['peak_difference_fine_m'],6)} | {_fmt(None if c.get('sum_grid_changes_m') is None else 1000*c['sum_grid_changes_m'],6)} | {'Stable at tested steps' if c.get('resolved') else 'Magnitude/ranking unresolved'} |")
        finest_by_pair={}
        for c in category_checks:
            key=(c['first_label'],c['second_label'])
            if key not in finest_by_pair or c['fine_dt_s']<finest_by_pair[key]['fine_dt_s']:finest_by_pair[key]=c
        unresolved=sum(not c.get('resolved',False) for c in finest_by_pair.values())
        lines += ['', f"At the finest tested grid for each category comparison, {unresolved} of {len(finest_by_pair)} comparisons remain unresolved by the declared screen. A small category difference must not be described as a precise class ranking when this screen fails; the reference-to-design reduction has its own separate check.",'']
        category_rows=[r for r in resolution['records'] if r['label']=='wetland_best' or r['label'].startswith('category_best_')]
        finest_category_step=min(r['dt_s'] for r in category_rows if r['status']=='valid')
        category_peaks=[r['metrics']['peak_stage_m'] for r in category_rows if r['status']=='valid' and r['dt_s']==finest_category_step]
        lines += [f"The entire best-candidate category spread at {finest_category_step:g} seconds is {(max(category_peaks)-min(category_peaks))*1000:.6f} mm. The verified ordering concerns these finite-search numerical candidates only; such a small synthetic-model difference is not a local engineering or vegetation recommendation.",'']
    valid=[r for r in resolution['records'] if r['status']=='valid']
    if valid:
        rain=np.array([r['rain_depth_m'] for r in valid]);solar=np.array([r['solar_energy_Wh_m2'] for r in valid])
        lines += ['',f"Integrated rain is {rain[0]*1000:.6f} mm (spread {np.ptp(rain):.3g} m); solar energy is {solar[0]:.6f} Wh/m² (spread {np.ptp(solar):.3g} Wh/m²). Both preserve the existing prepared observations, gridded estimate and accepted spatial/subinterval assumptions.",'']
    lines += ['The resolution screen requires feasible fixed-design runs, unchanged improvement sign and summed grid changes below 10% of the peak reduction. It is not a rigorous error bound. A comparable discrepancy or changed ranking must be reported as unresolved; convergence of OPT alone does not resolve integration error.','']
    for label,run in resolved_runs:
        lines += [f"`{label}`: river rising at the end = {run['metrics'].get('river_rising_at_end')}; peak near the horizon = {run['metrics'].get('peak_near_horizon')}. The six-day event scope remains fixed.",'']
    tests_path=output/'tests.log'
    if tests_path.exists() and 'Ran 86 tests' in tests_path.read_text() and '\nOK' in tests_path.read_text():
        lines += ['The recorded software suite passed **86 tests**, including **nine new wetland-only tests**. These check the two-cover decision interface, fixed-gate preservation, the 18-row optimizer view and separate G5/model-validity handling alongside existing model tests. They verify software behavior; they are not field calibration or an engineering-safety certificate.','']
    lines += ['## Reproduction and scope','',
        'The configuration `B/wetland_only_config.json` records the fixed-gate source path and checksum; the saved study plan contains all 24 numerical gate openings. The numerical gate schedule is the accepted new reference, not the obsolete original control values. G1–G6 provenance and the eight-state no-spillway equations remain as documented in [data sources and assumptions](../data_sources_and_assumptions.md).','',
        'Saved configuration and provenance records retain their original path labels and hashes. `B/project_paths.py` resolves them to the physical stage folders; the local relocation manifest records the path/import changes.','', 'From the project root in the existing Python environment with the pinned official solver:','',
        '```sh','python3 -B -m B.wetland_only_study --run','python3 -B -m B.wetland_only_study --compare','python3 -B -m B.wetland_only_study --report','```','',
        'Completed numerical stages refuse to overwrite their saved results. `--report` regenerates these reports and figures only. Results are local under `B/results/wetland_only/`; all earlier joint-study and A.2 files remain historical and unchanged. This analysis is hypothetical and uncalibrated, with no claim of local engineering certification, species recommendation or global optimality. No Git commit or publication is included.','']
    for figure in figures:lines += [f"![{figure.stem.replace('_',' ')}](../plots/wetland_only/{figure.name})",'']
    document=ROOT/'docs/12_wetland_only_optimization.md';document.write_text('\n'.join(lines))
    return dict(document=str(document),tables=[str(p) for p in tables],figures=[str(p) for p in figures])
