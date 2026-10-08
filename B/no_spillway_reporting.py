"""Reports for the separate B.2 revision without pool spillways.

Reads completed numerical artifacts. No optimizer, solver or scientific parameter
is changed or executed by this module's reporting functions.
"""
from __future__ import annotations

from collections import Counter
import csv
import json
from pathlib import Path
from B.project_paths import resolve_path
import numpy as np
from A import wetland_flood_model as a2

ROOT=Path(__file__).resolve().parent
OUTPUT=ROOT/'results/no_spillway_revision'


def mechanism_metrics(run):
    """Physical routing diagnostics for one complete trajectory; no simulation.

    Integral-based diagnostics use reporting-grid trapezoidal quadrature. Flow
    centroid is an event-time first moment, not a parcel travel/residence time.
    Internal gate, drainage, overflow and groundwater transfers are excluded
    from the whole-system external-water budget.
    """
    t=np.asarray(run['time'],dtype=float);x=np.asarray(run['x']);y=np.asarray(run['y']);u=np.asarray(run['u'])
    if t.ndim!=1 or len(t)<2 or np.any(np.diff(t)<=0):
        raise ValueError('Mechanism diagnostics require a complete increasing time grid.')
    if any(a.shape[-1]!=len(t) or not np.all(np.isfinite(a)) for a in (x,y,u)):
        raise ValueError('Mechanism diagnostics cannot complete an invalid trajectory.')
    integrate=lambda q: float(np.trapezoid(q,t))
    metric={}
    def flow(prefix,q):
        q=np.asarray(q);volume=integrate(q);k=int(np.argmax(q))
        metric[prefix+'_volume_m3']=volume
        metric[prefix+'_peak_m3_s']=float(q[k])
        metric[prefix+'_peak_time_h']=float(t[k]/3600) if q[k]>0 else None
        metric[prefix+'_centroid_time_h']=integrate(t*q)/volume/3600 if volume>0 else None
    for i in (1,2):
        for prefix,name in [('gate','Q_G'),('normal_drainage','Q_wr'),('wetland_overflow','Q_wsp'),
                            ('wetland_seepage','Q_wg'),('wetland_et','E_w'),('dam_spillway','Q_dsp')]:
            flow(f'{prefix}{i}',y[a2.Y[f'{name}{i}']])
        sw=x[1+i];k=int(np.argmax(sw));stage=y[a2.Y[f'H_w{i}']]
        metric[f'wetland{i}_peak_storage_m3']=float(sw[k]);metric[f'wetland{i}_peak_storage_time_h']=float(t[k]/3600)
        metric[f'wetland{i}_peak_level_m']=float(np.max(stage))
        metric[f'wetland{i}_initial_storage_m3']=float(sw[0]);metric[f'wetland{i}_final_storage_m3']=float(sw[-1])
        metric[f'wetland{i}_storage_change_m3']=float(sw[-1]-sw[0])
        metric[f'wetland_overflow{i}_positive_duration_h']=integrate((y[a2.Y[f'Q_wsp{i}']]>0).astype(float))/3600
        metric[f'pool{i}_peak_level_m']=float(np.max(y[a2.Y[f'H_d{i}']]))
        metric[f'pool{i}_storage_change_m3']=float(x[3+i,-1]-x[3+i,0])
    for prefix,names in [('gate_total',['Q_G1','Q_G2']),('normal_drainage_total',['Q_wr1','Q_wr2']),
                         ('wetland_overflow_total',['Q_wsp1','Q_wsp2']),
                         ('wetland_to_river_total',['Q_wr1','Q_wr2','Q_wsp1','Q_wsp2']),
                         ('river_outlet',['Q_out']),('upland_bypass',['Q_ur'])]:
        flow(prefix,sum(y[a2.Y[name]] for name in names))
    hr=y[a2.Y['H_r']];k=int(np.argmax(hr))
    metric['river_peak_m']=float(hr[k]);metric['river_peak_time_h']=float(t[k]/3600)
    metric['river_peak_near_horizon']=bool(t[k]>=.95*t[-1])
    metric['river_rising_at_end']=bool(hr[-1]>hr[-2])
    metric['rain_volume_m3']=integrate(y[a2.Y['P_total']]);metric['external_surface_inflow_volume_m3']=integrate(u[a2.U['Q_ext']])
    metric['external_groundwater_net_volume_m3']=integrate(y[a2.Y['Q_gx']])
    metric['external_groundwater_in_volume_m3']=integrate(y[a2.Y['Q_gx_in']])
    metric['external_groundwater_out_volume_m3']=integrate(y[a2.Y['Q_gx_out']])
    metric['evapotranspiration_volume_m3']=integrate(y[a2.Y['E_tot']])
    metric['total_storage_change_m3']=float(x[:,-1].sum()-x[:,0].sum())
    net=y[a2.Y['P_total']]+u[a2.U['Q_ext']]+y[a2.Y['Q_gx']]-y[a2.Y['Q_out']]-y[a2.Y['E_tot']]
    residual=x.sum(axis=0)-x[:,0].sum()-a2.cumulative_trapezoid(net,t)
    metric['net_external_input_volume_m3']=integrate(net)
    metric['external_balance_end_residual_m3']=float(residual[-1])
    metric['external_balance_max_abs_residual_m3']=float(np.max(np.abs(residual)))
    metric['external_net_output_identity_max_abs_m3_s']=float(np.max(np.abs(net-y[a2.Y['external_net']])))
    metric['dam_spillway_max_abs_m3_s']=float(np.max(np.abs(y[[a2.Y['Q_dsp1'],a2.Y['Q_dsp2']]])))
    metric['rain_depth_m']=integrate(u[a2.U['P']])
    metric['solar_energy_Wh_m2']=integrate(u[a2.U['R_s']])/3600
    return metric


def _load_run(path):
    with np.load(path) as data:
        run={name:np.array(data[name]) for name in ['time','x','u','y','theta','residuals','balance'] if name in data}
        run.update(json.loads(str(data['metadata_json'])))
    return run


def _csv(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    fields=list(dict.fromkeys(k for row in rows for k in row))
    if not fields:fields=['status']
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    return path


def _pair(pair):return '/'.join(p.removeprefix('type_').upper() for p in pair)

def _json(value):return json.dumps(value,default=lambda v:v.tolist() if isinstance(v,np.ndarray) else v.item())

def _fmt(value,n=6):return '—' if value is None else f'{value:.{n}f}'


def export_study(config,output_dir=None):
    """Read completed revision artifacts and write only this revision's reports."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    output=OUTPUT if output_dir is None else Path(output_dir)
    tabledir=ROOT/'tables/no_spillway_revision';plotdir=ROOT/'plots/no_spillway_revision'
    tabledir.mkdir(parents=True,exist_ok=True);plotdir.mkdir(parents=True,exist_ok=True)
    def read(name,default):
        path=output/name
        return json.loads(path.read_text()) if path.exists() else default
    search=read('reference_search.json',{'records':[],'selected_label':None,'selection':'Not available'})
    screen=read('start_screening.json',{'records':[],'accepted':[]})
    batch=read('batch_status.json',{'runs':[],'status':'No completed optimization batch'})
    plan=read('batch_plan.json',{})
    resolution=read('fixed_design_timesteps.json',{'records':[],'comparisons':[],'thirty_second_required':False})
    controlled=read('controlled_comparisons.json',{'records':[],'blocks':[]})
    reference=_load_run(output/'reference.npz') if (output/'reference.npz').exists() else None
    best=_load_run(output/'batch_best.npz') if (output/'batch_best.npz').exists() else None
    reports=[]
    for path in sorted(output.glob('*_run.json')):
        record=json.loads(path.read_text())
        if 'physical_callback_log' in record:reports.append(record)
    run_rows=[];constraint_rows=[];mechanism_rows=[];controlled_rows=[];controlled_loaded=[]
    def run_is_feasible(run):return bool(run and all(r['pass_at_tolerance'] for r in run['audit']))
    refpeak=None if reference is None else reference['metrics']['peak_stage_m']
    identical_joint_reference=bool(reference is not None and best is not None and reference['pair']==best['pair'] and
                                   np.allclose(reference['theta'],best['theta'],rtol=0,atol=3e-15))
    named=[]
    if reference is not None:named.append(('new_reference',reference))
    if best is not None:named.append(('joint_best_found',best))
    def register(label,run,block='joint_comparison'):
        metric=mechanism_metrics(run)
        if metric['dam_spillway_max_abs_m3_s']!=0:
            raise ValueError('Nonzero dam spillway in purported no-spillway trajectory: '+label)
        mechanism_rows.append(dict(block=block,label=label,pair=_pair(run['pair']),
            feasible=run_is_feasible(run),dt_s=float(np.diff(run['time'])[0]),
            cover1=float(run['theta'][-2]),cover2=float(run['theta'][-1]),**metric))
        for row in run['audit']:
            constraint_rows.append(dict(block=block,run=label,pair=_pair(run['pair']),
                physical_residual=-row['margin'],**row))
        return metric
    named_mechanism={label:register(label,run) for label,run in named}
    for r in reports:
        metrics=r.get('best_metrics');audit={a['label']:a for a in (r.get('best_audit') or [])}
        initial=r['requested_initial_design'];first=r.get('first_optimizer_design',[]);final=r.get('best_design')
        callbacks=r['physical_callback_log']
        origin=next((item['evaluation'] for item in callbacks if final is not None and item['feasible'] and
                     np.array_equal(np.asarray(json.loads(item['decision_json'])),np.asarray(final))),None)
        row=dict(pair=_pair(r['category_pair']),start_id=r['start_id'],
            initial_design_json=_json(initial),first_design_json=_json(first),final_design_json=_json(final),
            initial_cover1=initial[-2],initial_cover2=initial[-1],
            final_cover1=None if final is None else final[-2],final_cover2=None if final is None else final[-1],
            first_design_feasible=r.get('first_design_feasible'),
            max_first_design_shift=float(np.max(np.abs(np.asarray(initial)-first))) if first else None,
            physical_peak_m=None if metrics is None else metrics['peak_stage_m'],fixed_new_reference_peak_m=refpeak,
            peak_reduction_m=None if metrics is None or refpeak is None else refpeak-metrics['peak_stage_m'],
            max_scaled_violation=None if not audit else max(0.,max(a['scaled_residual'] for a in audit.values())),
            callbacks=r['callbacks'],cached_callbacks=r['cached_callbacks'],actual_ode_attempts=r['actual_ode_attempts'],
            valid_ode_completions=sum(not x['cached'] and x['raw_peak_stage_m'] is not None for x in callbacks),
            hydraulically_invalid_callbacks=sum(x.get('hydraulically_invalid',False) for x in callbacks),
            valid_but_infeasible_callbacks=sum(not x.get('failure') and not x['feasible'] for x in callbacks),
            best_origin_callback=origin,best_is_initial_optimizer_design=(origin==1),
            rejection_counts_json=_json(r['rejection_counts']),termination=r['termination'],elapsed_seconds=r['seconds'])
        for label in ['H_d1_safe','H_d2_safe','H_w1_safe','H_w2_safe','gate1_head_opening','gate2_head_opening','soil_capacity']:
            row[label+'_margin']=audit.get(label,{}).get('margin')
        run_rows.append(row)
    # A controlled category block must preserve every gate knot AND both covers.
    block_designs={};block_runs={}
    for r in controlled['records']:
        block=r['block'];theta=np.asarray(r['theta'],dtype=float)
        if block in block_designs and not np.allclose(theta,block_designs[block],rtol=0,atol=1e-14):
            raise ValueError('Controlled comparison changes gates/covers within block '+block)
        block_designs.setdefault(block,theta.copy())
        row=dict(block=block,label=r['label'],pair=_pair(r['pair']),cover1=r['covers'][0],cover2=r['covers'][1],
                 status=r['status'],feasible=r['feasible'],theta_json=_json(r['theta']),
                 failure=r.get('failure',''),diagnostic_json=_json(r.get('diagnostic')))
        run=None
        trajectory=r.get('trajectory_file')
        if trajectory:
            candidates=[resolve_path(trajectory),Path(trajectory),output/trajectory,ROOT/trajectory,output/Path(trajectory).name]
            path=next((p for p in candidates if p.is_file()),None)
            if path is None:raise FileNotFoundError('Reported controlled trajectory is missing: '+str(trajectory))
            run=_load_run(path)
            if not np.allclose(run['theta'],theta,rtol=0,atol=1e-14):raise ValueError('Controlled metadata/design mismatch.')
            metric=register(r['label'],run,block);row.update(metric)
            block_runs.setdefault(block,[]).append((r,run,metric))
        elif r['status']=='valid':
            raise ValueError('Valid controlled comparison must retain its actual trajectory.')
        controlled_loaded.append((r,run));controlled_rows.append(row)
    for row in controlled_rows:
        anchor=next((x for x in controlled_rows if x['block']==row['block'] and x['pair']=='A/A' and 'river_peak_m' in x),None)
        row['peak_change_from_same_block_AA_m']=None if anchor is None or 'river_peak_m' not in row else row['river_peak_m']-anchor['river_peak_m']
    search_rows=[]
    for r in search['records']:
        row={k:v for k,v in r.items() if k not in ['theta','pair','metrics','audit','diagnostic']}
        row.update(theta_json=_json(r.get('theta')),pair=_pair(r.get('pair',r.get('category_pair',[]))),
            physical_peak_m=r.get('metrics',{}).get('peak_stage_m'),diagnostic_json=_json(r.get('diagnostic')))
        for a in r.get('audit',[]):
            if a['label'].endswith('_safe') or 'head_opening' in a['label']:row[a['label']+'_margin']=a['margin']
        search_rows.append(row)
    timestep_rows=[]
    for r in resolution['records']:
        row={k:v for k,v in r.items() if k not in ['metrics','audit','design','pair']}
        row.update(pair=_pair(r.get('pair',[])),design_json=_json(r.get('design')),**r.get('metrics',{}))
        row['max_scaled_violation']=max(0.,max(a['scaled_residual'] for a in r['audit'])) if r.get('audit') else None
        timestep_rows.append(row)
    tables=[_csv(tabledir/'reference_search.csv',search_rows),_csv(tabledir/'optimization_summary.csv',run_rows),
            _csv(tabledir/'mechanism_metrics.csv',mechanism_rows),_csv(tabledir/'constraint_residuals.csv',constraint_rows),
            _csv(tabledir/'controlled_comparisons.csv',controlled_rows),_csv(tabledir/'fixed_design_timesteps.csv',timestep_rows)]

    figures=[]
    def finish(fig,name):
        narrow=fig.get_size_inches()[0]<8
        note=('No dam spillways; wetland overflow retained.\nHypothetical installation and synthetic response classes.' if narrow else
              'No dam spillways; wetland overflow retained. Hypothetical installation and synthetic response classes.')
        fig.text(.5,.006,note,ha='center',fontsize=8)
        fig.tight_layout(rect=(0,.04 if narrow else .022,1,.945 if narrow else .97));path=plotdir/name;fig.savefig(path,dpi=170);plt.close(fig);figures.append(path)
    if named:
        fig,axes=plt.subplots(3,2,figsize=(12,10),sharex=True)
        colors=['#2166ac','#b35806'];styles=['-','--']
        for ri,(label,run) in enumerate(named):
            t=run['time']/3600;y=run['y'];x=run['x'];color=colors[ri];title=label.replace('_',' ')
            axes[0,0].plot(t,y[a2.Y['H_r']],color=color,label=title)
            for i in range(2):
                suffix=f' / branch {i+1}';style=styles[i]
                axes[0,1].plot(t,run['u'][i],color=color,linestyle=style,label=title+suffix)
                axes[1,0].plot(t,y[a2.Y[f'H_d{i+1}']],color=color,linestyle=style,label=title+suffix)
                axes[1,1].plot(t,x[2+i],color=color,linestyle=style,label=title+suffix)
            axes[2,0].plot(t,y[a2.Y['Q_G1']]+y[a2.Y['Q_G2']],color=color,label=title)
            axes[2,1].plot(t,y[a2.Y['Q_wr1']]+y[a2.Y['Q_wr2']],color=color,label=title+' / normal drainage')
            axes[2,1].plot(t,y[a2.Y['Q_wsp1']]+y[a2.Y['Q_wsp2']],color=color,linestyle='--',label=title+' / wetland overflow')
        for i,limit in enumerate(config['safety']['dam_m']):axes[1,0].axhline(limit,color='#777777',linestyle=styles[i],linewidth=.8)
        labels=['River level (m)','Gate opening fraction','Pool levels (m; gray lines = limits)',
                'Wetland storage (m³)','Combined gate release (m³/s)','Wetland-to-river flows (m³/s)']
        for ax,label in zip(axes.flat,labels):
            ax.set_ylabel(label);ax.grid(alpha=.2);ax.legend(fontsize=7);ax.ticklabel_format(axis='y',useOffset=False,style='plain')
        for ax in axes[-1]:ax.set_xlabel('Elapsed time (h)')
        fig.suptitle('Reference and retained best coincide' if identical_joint_reference else 'New reference and joint best: actual trajectories and routing mechanisms')
        finish(fig,'01_joint_trajectories_and_mechanisms.png')
    if block_runs:
        block_ids=list(block_runs);fig,axes=plt.subplots(3,len(block_ids),figsize=(6.2*len(block_ids),9),sharex=True,squeeze=False)
        colors={'A/A':'#2166ac','A/B':'#1b9e77','B/A':'#e6ab02','B/B':'#b35806'}
        for col,block in enumerate(block_ids):
            for record,run,_ in block_runs[block]:
                pair=_pair(run['pair']);label=pair+(' (infeasible)' if not record['feasible'] else '')
                for row,name in enumerate(['H_r','H_w1','H_w2']):
                    axes[row,col].plot(run['time']/3600,run['y'][a2.Y[name]],color=colors.get(pair),
                        linestyle='-' if record['feasible'] else '--',label=label)
            theta=block_designs[block];axes[0,col].set_title(f'{block}\nFixed covers {theta[-2]:.4f}, {theta[-1]:.4f}')
            for row,name in enumerate(['River level (m)','Wetland 1 level (m)','Wetland 2 level (m)']):
                ax=axes[row,col];ax.set_ylabel(name);ax.grid(alpha=.2);ax.legend(fontsize=8)
                ax.ticklabel_format(axis='y',useOffset=False,style='plain')
                if row:ax.axhline(config['safety']['wetland_m'][row-1],color='gray',linestyle=':',linewidth=.8)
            axes[-1,col].set_xlabel('Elapsed time (h)')
        fig.suptitle('Controlled wetland-class comparison\nIdentical gate schedules and covers within each block',fontsize=13)
        finish(fig,'02_controlled_wetland_and_river_levels.png')
    if reports:
        fig,axes=plt.subplots(1,2,figsize=(13,4.8));colors={'A/A':'#2166ac','A/B':'#1b9e77','B/A':'#e6ab02','B/B':'#b35806'}
        paircounts=Counter();ticks=[];labels=[]
        for i,r in enumerate(reports):
            pair=_pair(r['category_pair']);paircounts[pair]+=1;label=f'{pair} S{paircounts[pair]}'
            ticks.append(i);labels.append(label)
            if r.get('best_metrics'):axes[0].scatter(i,r['best_metrics']['peak_stage_m'],color=colors.get(pair),s=40)
            count=0;xs=[];ys=[]
            for record in r['physical_callback_log']:
                if record['cached']:continue
                count+=1
                if record['best_feasible_peak_m'] is not None:xs.append(count);ys.append(record['best_feasible_peak_m'])
            axes[1].step(xs,ys,where='post',color=colors.get(pair),linestyle=['-','--',':'][(paircounts[pair]-1)%3],label=label)
        if refpeak is not None:
            for ax in axes:ax.axhline(refpeak,color='black',linestyle='--',linewidth=.8)
        axes[0].set_xticks(ticks);axes[0].set_xticklabels(labels,rotation=45,ha='right',fontsize=8)
        axes[0].set_title('Best feasible physical objective by start');axes[0].set_xlabel('Actual starts in each synthetic category')
        axes[1].set_title('Feasible-best progress');axes[1].set_xlabel('Actual ODE attempts; cache hits excluded');axes[1].legend(fontsize=7,ncol=2)
        for ax in axes:
            ax.set_ylabel('Physical river peak (m)');ax.ticklabel_format(axis='y',useOffset=False,style='plain');ax.grid(alpha=.2)
        fig.suptitle('Bounded optimization robustness under the revised topology')
        finish(fig,'03_search_robustness.png')

    # Narrative records physical findings separately from algorithm penalties.
    lines=['# B.2 revision without dam spillways','',
        'This new study removes the two pool-to-wetland dam spillways from B.2. Each pool releases to its wetland through its gate only. Wetland overflow to the river, normal wetland drainage, the direct upland-to-river bypass, soil and groundwater pathways remain. There is no replacement emergency spill, automatic gate control or water-level clipping. The unchanged safety limits and hydraulic-applicability constraints determine feasibility. A.2 and all historical B.2 results are separate and preserved.','',
        'The installation remains hypothetical and uncalibrated, with regional event weather and synthetic A/B vegetation-response classes. Read [data sources and assumptions](../data_sources_and_assumptions.md) for the active topology, unchanged provenance, design assumptions and interpretation limits. Older studies with dam spillways cannot supply numerical results for this topology.','',
        '## Reference selection and joint optimization','',
        f"Reference-selection method: {search.get('selection','Not recorded')}. Selected label: `{search.get('selected_label')}`. Tested candidates: {len(search['records'])}. The complete search outcomes and first failure diagnostics are in [reference_search.csv](../tables/no_spillway_revision/reference_search.csv).",'']
    rejected_path=output/'rejected_reference.npz'
    if rejected_path.exists():
        old=_load_run(rejected_path)
        lines += ['The unchanged previous small-opening reference was first rerun under the **new** topology. It completed the hydraulic integration but exceeded pool safety limits; those values were not clipped. First reporting-grid exceedances are shown below, rather than interpreted as exact continuous-time crossing times.','',
                  '| Pool | First violating time (h) | Level (m) | Limit (m) | Residual (m) | Maximum exceedance (m) |',
                  '|---|---:|---:|---:|---:|---:|']
        for i,limit in enumerate(config['safety']['dam_m']):
            h=old['y'][a2.Y[f'H_d{i+1}']];indices=np.flatnonzero(h>limit)
            if len(indices):
                k=indices[0]
                lines.append(f'| {i+1} | {old["time"][k]/3600:.6f} | {h[k]:.9f} | {limit:.2f} | {h[k]-limit:.9f} | {np.max(h)-limit:.6f} |')
        lines += ['']
    first_hydraulic=next((r for r in search['records'] if r.get('diagnostic') and r['diagnostic'].get('branch') is not None),None)
    if first_hydraulic:
        d=first_hydraulic['diagnostic']
        lines += [f"A distinct failed candidate, `{first_hydraulic['label']}`, stopped at its first invalid gate evaluation: branch {d['branch']}, {d['time_s']/3600:.6f} h, depth {d['water_depth_m']:.9f} m, opening height {d['opening_height_m']:.9f} m and G5 residual +{d['residual_m']:.9f} m. It is hydraulically invalid, unlike the complete but unsafe old reference.",'']
    if reference is None:
        lines += ['**No complete selected reference trajectory is available in this revision’s results.** This report does not invent a baseline, optimized improvement or proof of infeasibility. Candidate failures in the recorded bounded search remain the outcome.','']
    else:
        dt=float(np.diff(reference['time'])[0])
        lines += [f"The newly evaluated reference ({_pair(reference['pair'])}, covers {reference['theta'][-2]:.4f}/{reference['theta'][-1]:.4f}) has physical river peak **{refpeak:.9f} m** at {reference['metrics']['peak_time_s']/3600:.4f} h on its {dt:g}-second grid; feasibility: **{run_is_feasible(reference)}**. It is identified by its selection method above, not relabeled as the original baseline.",'']
    if best is None:
        lines += ['No best-feasible optimized trajectory is present. This means none was retained by the recorded search, not that the design problem has been proven infeasible. Safety limits were not raised to obtain a result.','']
    else:
        peak=best['metrics']['peak_stage_m'];delta=None if refpeak is None else refpeak-peak
        lines += [f"The best feasible joint design found is {_pair(best['pair'])}, covers {best['theta'][-2]:.6f}/{best['theta'][-1]:.6f}, with physical river peak **{peak:.9f} m**. Its same-study reference reduction is **{_fmt(None if delta is None else delta*1000,4)} mm** on the saved reporting grid. This joint change can include gate operations, covers and category choices; it is not a vegetation-only effect.",'']
    if identical_joint_reference:
        lines += ['**No improvement over the new reference was found.** The best feasible design retained by the bounded optimizer is the same reference design (up to coordinate round-trip precision), so their trajectories coincide. This is a completed search attempt with a zero-improvement outcome, not a newly improved operating policy.','']
    if best is not None:
        owners=[r for r in reports if r.get('best_design') is not None and r['category_pair']==best['pair']
                and np.array_equal(np.asarray(r['best_design']),np.asarray(best['theta']))]
        if owners:
            owner=owners[0];entry=next(row for row in run_rows if row['pair']==_pair(owner['category_pair']) and row['start_id']==owner['start_id'])
            origin=entry['best_origin_callback']
            if origin==1:
                lines += [f"**The selected joint candidate originated at OPT evaluation 1**, the already screened `{owner['start_id']}` initial design. This run retained that point rather than improving it. Its difference from the new reference is evidence about that candidate design; it must not be reported as an objective reduction produced by successful OPT steps.",'']
            else:
                lines += [f"The selected joint candidate first appeared at callback {origin} of `{owner['start_id']}`. The CSV distinguishes this origin from the requested and first-evaluated design.",'']
    lines += ['| Classes | Start | Final covers | Physical peak (m) | Reduction from new reference (mm) | ODE attempts | Termination |',
              '|---|---|---|---:|---:|---:|---|']
    for row in run_rows:
        lines.append(f"| {row['pair']} | {row['start_id']} | {_fmt(row['final_cover1'],4)}/{_fmt(row['final_cover2'],4)} | {_fmt(row['physical_peak_m'],9)} | {_fmt(None if row['peak_reduction_m'] is None else row['peak_reduction_m']*1000,4)} | {row['actual_ode_attempts']} | {row['termination']} |")
    callbacks=sum(r['callbacks'] for r in reports);cache=sum(r['cached_callbacks'] for r in reports);odes=sum(r['actual_ode_attempts'] for r in reports)
    failures=Counter();termination=Counter(r['termination'] for r in reports)
    for r in reports:failures.update(r['rejection_counts'])
    retained_initial=sum(row['best_is_initial_optimizer_design'] for row in run_rows)
    only_first_feasible=sum(bool(r['physical_callback_log']) and r['physical_callback_log'][0]['feasible'] and
                            not any(item['feasible'] for item in r['physical_callback_log'][1:]) for r in reports)
    lines += [f"{retained_initial} of {len(reports)} runs retained their initial physical design as their best feasible point; {only_first_feasible} had no feasible callback after evaluation 1. This distinguishes a search that retained screened starts from an optimizer that found additional feasible improvements.",'']
    hydraulic_diagnostics=[item['diagnostic'] for r in reports for item in r['physical_callback_log']
                           if item.get('hydraulically_invalid') and item.get('diagnostic')]
    if hydraulic_diagnostics:
        times=[d['time_s']/3600 for d in hydraulic_diagnostics]
        branches=Counter(d.get('branch') for d in hydraulic_diagnostics)
        lines += [f"The recorded hydraulic rejections occurred at {min(times):.4f}–{max(times):.4f} h, with branch counts {dict(branches)}. These first-invalid times identify where the trial trajectories leave the gate-applicability domain; they are distinct from a complete-event safety failure or a demonstrated global feasibility boundary.",'']
    reference_odes=search.get('actual_ode_attempts',len(search['records']))
    screening_odes=screen.get('actual_ode_attempts',sum(r.get('ode_attempted',True) for r in screen['records']))
    fixed_odes=sum(r.get('ode_attempted',not bool(r.get('simulation_reused_from'))) for r in resolution['records'])
    controlled_odes=sum(r.get('ode_attempted',True) for r in controlled['records'])
    total_study_odes=reference_odes+screening_odes+odes+fixed_odes+controlled_odes
    lines += [f"The study used **{total_study_odes} actual ODE attempts**: {reference_odes} reference-search, {screening_odes} start-screening, {odes} optimizer, {fixed_odes} distinct fixed-design timestep and {controlled_odes} controlled-class runs. Reused identical-design rows are explicitly labeled and are not additional integrations. Unit-test calls are separate from these study simulations.",'']
    lines += ['',f"Optimizer counts: {callbacks} callbacks, {cache} cache hits, {odes} ODE attempts. Terminations: {dict(termination)}. Rejected callbacks by diagnostic: {dict(failures)}. Full designs, requested/first-evaluated starts, physical margins and counts are in [optimization_summary.csv](../tables/no_spillway_revision/optimization_summary.csv). An invalid hydraulic evaluation stops that candidate; it has no complete physical peak. A valid trajectory exceeding a safety limit is separately infeasible.",'',
        f"Recorded batch status: `{batch.get('status')}`. Declared optimizer settings: `{plan.get('settings','See individual run JSON records')}`. Screening records: {len(screen['records'])}, accepted starts: {len(screen['accepted'])}. Budget/stagnation termination is not convergence or global optimality; the official 26-dimensional NMS convergence test additionally requires more than 676 iterations. Multiple starts provide evidence about this bounded search only.",'',
        'The decision vector remains 24 independently variable gate knots and two covers. The official native bounds are [−1,1], mapped onto the full physical [0,1] box. The official initial-point clip is addressed with one fixed map for this entire new batch: `g = q³ / [q³ + (1−q)³]`, while cover equals q. This bijection accommodates both small dry-period openings and larger storm releases without narrowing physical bounds. Requested and first-evaluated physical designs are recorded separately; the transformation changes search conditioning and is not a physical model parameter. It differs from the historical low-opening study’s rational map, so cross-topology optimizer outcomes are not a controlled coordinate-map comparison.','',
        '## Routing mechanisms and conservation','',
        'The following diagnostics use complete simulated trajectories. Centroid time is the first temporal moment of release, not a water-parcel residence time. Volumes and moments use reporting-grid trapezoidal quadrature.','',
        '| Design | Gate volume (m³) | Gate centroid (h) | Wetland normal-drainage volume (m³) | Wetland overflow volume (m³) | River peak time (h) | Δ total storage (m³) | Max external-balance residual (m³) |',
        '|---|---:|---:|---:|---:|---:|---:|---:|']
    for label,metric in named_mechanism.items():
        lines.append(f"| {label} | {metric['gate_total_volume_m3']:.3f} | {_fmt(metric['gate_total_centroid_time_h'],4)} | {metric['normal_drainage_total_volume_m3']:.3f} | {metric['wetland_overflow_total_volume_m3']:.3f} | {metric['river_peak_time_h']:.4f} | {metric['total_storage_change_m3']:.3f} | {metric['external_balance_max_abs_residual_m3']:.6f} |")
    for label,run in named:
        lines += [f"For `{label}`, river rising at the final time: **{run['metrics'].get('river_rising_at_end')}**; event peak near the horizon: **{run['metrics'].get('peak_near_horizon')}**. The six-day horizon remains fixed; a rising endpoint warrants that scope limitation even when the recorded maximum occurred earlier."]
    lines += ['', 'Branch-specific gate volumes/centroids/peak times, wetland peak storage and timing, normal drainage, overflow, seepage, ET and final-storage changes are in [mechanism_metrics.csv](../tables/no_spillway_revision/mechanism_metrics.csv). Dam-spillway outputs are verified identically zero for every reported complete trajectory. Wetland overflow remains in the model even when its simulated event volume happens to be zero.','',
        'Whole-system conservation uses only `P_total + Q_ext + Q_gx − Q_out − E_tot`. Gate releases, wetland drainage/overflow, upland bypass and internal soil/groundwater transfers cancel and are not counted again as external water. A change in final storage indicates redistribution within the finite event; it does not establish reduced discharge over an indefinitely extended horizon.','',
        '## Controlled wetland-class comparisons','',
        'Within each block, every gate-knot value, both vegetation covers, initial conditions, forcing and fixed parameters are held identical. Only A/A, A/B, B/A or B/B changes. This separates the synthetic class package effect from joint gate/cover optimization. A second block, if present, uses the joint-best gates/covers and remains distinct from the reference block. All four class combinations are checked independently; an invalid or unsafe comparison is not silently replaced.','',
        '| Fixed-control block | Classes | Covers | Status / feasible | Physical peak (m) | Peak change from block A/A (mm) | Normal drainage (m³) | Wetland overflow (m³) |',
        '|---|---|---|---|---:|---:|---:|---:|']
    for row in controlled_rows:
        lines.append(f"| {row['block']} | {row['pair']} | {row['cover1']:.4f}/{row['cover2']:.4f} | {row['status']} / {row['feasible']} | {_fmt(row.get('river_peak_m'),9)} | {_fmt(None if row['peak_change_from_same_block_AA_m'] is None else 1000*row['peak_change_from_same_block_AA_m'],4)} | {_fmt(row.get('normal_drainage_total_volume_m3'),3)} | {_fmt(row.get('wetland_overflow_total_volume_m3'),3)} |")
    lines += ['']
    for block,group in block_runs.items():
        feasible=[(r,run,m) for r,run,m in group if r['feasible']]
        if len(feasible)>1:
            levels=[m['river_peak_m'] for _,_,m in feasible];volumes=[m['gate_total_volume_m3'] for _,_,m in feasible]
            peakrange=max(levels)-min(levels);release_range=max(volumes)-min(volumes)
            low=min(feasible,key=lambda x:x[2]['river_peak_m']);high=max(feasible,key=lambda x:x[2]['river_peak_m'])
            lines += [f"For block `{block}`, the feasible class-substitution peak range is {1000*peakrange:.4f} mm ({_pair(low[1]['pair'])} lowest; {_pair(high[1]['pair'])} highest). The range in integrated gate release is {release_range:.6g} m³ despite fixed openings; the table records the resulting drainage/overflow differences rather than assigning the entire joint benefit to vegetation.",'']
    for block,group in block_runs.items():
        aa=next(((r,run,m) for r,run,m in group if _pair(run['pair'])=='A/A' and r['feasible']),None)
        bb=next(((r,run,m) for r,run,m in group if _pair(run['pair'])=='B/B' and r['feasible']),None)
        if aa is not None and bb is not None:
            _,ar,am=aa;_,br,bm=bb
            gate_diff=float(np.max(np.abs(ar['y'][[a2.Y['Q_G1'],a2.Y['Q_G2']]]-br['y'][[a2.Y['Q_G1'],a2.Y['Q_G2']]])))
            lines += [f"**Mechanism in `{block}` (A/A versus B/B, fixed gates and covers):** the maximum pointwise difference in either gate discharge is {gate_diff:.3g} m³/s, and the common gate-release volume is {am['gate_total_volume_m3']:.3f} m³ with centroid {am['gate_total_centroid_time_h']:.4f} h. Thus this comparison changes the downstream wetland response while holding the modeled release histories fixed.",'',
                f"Changing A/A to B/B raises the branch-1 peak wetland storage from {am['wetland1_peak_storage_m3']:.2f} to {bm['wetland1_peak_storage_m3']:.2f} m³ and branch 2 from {am['wetland2_peak_storage_m3']:.2f} to {bm['wetland2_peak_storage_m3']:.2f} m³. Normal-drainage volume changes from {am['normal_drainage_total_volume_m3']:.2f} to {bm['normal_drainage_total_volume_m3']:.2f} m³, while retained wetland overflow changes from {am['wetland_overflow_total_volume_m3']:.2f} to {bm['wetland_overflow_total_volume_m3']:.2f} m³.",'',
                f"Normal-drainage peak decreases from {am['normal_drainage_total_peak_m3_s']:.6f} to {bm['normal_drainage_total_peak_m3_s']:.6f} m³/s, but the overflow peak increases from {am['wetland_overflow_total_peak_m3_s']:.6f} to {bm['wetland_overflow_total_peak_m3_s']:.6f} m³/s. Consequently the combined wetland-to-river peak increases from {am['wetland_to_river_total_peak_m3_s']:.6f} to {bm['wetland_to_river_total_peak_m3_s']:.6f} m³/s. The stronger synthetic resistance/ET package therefore does not automatically reduce the event river maximum; the overflow pathway offsets the reduction in normal drainage.",'',
                f"The combined wetland outflow peak time shifts from {am['wetland_to_river_total_peak_time_h']:.4f} to {bm['wetland_to_river_total_peak_time_h']:.4f} h, while its centroid shifts from {am['wetland_to_river_total_centroid_time_h']:.4f} to {bm['wetland_to_river_total_centroid_time_h']:.4f} h, and its total event volume changes from {am['wetland_to_river_total_volume_m3']:.2f} to {bm['wetland_to_river_total_volume_m3']:.2f} m³. Nevertheless, the river maximum changes by {(bm['river_peak_m']-am['river_peak_m'])*1000:+.4f} mm and its time shifts from {am['river_peak_time_h']:.4f} to {bm['river_peak_time_h']:.4f} h. A later event-volume centroid and a smaller integrated volume do not imply a smaller or later instantaneous river peak. These are model-specific class-package responses, not a general species ranking.",'']
    skipped=controlled.get('skipped_duplicate_blocks',[])
    if skipped:
        lines += [f"Duplicate controlled blocks were omitted because their gates and covers matched an already evaluated block: `{_json(skipped)}`. Their absence does not represent an unreported comparison.",'']
    lines += ['The class package changes the existing resistance and ET endpoints together, with consequent wetland storage, drainage, overflow and groundwater responses. It does not identify a separate causal effect of resistance alone, a real species response or a Milwaukee vegetation recommendation. See [controlled_comparisons.csv](../tables/no_spillway_revision/controlled_comparisons.csv) for all outcomes and diagnostics.','',
        '## Constraints and numerical resolution','',
        'The unchanged optimizer receives 20 scaled residuals, with eight further rows independently auditing native/reconstructed gate and cover bounds. Water-level and G5 residuals are physical level differences divided by a fixed 1 m scale. Nonnegative storages use the eight frozen positive storage scales; soil capacity remains `(S_m−120000 m³)/90000 m³`. No absolute elevation is used as a divisor, and large residuals are not clipped. The audit tolerance of 1e-6 corresponds to 1e-6 m for level rows and 0.09 m³ for soil capacity; strict runtime domain checks may be stronger.','',
        'The callback objective remains `(physical H_peak−97 m)/1 m`; the separate penalty is `100 × [sum(g_i for g_i>1e-6)]²`. The relative objective tolerance `tol_f=1e-6` has a nominal feasible-objective spread of about 0.00056 mm near the current river peak. That stopping scale is much smaller than a millimeter but does not supply numerical integration accuracy or establish convergence. Physical peaks and feasible designs, not penalty scores, are reported here.','',
        'Safety levels remain 101.20/101.00 m for pools and 98.90/98.80 m for contained wetland units. Physical residual ≤0 denotes satisfaction; scales and physical-unit tolerances are recorded in [constraint_residuals.csv](../tables/no_spillway_revision/constraint_residuals.csv). G5 gate applicability is a strict integration-stage restriction; finite infeasible penalties are not physical river peaks. The inherited A.2 state checks already reject any negative storage or soil storage above capacity at every RHS/RK evaluation. G5 and hard geometry limits also execute at those stages. The four design safety limits and the extra 0.001 m geometry margin are audited at reporting nodes, including final time. This corrects the earlier report’s overly broad statement that ordinary storage constraints were reporting-node-only; historical files are preserved. These finite checks do not establish continuous-time safety or site certification.','',
        '| Fixed design | Step (s) | Physical peak (m) | Feasible | Max balance residual (m³) | Relative balance residual | Execution |',
        '|---|---:|---:|---|---:|---:|---|']
    for r in resolution['records']:
        metric=r.get('metrics',{})
        lines.append(f"| {r['label']} | {r['dt_s']} | {_fmt(metric.get('peak_stage_m'),9)} | {r['feasible']} | {_fmt(metric.get('balance_max_m3'),6)} | {metric.get('balance_relative','—')} | {'reused: '+str(r['simulation_reused_from']) if r.get('simulation_reused_from') else 'ODE evaluated'} |")
    lines += ['', '`historical_design_recomputed` means the previous saved decision vector was simulated again under this new no-spillway model. It is not a historical output relabeled as current, and an unsafe recomputed design is not an admissible performance comparator. A stable grid result for an unsafe design establishes neither feasibility nor a usable improvement.','']
    valid=[r for r in resolution['records'] if r.get('status')=='valid']
    if valid:
        rain=np.asarray([r['rain_depth_m'] for r in valid]);solar=np.asarray([r['solar_energy_Wh_m2'] for r in valid])
        lines += ['',f"Fixed-design input checks: rainfall {rain[0]*1000:.6f} mm (spread {np.ptp(rain):.3g} m), solar energy {solar[0]:.6f} Wh/m² (spread {np.ptp(solar):.3g} Wh/m²). The same underlying prepared forcing knots and gate functions are used at each step, with no reoptimization.",'']
    lines += ['| Difference / grid comparison | Coarser difference (mm) | Finer difference (mm) | Summed grid changes (mm) | Resolution screen |',
              '|---|---:|---:|---:|---|']
    for r in resolution['comparisons']:
        fine30=r.get('refinement')=='60_vs_30';coarse_key='difference_60s_m' if fine30 else 'difference_120s_m';fine_key='difference_30s_m' if fine30 else 'difference_60s_m'
        grid_label='60/30 s' if fine30 else '120/60 s'
        if r.get('identical_design'):grid_label+='; identical design, exact zero comparison'
        lines.append(f"| {r['second']} − {r['first']} / {grid_label} | {_fmt(None if r.get(coarse_key) is None else r[coarse_key]*1000,6)} | {_fmt(None if r.get(fine_key) is None else r[fine_key]*1000,6)} | {_fmt(None if r.get('sum_grid_changes_m') is None else r['sum_grid_changes_m']*1000,6)} | {'Stable at tested steps' if r.get('resolved') else 'Unresolved'} |")
    lines += ['',f"Thirty-second check requested by the recorded screen: {resolution.get('thirty_second_required',False)}. The comparison is an empirical timestep check, not a rigorous error bound. If a reported difference is comparable with discretization discrepancies or changes rank, its millimeter-scale magnitude/ranking is unresolved. The controlled-class table uses its declared {controlled.get('dt_s','unrecorded')}-second grid; near-tied class rankings are not automatically validated by a timestep check of different joint designs.",'']
    tests_path=output/'tests.log'
    if tests_path.exists() and 'Ran 77 tests' in tests_path.read_text() and '\nOK' in tests_path.read_text():
        lines += ['The recorded test run passed **77 tests**, including **10 structural tests** for this revision: gate-only pool routing, retained wetland overflow, original units/output schema, absence of safety clipping and internal/external conservation checks. These software tests verify implementation properties; they do not calibrate or physically certify the model. The original A.2 checks and historical B.2 tests are included in the preserved suite.','']
    if reference is not None:
        lines += ['The accepted new reference uses covers 0.65/0.60 and the existing gate knots. Its full opening vector is recorded in `reference_search.json` and the optimization CSV; the non-storm value is 0.005, with larger releases only at existing storm knots. The exact schedule is reproduced below.','',
                  '| Gate | Existing knot hours | Opening fractions |','|---|---|---|']
        n1=len(config['controls']['gate1_hours'])
        for i,(start,end) in enumerate([(0,n1),(n1,len(reference['theta'])-2)]):
            lines.append(f"| {i+1} | `{_json(config['controls'][f'gate{i+1}_hours'])}` | `{_json(reference['theta'][start:end])}` |")
        lines += ['']
    lines += ['## Reproduction','',
        'Use the existing project Python environment and pinned official `multivarious` checkout. From the project root:','',
        '```sh','python3 -B -m B.work.no_spillway_reference_search','python3 -B -m B.no_spillway_study --screen',
        'python3 -B -m B.no_spillway_study --run','python3 -B -m B.no_spillway_study --compare',
        'python3 -B -m B.no_spillway_study --report','```','',
        'Completed execution stages refuse to overwrite their existing results. For this completed case, the `--report` command regenerates reports only and does not run the ODE or optimizer. Preserve the recorded configuration lock and source checksums; a fresh numerical study requires a separate output location rather than deleting or relabeling this one.','']
    for fig in figures:
        lines += [f"![{fig.stem.replace('_',' ')}](../plots/no_spillway_revision/{fig.name})",'']
    lines += ['All new artifacts belong to the no-spillway revision folders. Historical A.2/B.2 files and results, approved event-data provenance and the official solver implementation are preserved. No Git commit or publication is part of this study.','']
    document=ROOT/'docs/11_no_spillway_revision.md';document.write_text('\n'.join(lines))
    return dict(document=str(document),tables=[str(p) for p in tables],figures=[str(p) for p in figures])
