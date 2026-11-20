"""Saved-result system-storage and paired-volume accounting; never runs an ODE.

Outputs are confined to this new robustness_checks directory. No initial state,
model, controller, parameter draw, forcing or historical artifact is altered.
"""
from __future__ import annotations
from pathlib import Path
import hashlib
import json
import os
import sys
sys.dont_write_bytecode = True
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/cee690-storage-checks-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
HIST = ROOT/'D/d2_analysis/historical_storm_robustness'
PARAM = ROOT/'D/d2_analysis/parameter_uncertainty'
CASES = ('wetland', 'no_wetland')
COLORS = ['#0072B2', '#009E73', '#D55E00', '#777777']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def clean(value):
    if isinstance(value, dict): return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value, (list, tuple)): return [clean(v) for v in value]
    if isinstance(value, np.generic): value=value.item()
    if isinstance(value, float) and not np.isfinite(value): return None
    return value


def label(event):
    return pd.Timestamp(event['start_utc']).strftime('%b %Y') + (' *' if event['benchmark'] else '')


def trajectory_path(event, sample_id, case):
    if event['benchmark']:
        path=PARAM/'results/pilot'/f'pilot_{sample_id}__{case}__15s.npz'
    else:
        path=HIST/'results/ensemble'/event['event_id']/f'{sample_id}__{case}__15s.npz'
    return path if path.exists() else None


def last_rain_end(event):
    frame=pd.read_csv(event['forcing_file'])
    if event['benchmark']:
        t=frame.time_s.to_numpy(float);p=frame.P.to_numpy(float)
        active=np.flatnonzero(p>0)
        return None if not len(active) else float(t[min(active[-1]+1,len(t)-1)])
    active=np.flatnonzero(frame.precipitation_mm_per_h.to_numpy(float)>0)
    return None if not len(active) else float((active[-1]+1)*3600)


def integrate_after(t, y, start):
    """Integral of the saved piecewise-linear curve, including an exact start."""
    t=np.asarray(t);y=np.asarray(y);mask=t>start
    return float(np.trapezoid(np.r_[np.interp(start,t,y),y[mask]],np.r_[start,t[mask]]))


def final_entry_below(t, value, threshold):
    """Last downward crossing followed by values <= threshold to the endpoint."""
    if value[-1]>threshold: return None
    above=np.flatnonzero(value>threshold)
    if not len(above): return float(t[0])
    k=int(above[-1]);fraction=(value[k]-threshold)/(value[k]-value[k+1])
    return float(t[k]+fraction*(t[k+1]-t[k]))


def main():
    manifest_path=HIST/'events/manifest.json'
    events=json.loads(manifest_path.read_text())['events']
    plan_path=HIST/'study_plan.json';plan=json.loads(plan_path.read_text())
    expected={p['sample_id']:p for p in plan['points']}
    assert len(expected)==64 and len(events)==5
    assert plan['events_manifest_sha256']==sha(manifest_path)
    tables=HERE/'tables';plots=HERE/'plots';tables.mkdir(exist_ok=True);plots.mkdir(exist_ok=True)
    models=[];pairs=[];timing=[];sources=[];compartments=[]
    worst_closure=worst_outlet=worst_npz_net=0.
    trajectory_pairs=0;unavailable=[]
    for event in events:
        eid=event['event_id'];directory=HIST/'results/ensemble'/eid
        paired_path=directory/'paired_results.csv'
        frame=pd.read_csv(paired_path,float_precision='round_trip')
        assert len(frame)==64 and set(frame.sample_id)==set(expected)
        assert sha(event['forcing_file'])==event['forcing_sha256']
        rain_end=last_rain_end(event)
        sources.append(dict(event_id=eid,paired_results=str(paired_path.relative_to(ROOT)),
            paired_results_sha256=sha(paired_path),forcing_sha256=event['forcing_sha256'],last_rain_interval_end_h=rain_end/3600))
        for row in frame.itertuples(index=False):
            sample=row.sample_id;by_case={};arrays={}
            for case in CASES:
                path=directory/f'{sample}__{case}__15s.json'
                d=json.loads(path.read_text());valid=d['full_valid_trajectory_available']
                assert valid==getattr(row,case+'_valid')
                assert d['parameter_evaluation']['effective_physical_values']==expected[sample]['physical_values']
                rec=dict(event_id=eid,sample_id=sample,case=case,valid=valid,
                    unsafe=d.get('safety_violation'),classification=d['classification'],
                    diagnostic_source=str(path.relative_to(ROOT)),diagnostic_sha256=sha(path))
                if not valid or not d['full_integration_horizon_completed']:
                    rec.update(endpoint_available=False,unavailable_reason=d.get('failure_reason') or 'Incomplete or invalid trajectory')
                    models.append(rec);by_case[case]=rec;continue
                names=list(d['state_names']);x0=np.asarray(d['initial_storage_m3'],float)
                final=d.get('last_completed_state_m3');mass=d.get('completed_blocks_mass_diagnostics',{})
                needed=('external_net_volume_m3','final_balance_residual_m3','through_time_s')
                if final is None or not all(k in mass for k in needed) or mass['through_time_s']!=518400.:
                    rec.update(endpoint_available=False,unavailable_reason='Saved complete final state or complete block balance absent')
                    models.append(rec);by_case[case]=rec;unavailable.append([eid,sample,case]);continue
                xT=np.asarray(final,float);assert len(xT)==len(x0)==len(names)
                assert np.isfinite(xT).all() and d['completed_until_s']==518400.
                s0=float(x0.sum());sT=float(xT.sum());net=float(mass['external_net_volume_m3'])
                residual=float(mass['final_balance_residual_m3'])
                outlet=float(getattr(row,case+'_total_outlet_volume_m3'))
                closure=sT-s0-net-residual;worst_closure=max(worst_closure,abs(closure))
                assert abs(closure)<1e-7,(eid,sample,case,'endpoint balance mismatch',closure)
                rec.update(endpoint_available=True,horizon_s=518400.,initial_total_storage_m3=s0,
                    final_total_storage_m3=sT,system_storage_change_m3=sT-s0,
                    cumulative_outlet_volume_m3=outlet,external_net_integrated_m3=net,
                    net_nonoutlet_external_integral_m3=net+outlet,final_balance_residual_m3=residual,
                    endpoint_identity_error_m3=closure,endpoint_source='Saved complete final-state and block-integral JSON fields',
                    nonoutlet_definition='N = external_net + V_out = integral(P*A_rain + Q_ext + Q_gx - E_total) dt',
                    nonoutlet_method='Recovered accounting quantity; ET and groundwater-boundary exchange are not separated')
                for name,initial,final_value in zip(names,x0,xT):
                    rec['initial_'+name+'_m3']=float(initial);rec['final_'+name+'_m3']=float(final_value)
                tr=trajectory_path(event,sample,case);rec['trajectory_available']=tr is not None
                rec['trajectory_source']=str(tr.relative_to(ROOT)) if tr else None
                if tr:
                    with np.load(tr,allow_pickle=False) as z:
                        a={k:z[k].copy() for k in ('time','x','balance','downstream_discharge','state_names')}
                        if 'external_net_integrated_m3' in z:a['external_net_integrated_m3']=z['external_net_integrated_m3'].copy()
                    arrays[case]=a;t=a['time'];assert t[0]==0 and t[-1]==518400. and np.all(np.diff(t)==15.)
                    assert np.array_equal(a['x'][:,0],x0) and np.array_equal(a['x'][:,-1],xT)
                    assert list(a['state_names'])==names
                    volume=float(np.trapezoid(a['downstream_discharge'],t))
                    worst_outlet=max(worst_outlet,abs(volume-outlet));assert abs(volume-outlet)<1e-7
                    derived_net=sT-s0-float(a['balance'][-1]);worst_npz_net=max(worst_npz_net,abs(derived_net-net))
                    assert abs(derived_net-net)<1e-7 and abs(float(a['balance'][-1])-residual)<1e-7
                    if 'external_net_integrated_m3' in a:assert abs(a['external_net_integrated_m3'][-1]-net)<1e-7
                    rec['cumulative_outlet_after_last_rain_m3']=integrate_after(t,a['downstream_discharge'],rain_end)
                    rec['trajectory_sha256']=sha(tr)
                models.append(rec);by_case[case]=rec
            if not all(by_case[c].get('endpoint_available',False) for c in CASES):
                pairs.append(dict(event_id=eid,sample_id=sample,endpoint_pair_available=False));continue
            w,n=by_case['wetland'],by_case['no_wetland']
            component_final={}
            for state in ('S_u','S_m','S_w1','S_w2','S_d1','S_d2','S_g','S_r'):
                present=state not in ('S_w1','S_w2')
                wet_initial=w['initial_'+state+'_m3'];wet_final=w['final_'+state+'_m3']
                no_initial=n['initial_'+state+'_m3'] if present else None
                no_final=n['final_'+state+'_m3'] if present else None
                if present: assert wet_initial==no_initial, (eid,sample,state,'shared initial stock differs')
                difference=wet_final-(no_final if present else 0.)
                component_final[state]=difference
                compartments.append(dict(event_id=eid,sample_id=sample,state=state,
                    no_wetland_compartment_present=present,initial_wetland_case_m3=wet_initial,
                    initial_no_wetland_case_m3=no_initial,final_wetland_case_m3=wet_final,
                    final_no_wetland_case_m3=no_final,
                    initial_stock_difference_wet_minus_no_m3=wet_initial-(no_initial if present else 0.),
                    final_stock_difference_wet_minus_no_m3=difference,
                    absence_convention='Absent compartment contributes zero to the system sum; it is not a simulated zero state' if not present else None))
            initial_term=n['initial_total_storage_m3']-w['initial_total_storage_m3']
            assert initial_term==-5500.
            final_term=w['final_total_storage_m3']-n['final_total_storage_m3']
            assert abs(sum(component_final.values())-final_term)<1e-7
            nonoutlet_term=n['net_nonoutlet_external_integral_m3']-w['net_nonoutlet_external_integral_m3']
            residual_term=n['final_balance_residual_m3']-w['final_balance_residual_m3']
            observed=n['cumulative_outlet_volume_m3']-w['cumulative_outlet_volume_m3']
            reconstructed=initial_term+final_term+nonoutlet_term+residual_term
            assert abs(observed-reconstructed)<1e-7
            assert abs(observed-row.delta_total_outlet_volume_m3)<1e-7
            pair=dict(event_id=eid,sample_id=sample,endpoint_pair_available=True,
                delta_outlet_no_minus_wet_m3=observed,initial_stock_term_m3=initial_term,
                final_stock_term_m3=final_term,nonoutlet_external_term_m3=nonoutlet_term,
                final_groundwater_stock_difference_m3=component_final['S_g'],
                final_shared_non_groundwater_stock_difference_m3=sum(component_final[k] for k in ('S_u','S_m','S_d1','S_d2','S_r')),
                residual_term_m3=residual_term,reconstructed_delta_outlet_m3=reconstructed,
                decomposition_error_m3=observed-reconstructed,
                wetland_compartments_initial_m3=w['initial_S_w1_m3']+w['initial_S_w2_m3'],
                wetland_compartments_final_m3=w['final_S_w1_m3']+w['final_S_w2_m3'],
                full_trajectory_pair_available=all(c in arrays for c in CASES),
                interpretation='Exact saved finite-window accounting, not a causal partition of initial water, rainfall, ET or wetland benefit')
            if pair['full_trajectory_pair_available']:
                trajectory_pairs+=1
                pair['delta_outlet_after_last_rain_m3']=n['cumulative_outlet_after_last_rain_m3']-w['cumulative_outlet_after_last_rain_m3']
            pairs.append(pair)
            tm=dict(event_id=eid,sample_id=sample,trajectory_available=pair['full_trajectory_pair_available'],
                endpoint_wetland_storage_m3=pair['wetland_compartments_final_m3'],last_rain_interval_end_h=rain_end/3600,
                remaining_window_after_last_rain_h=(518400.-rain_end)/3600)
            if 'wetland' in arrays:
                a=arrays['wetland'];t=a['time'];names=list(a['state_names'])
                sw=a['x'][[names.index('S_w1'),names.index('S_w2')]];total=sw.sum(axis=0)
                k=int(np.argmax(total));at_rain_end=float(np.interp(rain_end,t,total))
                entry=final_entry_below(t,total,55.)
                tm.update(wetland_storage_peak_m3=float(total[k]),wetland_storage_peak_time_h=float(t[k]/3600),
                    storage_at_last_rain_end_m3=at_rain_end,net_wetland_storage_decline_after_last_rain_m3=at_rain_end-float(total[-1]),
                    final_entry_below_55m3_h=None if entry is None else entry/3600,
                    below_55m3_at_end=bool(total[-1]<=55.),
                    threshold_meaning='55 m3 = 1% of the primary 5500 m3 wetland initial stock; descriptive storage threshold, not tagged initial water or a safety threshold',
                    storage_decline_meaning='Net stock decline; cannot separate normal drainage, overflow, seepage, ET or later gate inflow from compact arrays')
                for j in range(2):
                    idx=int(np.argmax(sw[j]));tm[f'wetland_{j+1}_peak_storage_m3']=float(sw[j,idx]);tm[f'wetland_{j+1}_peak_time_h']=float(t[idx]/3600)
            else:tm['unavailable_reason']='No saved full trajectory for this benchmark draw; endpoint state is still directly available in JSON'
            timing.append(tm)
    models=pd.DataFrame(models);pairs=pd.DataFrame(pairs);timing=pd.DataFrame(timing)
    summaries=[]
    for event in events:
        eid=event['event_id'];p=pairs[pairs.event_id.eq(eid)];tm=timing[timing.event_id.eq(eid)];available=tm[tm.trajectory_available.eq(True)]
        s=dict(event_id=eid,label=label(event),endpoint_pairs=len(p),full_trajectory_pairs=len(available),
            positive_outlet_difference_count=int(p.delta_outlet_no_minus_wet_m3.gt(0).sum()),
            zero_outlet_difference_count=int(p.delta_outlet_no_minus_wet_m3.eq(0).sum()),
            negative_outlet_difference_count=int(p.delta_outlet_no_minus_wet_m3.lt(0).sum()))
        for col in ['delta_outlet_no_minus_wet_m3','initial_stock_term_m3','final_stock_term_m3','nonoutlet_external_term_m3','residual_term_m3','wetland_compartments_final_m3','final_groundwater_stock_difference_m3','final_shared_non_groundwater_stock_difference_m3']:
            s['mean_'+col]=float(p[col].mean());s['min_'+col]=float(p[col].min());s['max_'+col]=float(p[col].max())
        if len(available):
            s.update(mean_net_storage_decline_after_last_rain_m3=float(available.net_wetland_storage_decline_after_last_rain_m3.mean()),
                net_storage_decline_after_last_rain_count=int(available.net_wetland_storage_decline_after_last_rain_m3.gt(0).sum()),
                mean_postrain_delta_outlet_m3=float(p.delta_outlet_after_last_rain_m3.mean()),
                more_wetland_outlet_after_last_rain_count=int(p.delta_outlet_after_last_rain_m3.lt(0).sum()),
                final_below_55m3_count=int(available.below_55m3_at_end.sum()),
                median_final_entry_below_55m3_h=float(available.final_entry_below_55m3_h.median()) if available.final_entry_below_55m3_h.notna().any() else None)
        summaries.append(s)
    summaries=pd.DataFrame(summaries)
    compartments=pd.DataFrame(compartments)
    compartments.to_csv(tables/'storage_compartment_differences.csv',index=False)
    compartments.groupby(['event_id','state'],sort=False)['final_stock_difference_wet_minus_no_m3'].agg(['count','mean','min','max']).reset_index().to_csv(tables/'storage_compartment_summary.csv',index=False)
    models.to_csv(tables/'storage_model_balances.csv',index=False)
    pairs.to_csv(tables/'storage_paired_decomposition.csv',index=False)
    timing.to_csv(tables/'storage_wetland_timing.csv',index=False)
    summaries.to_csv(tables/'storage_event_summary.csv',index=False)
    audit=dict(status='pass',physical_model_runs=0,new_initial_condition_experiments=0,
        planned_pairs=320,endpoint_pairs=int(pairs.endpoint_pair_available.sum()),full_trajectory_pairs=trajectory_pairs,
        benchmark_endpoint_pairs=int(pairs.loc[pairs.event_id.eq('benchmark_september_2022'),'endpoint_pair_available'].sum()),
        benchmark_full_trajectory_pairs=int(pairs.loc[pairs.event_id.eq('benchmark_september_2022'),'full_trajectory_pair_available'].sum()),
        unavailable_endpoint_cases=unavailable,maximum_endpoint_identity_error_m3=worst_closure,
        maximum_saved_outlet_integral_difference_m3=worst_outlet,maximum_npz_block_net_difference_m3=worst_npz_net,
        maximum_paired_decomposition_error_m3=float(pairs.decomposition_error_m3.abs().max()),
        maximum_absolute_residual_term_m3=float(pairs.residual_term_m3.abs().max()),
        additional_initial_wetland_water_m3=5500.,sign='All delta terms mean no-wetland minus wetland, except final_stock_term=wetland final minus no-wetland final',
        identity='Delta V_out = Delta S_initial - Delta S_final + Delta N + Delta residual; N=external_net+V_out',
        balance_definition='residual = S_final - S_initial - external_net',
        no_internal_transfer_sum=True,no_water_tagging=True,no_ET_specific_attribution=True,
        benchmark_endpoint_note='All 64 benchmark endpoint pairs are directly recorded in complete diagnostic JSONs, although the previous summary CSV omitted final-storage columns.',
        temporal_scope_note='Timing and post-rain integrals use 256 new-event pairs plus benchmark mc_00000/mc_00001 only; no unsupported benchmark trajectories are reconstructed.',
        data_sources=sources,antecedent_window_audit_sha256=sha(tables/'antecedent_windows.csv'),
        antecedent_continuous_coverage_sha256=sha(tables/'antecedent_continuous_coverage.csv'),
        study_plan_sha256=sha(plan_path),analysis_source_sha256=sha(__file__))
    (tables/'storage_audit.json').write_text(json.dumps(clean(audit),indent=2,allow_nan=False)+'\n')
    make_plot(events,pairs,summaries,plots,audit)
    write_findings(events,summaries,audit)
    print(json.dumps(clean(dict(audit=audit,summary=summaries.to_dict('records'))),indent=2,allow_nan=False))


def make_plot(events,pairs,summary,plots,audit):
    fig=plt.figure(figsize=(14,10),layout='constrained');grid=fig.add_gridspec(2,2,height_ratios=[1.18,1])
    ax=fig.add_subplot(grid[0,:]);low=fig.add_subplot(grid[1,0]);right=fig.add_subplot(grid[1,1])
    y=np.arange(len(events));pos=np.zeros(len(events));neg=pos.copy()
    cols=['initial_stock_term_m3','final_stock_term_m3','nonoutlet_external_term_m3','residual_term_m3']
    names=['Initial stock: no − wet','Final stock: wet − no','Net non-outlet external flux: no − wet','Numerical residual: no − wet']
    for col,name,color in zip(cols,names,COLORS):
        values=summary['mean_'+col].to_numpy();left=np.where(values>=0,pos,neg)
        ax.barh(y,values,left=left,color=color,label=name,height=.57)
        pos+=np.maximum(values,0);neg+=np.minimum(values,0)
    ax.set_xlim(float(neg.min())-400.,float(pos.max())+400.)
    ax.scatter(summary.mean_delta_outlet_no_minus_wet_m3,y,marker='D',c='black',s=34,label='Saved model mean outlet-volume difference',zorder=5)
    for k,v in enumerate(summary.mean_delta_outlet_no_minus_wet_m3):
        ax.annotate(f'{v:+,.1f}',(v,k),xytext=(0,10),textcoords='offset points',ha='center',fontsize=8)
    labels=[label(e) for e in events]
    ax.set(yticks=y,yticklabels=labels,xlabel='Mean accounting terms (m³); no-wetland − wetland outlet volume',title='A. Exact finite-window accounting — not a causal partition of initial water')
    ax.invert_yaxis();ax.axvline(0,color='black',lw=.6);ax.grid(axis='x',alpha=.18)
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.17),ncol=3,fontsize=8)
    groups=[pairs.loc[pairs.event_id.eq(e['event_id'])] for e in events]
    low.boxplot([g.final_stock_term_m3 for g in groups],tick_labels=labels,showfliers=True,widths=.55)
    low.axhline(5500,color='#0072B2',ls='--',lw=1,label='Initial system-storage difference: +5,500 m³')
    low.set(ylabel='Wetland − no-wetland final system storage (m³)',title='B. Final total-system storage difference; 64 pairs/event')
    low.tick_params(axis='x',rotation=18);low.grid(axis='y',alpha=.2);low.legend(fontsize=8,loc='best')
    right.boxplot([g.wetland_compartments_final_m3 for g in groups],tick_labels=labels,showfliers=True,widths=.55)
    assert all((g.wetland_compartments_final_m3>0).all() for g in groups), 'Log plot requires the actual positive recorded endpoint stocks'
    right.set_yscale('log')
    right.set(ylabel='Final S_w1 + S_w2 (m³; log scale)',title='C. Wetland-compartment endpoint; initially 5,500 m³')
    right.tick_params(axis='x',rotation=18);right.grid(axis='y',alpha=.2)
    fig.suptitle('Initial water, endpoint storage and cumulative outlet volume',fontsize=15)
    fig.get_layout_engine().set(h_pad=.12,w_pad=.08,rect=(0,.065,1,.935))
    fig.text(.5,.008,'* Original benchmark endpoints are directly saved for all 64 draws. Boxes show median/IQR and 1.5×IQR whiskers.\nOnly bookkeeping is decomposed; stored water is not tagged by origin. All horizons remain 144 h.',ha='center',fontsize=9)
    fig.savefig(plots/'03_storage_balance.png',dpi=180,bbox_inches='tight');plt.close(fig)


def write_findings(events,summary,audit):
    rows=[]
    for r in summary.itertuples(index=False):
        rows.append(f'| {r.label} | {r.mean_initial_stock_term_m3:+,.2f} | {r.mean_final_stock_term_m3:+,.2f} | {r.mean_nonoutlet_external_term_m3:+,.2f} | {r.mean_residual_term_m3:+.6f} | {r.mean_delta_outlet_no_minus_wet_m3:+,.2f} |')
    components=[]
    for r in summary.itertuples(index=False):
        components.append(f'| {r.label} | {r.mean_final_groundwater_stock_difference_m3:,.2f} | {r.mean_final_shared_non_groundwater_stock_difference_m3:,.2f} | {r.mean_wetland_compartments_final_m3:.6g} | {r.mean_final_stock_term_m3:,.2f} |')
    temporal=[]
    for r in summary.itertuples(index=False):
        temporal.append(f'| {r.label} | {r.full_trajectory_pairs} | {r.mean_wetland_compartments_final_m3:.6g} | {r.mean_net_storage_decline_after_last_rain_m3:,.2f} | {r.more_wetland_outlet_after_last_rain_count}/{r.full_trajectory_pairs} |')
    text='''# Saved-result storage and initial-water accounting

This check reads existing results only. No ODE, new initial-condition experiment, parameter sampling, forcing adjustment or optimization was performed. All primary comparisons retain the accepted extra **5,500 m³** in the wetland compartments and unchanged shared initial states. The analysis does not identify the trajectory of those initial water molecules.

## Coverage and source distinction

All **320 paired endpoint balances** are available. In particular, all 64 benchmark pairs retain `initial_storage_m3`, `last_completed_state_m3` and complete `completed_blocks_mass_diagnostics` in their diagnostic JSONs. The earlier historical summary CSV left final storage blank; that omission does not mean the underlying endpoint data are absent. This new analysis reads those recorded endpoints directly and does not rewrite the previous tables.

Full saved waveforms are available for **258 pairs**: all 256 additional-event pairs and benchmark `mc_00000`/`mc_00001` from the original parameter pilot. Benchmark timing summaries based on two predetermined retained draws are descriptive examples, not estimates representative of all 64 draws. No other benchmark waveform is reconstructed.

## Paired balance identity

For each case, let S0 and ST denote the sum of all water-storage states at the beginning and end; V the integrated river outlet discharge; N the net non-outlet external flux integral; and r the saved numerical balance residual. The unchanged model identity is

```
ST - S0 = N - V + r
N = integral(P*A_rain + Q_ext + Q_gx - E_total) dt
  = saved_external_net + V

Delta V = Delta S0 - Delta ST + Delta N + Delta r
Delta = no-wetland minus wetland
Delta S0 = -5,500 m3
```

The positive sign of Delta r follows from the stored definition `r=ST-S0-external_net`. Internal gate releases, wetland drainage/overflow, soil transfers and other transfers between compartments are not added to external flow. N is recovered from already recorded total external-net accounting and outlet volume; it is not an independent flux measurement or a second conservation test. Its paired contrast combines groundwater-boundary exchange and evapotranspiration differences after common rain/external inflow cancel. These components cannot be separated from the compact arrays used here, so no ET-only attribution is made.

Mean values across 64 pairs per fixed event are:

| Event | Initial term (m³) | Final term, wet−no (m³) | Net non-outlet term (m³) | Residual term (m³) | Outlet no−wet (m³) |
|---|---:|---:|---:|---:|---:|
'''+ '\n'.join(rows)+f'''

The largest absolute paired-identity discrepancy is **{audit['maximum_paired_decomposition_error_m3']:.3g} m³**. The largest case endpoint-balance discrepancy is **{audit['maximum_endpoint_identity_error_m3']:.3g} m³**. Recorded outlet integrals match available saved discharge arrays within **{audit['maximum_saved_outlet_integral_difference_m3']:.3g} m³**. These small arithmetic discrepancies check consistency of saved quantities; numerical timestep validation remains the separately documented original pilot.

The explicit initial-stock term is the same −5,500 m³ in every pair. Final storage and net external fluxes respond dynamically and differ by case. Therefore the four-new-event negative outlet-volume differences cannot be assigned entirely to release of the initial 5,500 m³, nor can the benchmark positive difference be called permanent storm-water removal. No exact causal percentage of a volume contrast is identified by rearranging the balance.

## Where the final system-storage contrast resides

The complete saved endpoint state vectors support a compartmentwise comparison. Mean final stocks below use wetland minus no-wetland; the no-wetland model has no wetland compartments, so their absent contribution to the total sum is zero, without inventing simulated zero-valued states.

| Event | Groundwater difference (m³) | Other shared-state difference (m³) | Wetland-unit stock (m³) | Total system difference (m³) |
|---|---:|---:|---:|---:|
'''+ '\n'.join(components)+'''

The final total-system contrast is concentrated in groundwater, while the additional-event wetland compartments are numerically near empty. Final system storage therefore must not be described as water still held in the wetland units. This is a state-location accounting result, not tagged-water tracking: the final groundwater contrast divided by 5,500 is **not** a fraction of original wetland water remaining. Rainfall, seepage, external exchange and evapotranspiration act together during the event. The full [compartment table](tables/storage_compartment_differences.csv) preserves every pair, shared state and absent-compartment flag; its [summary](tables/storage_compartment_summary.csv) retains ranges.

## Wetland stock and late recession

The last positive-rain interval is defined by the saved forcing, not an outcome-selected hydrograph threshold. The new-event hourly rain is held through the end of its last positive hour; the benchmark retains its original linear 120 s forcing, including the final positive-to-zero segment. Stock is interpolated on the saved 15 s state grid. Net stock decline after rain ends is `S_w(t_rain_end)-S_w(T)`; it is not gross wetland drainage, because gate inflow, drainage, overflow, seepage and evaporation can continue together.

| Event | Timing pairs | Mean endpoint wetland stock across 64 (m³) | Mean post-rain net stock decline, timing subset (m³) | Timing pairs with more post-rain wetland-case outlet |
|---|---:|---:|---:|---:|
'''+ '\n'.join(temporal)+'''

`storage_wetland_timing.csv` also retains the combined and branch storage peaks/times and, where achieved, the final downward crossing into combined stock≤55 m³. The 55 m³ mark is a declared descriptive threshold equal to 1% of the primary initial wetland stock. It is neither a safety threshold nor evidence that 99% of the tagged initial water has left. Very small endpoint volumes reflect the numerical model's near-empty asymptote, not measured water volumes at the displayed precision. The terminal storage and post-rain volume patterns describe late recession; they do not reconstruct unrecorded normal-drainage or overflow hydrographs.

## Proposed alternative initial-condition experiment — not executed

A defensible alternative is **common antecedent-weather conditioning**, preserving water through the unchanged equations. The existing [coverage audit](tables/antecedent_continuous_coverage.csv) finds 351, 93, 0, 254 and 109 consecutive supported nonfreezing hours immediately before the benchmark, August 2023, July 2022, August 2019 and August 2024 event starts, respectively. July 2022 has a missing rainfall hour immediately before its start; missing rain must not be treated as zero. The usable record combines station-derived/interpolated weather, grid-estimated solar and retained trace/QC flags; it is not uniformly direct observation.

A small specific proposed check is **24 versus 72 hours of common antecedent conditioning for August 2019 and August 2024**, using only the existing three pilot vectors `mc_00016`, `mc_00024` and `mc_00018`. This would require 12 pairs / 24 model trajectories including conditioning and the original event, with no new parameter draw. Both durations are supported by the existing [window audit](tables/antecedent_windows.csv). They are exploratory durations, not assumed equilibration times. This proposal is not executed, and it does not silently extend to a 7- or 14-day antecedent record. The primary comparisons remain unchanged.

At the start of that antecedent window, retain the accepted physical storage values for each model, including the wetland model's 5,500 m³ extra water. Apply the same previously saved parameter vector, hypothetical geometry, boundary conditions and feedback policy to both models throughout the common antecedent weather. Carry every resulting water state and the actual gate positions/latch memory continuously into the original 144-hour event; do not reset river level, soil storage, groundwater, wetland water or actuators at the event boundary. No water is deleted, moved arbitrarily between compartments or inserted to equalize totals. The 120 s controller clock must remain aligned through the boundary, and the event forcing remains the saved unchanged input.

This would produce explicit, model-conditioned alternative event-start states. Initial values at the antecedent-window start would still be assumptions, and finite conditioning is not evidence of equilibrium, spin-up convergence or observed antecedent reconstruction. All initial applicability checks, water balances, safety outcomes and numerical criteria remain unchanged. Invalid or unsafe conditioning/event runs must remain reported. Available longer durations could support a separately approved initial-memory sensitivity check; no duration or feasibility outcome is assumed here.

The two models can develop different shared event-start storages and actual gate positions under the same antecedent weather. This experiment therefore tests sensitivity to antecedent state and the persistence of the original initial-water difference; it does **not** isolate a pure causal effect of exactly 5,500 m³ or create matched event-start states by construction. Event metrics must include each resulting event-start stage and retain any initial peak. No such conditioning or alternative trajectory is executed by this saved-result analysis.

## Files

- [System balances](tables/storage_model_balances.csv): all 640 model endpoints and source checksums; full-waveform availability is explicit.
- [Paired decomposition](tables/storage_paired_decomposition.csv): all 320 outlet-volume identities.
- [Wetland timing](tables/storage_wetland_timing.csv): full-waveform diagnostics, with absent benchmark timings left blank.
- [Event summaries](tables/storage_event_summary.csv), [audit](tables/storage_audit.json), and [figure](plots/03_storage_balance.png).

Figure 03 shows the mean accounting terms, final system-storage contrasts and wetland-compartment endpoint stocks. It is an accounting illustration, not causal source-water tracing or empirical engineering validation.
'''
    (HERE/'storage_findings.md').write_text(text)


if __name__=='__main__': main()
