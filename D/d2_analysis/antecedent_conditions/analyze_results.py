"""Saved-trajectory comparisons; no ODE execution or parameter sampling."""
from __future__ import annotations
import json
import os
from pathlib import Path
import sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/cee690-antecedent-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from D.feedback_control import run_policy as frozen
from D.d2_analysis.parameter_uncertainty import run_study as prior
from D.d2_analysis.analysis import exceedance
from D.d2_analysis.antecedent_conditions import preservation

TABLES=HERE/'tables';PLOTS=HERE/'plots'
LABELS={'historical_20190825T0200Z':'August 2019','historical_20240804T1500Z':'August 2024'}
COLORS={0:'#555555',24:'#16867a',72:'#cb663a'}


def write(path,value):frozen.write_json(path,value)


def load(plan,event,point,hours,case):
    name=prior.stem(point['sample_id'],case,15.)
    if hours==0:
        folder=ROOT/'D/d2_analysis/historical_storm_robustness/results/pilot'/event['event_id']
    else:folder=HERE/'results'/f"{event['event_id']}__{hours}h__{point['sample_id']}"
    path=folder/(name+'.json');d=json.loads(path.read_text())
    result=dict(case_name=case,diagnostics=d,valid=bool(d.get('full_valid_trajectory_available')),source=str(path.relative_to(ROOT)),source_sha256=frozen.sha(path))
    if result['valid']:
        with np.load(folder/(name+'.npz'),allow_pickle=False) as z:result.update({k:z[k].copy() for k in z.files})
        result['controller_dataframe']=pd.read_csv(folder/(name+'_controller.csv.gz'),float_precision='round_trip')
    return result


def window_metrics(r,event,point,hours,scope,left,right):
    d=r['diagnostics'];row=dict(event_id=event['event_id'],sample_id=point['sample_id'],antecedent_hours=hours,
        case_name=r['case_name'],scope=scope,valid=r['valid'],classification=d['classification'],
        combined_horizon_unsafe=d.get('safety_violation'),source=r['source'],source_sha256=r['source_sha256'],
        window_start_s=left,window_end_s=right,window_duration_h=(right-left)/3600,
        failure_reason=d.get('failure_reason') or d.get('failure_message'))
    if not r['valid']:return row,[],[]
    t=r['time'];ix=(t>=left)&(t<=right);tt=t[ix];h=r['downstream_stage'][ix];q=r['downstream_discharge'][ix]
    x=r['x'][:,ix]; bal=r['balance'][ix]-r['balance'][np.flatnonzero(ix)[0]]
    volume=float(np.trapezoid(q,tt));net=float(x[:,-1].sum()-x[:,0].sum()-bal[-1]);start=pd.Timestamp(event['start_utc'])-pd.Timedelta(hours=hours)
    ph=int(np.argmax(h));pq=int(np.argmax(q));dur,area=exceedance(tt,h,h[0])
    row.update(window_start_utc=(start+pd.Timedelta(seconds=left)).isoformat(),window_end_utc=(start+pd.Timedelta(seconds=right)).isoformat(),peak_stage_m=float(h[ph]),peak_discharge_m3_s=float(q[pq]),initial_stage_m=float(h[0]),final_stage_m=float(h[-1]),
        peak_stage_time_from_window_start_h=float((tt[ph]-left)/3600),peak_stage_time_from_event_start_h=float(tt[ph]/3600-hours),
        peak_discharge_time_from_event_start_h=float(tt[pq]/3600-hours),peak_stage_utc=(start+pd.Timedelta(seconds=float(tt[ph]))).isoformat(),
        peak_at_window_initial=bool(ph==0),peak_discharge_at_window_initial=bool(pq==0),maximum_rise_above_window_initial_m=float(h.max()-h[0]),
        duration_above_own_window_initial_h=dur/3600,initial_reference_area_m_h=area/3600,
        initial_total_storage_m3=float(x[:,0].sum()),final_total_storage_m3=float(x[:,-1].sum()),
        storage_change_m3=float(x[:,-1].sum()-x[:,0].sum()),outlet_volume_m3=volume,
        external_net_volume_m3=net,nonoutlet_external_net_volume_m3=net+volume,
        end_mass_residual_m3=float(bal[-1]),max_mass_residual_m3=float(np.abs(bal).max()),
        accounting_method='External net recovered from unchanged block-integrated mass residual; Q_out is continuous and integrated over saved nodes',
        flood_threshold_m=None,flood_metrics_status='Undefined without an accepted flood threshold; window-initial response is separate')
    phase=d.get('phase_diagnostics',{}).get('event' if scope=='event_window' else scope)
    if phase:
        row.update(direct_block_outlet_volume_m3=phase['outlet_volume_m3'],direct_block_external_net_volume_m3=phase['external_net_volume_m3'],
            own_phase_balance_relative=phase['relative_balance_residual'],own_phase_balance_passed=phase['balance_passed'],
            phase_forcing_conservation_passed=phase['forcing_conservation_passed'],phase_rhs_unsafe=phase['safety_exceeded_at_successful_rhs'],
            endpoint_vs_direct_external_net_error_m3=net-phase['external_net_volume_m3'],joined_vs_direct_outlet_error_m3=volume-phase['outlet_volume_m3'])
    elif scope=='full_horizon':
        row.update(own_phase_balance_relative=d.get('balance_relative'),own_phase_balance_passed=d.get('balance_relative',float('inf'))<d.get('balance_relative_tolerance',0),phase_rhs_unsafe=d.get('safety_violation_any_successful_rhs'))
    names=list(r['output_names']);safety=[]
    for variable,limit in [('H_d1',101.2),('H_d2',101.),('H_w1',98.9),('H_w2',98.8)]:
        if variable not in names:continue
        values=r['y'][names.index(variable),ix];k=int(np.argmax(values))
        safety.append({k:row[k] for k in ('event_id','sample_id','antecedent_hours','case_name','scope')}|dict(variable=variable,threshold_m=limit,
            node_peak_m=float(values[k]),node_margin_m=float(limit-values[k]),exceeded=bool(values[k]>limit),peak_time_from_event_h=float(tt[k]/3600-hours),
            evidence='Accepted output nodes; all successful RHS safety is separately retained in run/phase diagnostics'))
    row['window_node_unsafe']=any(v['exceeded'] for v in safety)
    states=[]
    for i,name in enumerate(r['state_names']):
        states.append({k:row[k] for k in ('event_id','sample_id','antecedent_hours','case_name','scope')}|dict(state_name=str(name),initial_m3=float(x[i,0]),final_m3=float(x[i,-1]),change_m3=float(x[i,-1]-x[i,0])))
    return row,safety,states


def main():
    TABLES.mkdir(exist_ok=True);PLOTS.mkdir(exist_ok=True)
    plan=json.loads((HERE/'study_plan.json').read_text());data={};rows=[];safety=[];states=[];controllers=[];rhs_rows=[];execution=[]
    for event in plan['events']:
        eid=event['event_id']
        for point in plan['points']:
            sid=point['sample_id']
            for hours in (0,24,72):
                for case in prior.CASES:
                    r=load(plan,event,point,hours,case);data[(eid,sid,hours,case)]=r
                    d=r['diagnostics']
                    execution.append(dict(event_id=eid,sample_id=sid,antecedent_hours=hours,case_name=case,reused_baseline=hours==0,valid=r['valid'],classification=d['classification'],completed=d.get('full_integration_horizon_completed'),unsafe=d.get('safety_violation'),full_balance_relative=d.get('balance_relative'),hydraulic_regimes=json.dumps(d.get('hydraulic_diagnostics',{}).get('regime_rhs_counts',{})),failure_reason=d.get('failure_reason') or d.get('failure_message')))
                    checks={'full_horizon':d['all_rhs_observations']} | {k:v['rhs_observations'] for k,v in d.get('phase_diagnostics',{}).items()}
                    for scope,obs in checks.items():
                        for safe in obs.get('safety_successful_rhs',[]):
                            rhs_rows.append(dict(event_id=eid,sample_id=sid,antecedent_hours=hours,case_name=case,scope=scope,**safe))
                    windows=[('full_horizon',0,(144+hours)*3600),('event_window',hours*3600,(hours+144)*3600)]
                    if hours:windows.append(('antecedent',0,hours*3600))
                    for scope,left,right in windows:
                        row,ss,xx=window_metrics(r,event,point,hours,scope,left,right);rows.append(row);safety.extend(ss);states.extend(xx)
                    if r['valid']:
                        ctrl=r['controller_dataframe'];boundary=hours*3600
                        for branch in (1,2):
                            c=ctrl[ctrl.branch==branch];at=c[c.time_s==boundary].iloc[0]
                            controllers.append(dict(event_id=eid,sample_id=sid,antecedent_hours=hours,case_name=case,branch=branch,
                                initial_opening_m=float(c.start_opening_m.iloc[0]),event_start_opening_m=float(at.start_opening_m),
                                event_start_enabled=bool(at.enabled),final_opening_m=float(c.attained_end_opening_m.iloc[-1]),
                                maximum_opening_m=float(max(c.start_opening_m.max(),c.attained_end_opening_m.max())),
                                maximum_abs_rate_m_s=float(c.rate_m_per_s.abs().max()),latch_changes=int(c.latch_changed.sum()),
                                completed_blocks=int(c.completed_block.sum()),controller_steps=len(c)))
    df=pd.DataFrame(rows);df.to_csv(TABLES/'model_window_metrics.csv',index=False)
    pd.DataFrame(safety).to_csv(TABLES/'safety_at_output_nodes.csv',index=False)
    pd.DataFrame(states).to_csv(TABLES/'compartment_storages.csv',index=False)
    pd.DataFrame(controllers).to_csv(TABLES/'controller_summary.csv',index=False)
    pd.DataFrame(rhs_rows).to_csv(TABLES/'safety_at_successful_rhs.csv',index=False)
    pd.DataFrame(execution).to_csv(TABLES/'execution_accounting.csv',index=False)
    pairs=[]
    for (eid,sid,hours,scope),group in df.groupby(['event_id','sample_id','antecedent_hours','scope'],sort=False):
        g=group.set_index('case_name');w=g.loc['wetland'];n=g.loc['no_wetland'];valid=bool(w.valid and n.valid)
        a,b=(data[(eid,sid,int(hours),case)] for case in prior.CASES)
        clock=weather=common=False
        if valid:
            clock=bool(np.array_equal(a['time'],b['time']))
            weather=bool(clock and np.array_equal(a['u'][2:],b['u'][2:]))
            ar=a['diagnostics']['parameter_evaluation'];br=b['diagnostics']['parameter_evaluation']
            common=all(ar[k]==br[k] for k in ('effective_physical_values','policy_sha256','prepared_forcing_sha256','boundaries','gate_parameters','design_categories','design_covers'))
            valid=bool(clock and weather and common)
        row=dict(event_id=eid,sample_id=sid,antecedent_hours=int(hours),scope=scope,paired_valid=valid,paired_clocks_exact=clock,paired_weather_exact=weather,paired_physical_provenance_equal=common)
        if valid:
            for col in ['peak_stage_m','peak_discharge_m3_s','initial_stage_m','maximum_rise_above_window_initial_m',
                        'duration_above_own_window_initial_h','outlet_volume_m3','initial_total_storage_m3','final_total_storage_m3',
                        'external_net_volume_m3','nonoutlet_external_net_volume_m3','end_mass_residual_m3']:
                row['delta_'+col]=float(n[col]-w[col])
            row.update(wetland_peak_stage_m=float(w.peak_stage_m),no_wetland_peak_stage_m=float(n.peak_stage_m),
                wetland_peak_time_h=float(w.peak_stage_time_from_event_start_h),no_wetland_peak_time_h=float(n.peak_stage_time_from_event_start_h),
                both_initial_peaks=bool(w.peak_at_window_initial and n.peak_at_window_initial),
                either_node_unsafe=bool(w.window_node_unsafe or n.window_node_unsafe))
            prediction=(row['delta_initial_total_storage_m3']-row['delta_final_total_storage_m3']+
                row['delta_nonoutlet_external_net_volume_m3']+row['delta_end_mass_residual_m3'])
            row['paired_outlet_identity_error_m3']=float(row['delta_outlet_volume_m3']-prediction)
            a,b=(data[(eid,sid,int(hours),case)] for case in prior.CASES)
            row['weather_arrays_exact']=bool(np.array_equal(a['u'][2:],b['u'][2:]))
            row['full_horizon_max_paired_gate_difference_m']=float(np.max(np.abs((a['u'][:2]-b['u'][:2])*.3)))
            row['full_horizon_paired_latch_mismatch_count']=int(np.count_nonzero(a['controller_dataframe'].enabled.to_numpy()!=b['controller_dataframe'].enabled.to_numpy()))
        pairs.append(row)
    pair=pd.DataFrame(pairs);pair.to_csv(TABLES/'paired_window_metrics.csv',index=False)
    small=pair[pair.scope.isin(['full_horizon','event_window'])].pivot(index=['event_id','sample_id'],columns=['scope','antecedent_hours'],values='delta_peak_stage_m')*1000
    small.columns=[f'{sc}__{h}h__reduction_mm' for sc,h in small.columns]
    small.to_csv(TABLES/'peak_reduction_comparison.csv')
    account_cols=['event_id','sample_id','antecedent_hours','scope','paired_valid','delta_initial_total_storage_m3','delta_final_total_storage_m3','delta_nonoutlet_external_net_volume_m3','delta_end_mass_residual_m3','delta_outlet_volume_m3','paired_outlet_identity_error_m3']
    pair[account_cols].to_csv(TABLES/'paired_water_accounting.csv',index=False)
    summaries=[]
    for (eid,h,sc),g in pair.groupby(['event_id','antecedent_hours','scope'],sort=False):
        summaries.append(dict(event_id=eid,antecedent_hours=h,scope=sc,selected_configurations=len(g),probability_interpretation=False,
            valid_pairs=int(g.paired_valid.sum()),positive_peak_contrasts=int((g.delta_peak_stage_m>0).sum()),zero_peak_contrasts=int((g.delta_peak_stage_m==0).sum()),negative_peak_contrasts=int((g.delta_peak_stage_m<0).sum()),
            minimum_peak_reduction_mm=g.delta_peak_stage_m.min()*1000,mean_peak_reduction_mm=g.delta_peak_stage_m.mean()*1000,maximum_peak_reduction_mm=g.delta_peak_stage_m.max()*1000,mean_outlet_difference_m3=g.delta_outlet_volume_m3.mean()))
    pd.DataFrame(summaries).to_csv(TABLES/'event_conditioning_summary.csv',index=False)
    comparisons=[]
    for event in plan['events']:
        for point in plan['points']:
            subset=pair[(pair.event_id==event['event_id'])&(pair.sample_id==point['sample_id'])&(pair.scope=='event_window')].set_index('antecedent_hours')
            base=subset.loc[0]
            for hours in (24,72):
                r=subset.loc[hours];out=dict(event_id=event['event_id'],sample_id=point['sample_id'],antecedent_hours=hours)
                for col in ['delta_peak_stage_m','delta_peak_discharge_m3_s','delta_outlet_volume_m3','delta_initial_stage_m']:
                    out['baseline_'+col]=base.get(col);out['conditioned_'+col]=r.get(col)
                    if pd.notna(base.get(col)) and pd.notna(r.get(col)):out['change_'+col]=r[col]-base[col]
                comparisons.append(out)
    pd.DataFrame(comparisons).to_csv(TABLES/'change_from_saved_baseline.csv',index=False)
    # All three predeclared vectors are plotted; none selected by response.
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    for event in plan['events']:
        eid=event['event_id'];fig,axes=plt.subplots(2,3,figsize=(14,7),sharex=True)
        for j,point in enumerate(plan['points']):
            sid=point['sample_id']
            for hours in (0,24,72):
                for case in prior.CASES:
                    r=data[(eid,sid,hours,case)]
                    if not r['valid']:continue
                    t=r['time']/3600-hours;style='-' if case=='wetland' else '--'
                    label=f'{hours} h prefix; '+('wetland' if case=='wetland' else 'no wetland')
                    axes[0,j].plot(t,r['downstream_stage'],style,color=COLORS[hours],lw=1.25,label=label)
                    axes[1,j].plot(t,r['downstream_discharge'],style,color=COLORS[hours],lw=1.25)
            axes[0,j].set_title(sid);axes[0,j].set_ylabel('River stage (m)');axes[1,j].set_ylabel('Outlet discharge (m³/s)')
            for ax in axes[:,j]:ax.axvline(0,color='k',lw=.7);ax.axvspan(-72,0,color='grey',alpha=.07);ax.grid(alpha=.18);ax.set_xlim(-72,144)
            axes[1,j].set_xlabel('Hours relative to original event start')
        handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3,bbox_to_anchor=(.5,0))
        fig.suptitle(f"{LABELS[eid]}: continuous antecedent + event trajectories\nOriginal no-prefix result retained; event starts at 0 h; initial peaks included")
        fig.tight_layout(rect=(0,.10,1,.93));fig.savefig(PLOTS/f'01_{eid}_hydrographs.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(13,7))
    for i,event in enumerate(plan['events']):
        eid=event['event_id']
        for sid in [p['sample_id'] for p in plan['points']]:
            s=pair[(pair.event_id==eid)&(pair.sample_id==sid)&(pair.scope=='event_window')].sort_values('antecedent_hours')
            for j,(metric,scale,label) in enumerate([('delta_peak_stage_m',1000,'Event peak reduction (mm)'),('delta_peak_discharge_m3_s',1,'Event peak discharge reduction (m³/s)'),('delta_outlet_volume_m3',1,'Event outlet volume difference (m³)')]):
                axes[i,j].plot(s.antecedent_hours,s[metric]*scale,'o-',label=sid)
                axes[i,j].axhline(0,color='k',lw=.6);axes[i,j].set_xticks([0,24,72]);axes[i,j].set_ylabel(label);axes[i,j].set_xlabel('Antecedent duration (h)');axes[i,j].grid(alpha=.18)
        axes[i,0].set_title(LABELS[eid]);axes[i,2].legend(fontsize=8)
    fig.suptitle('Original 144 h event-window contrasts: no wetland minus wetland\nThree predeclared configurations; not a new probabilistic ensemble')
    fig.tight_layout(rect=(0,0,1,.93));fig.savefig(PLOTS/'02_paired_event_outcomes.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(13,7),sharex=True)
    for i,event in enumerate(plan['events']):
        eid=event['event_id']
        for j,point in enumerate(plan['points']):
            for hours in (0,24,72):
                for case in prior.CASES:
                    r=data[(eid,point['sample_id'],hours,case)]
                    if not r['valid']:continue
                    t=r['time']/3600-hours
                    # Both branches included; paired-case agreement is audited numerically.
                    for branch in (0,1):
                        style='-' if branch==0 else ':'
                        axes[i,j].plot(t,r['u'][branch]*.3,style,color=COLORS[hours],lw=1,alpha=1 if case=='wetland' else .35,
                            label=f'{hours} h; Gate {branch+1}' if case=='wetland' else None)
            axes[i,j].axvline(0,color='k',lw=.7);axes[i,j].set_xticks([-72,-24,0,48,96,144]);axes[i,j].set_title(f"{LABELS[eid]} / {point['sample_id']}");axes[i,j].set_ylabel('Actual gate opening (m)');axes[i,j].set_xlabel('Hours relative to event start');axes[i,j].grid(alpha=.18)
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3)
    fig.suptitle('Unchanged feedback policy, actual gate histories\nWetland and no-wetland traces overlaid; numerical agreement reported separately')
    fig.tight_layout(rect=(0,.10,1,.93));fig.savefig(PLOTS/'03_gate_histories.png',dpi=180);plt.close(fig)
    # Storage contrasts use exact phase initial values, not a fixed -5500 m3 correction.
    fig,axes=plt.subplots(1,2,figsize=(13,5))
    for ax,event in zip(axes,plan['events']):
        subset=pair[(pair.event_id==event['event_id'])&(pair.scope=='event_window')]
        for ci,sid in enumerate([p['sample_id'] for p in plan['points']]):
            color=plt.get_cmap('tab10')(ci)
            s=subset[subset.sample_id==sid].sort_values('antecedent_hours')
            ax.plot(s.antecedent_hours,-s.delta_initial_total_storage_m3,'o-',color=color,label=sid+' event start')
            ax.plot(s.antecedent_hours,-s.delta_final_total_storage_m3,'x--',color=color,label=sid+' final')
        ax.set_title(LABELS[event['event_id']]);ax.set_xticks([0,24,72]);ax.set_xlabel('Antecedent duration (h)');ax.set_ylabel('Extra wetland-system water: wet − no (m³)');ax.grid(alpha=.18);ax.legend(fontsize=7)
    fig.suptitle('Whole-system storage contrast; water is not tagged by origin')
    fig.tight_layout(rect=(0,0,1,.93));fig.savefig(PLOTS/'04_storage_contrasts.png',dpi=180);plt.close(fig)
    summary=dict(original_baseline_pairs_reused=6,new_primary_pairs=12,new_primary_valid_pairs=int(pair[(pair.antecedent_hours>0)&(pair.scope=='full_horizon')].paired_valid.sum()),
        max_paired_balance_identity_error_m3=float(pair.paired_outlet_identity_error_m3.abs().max()),
        any_node_unsafe=bool(pd.DataFrame(safety).exceeded.any()),
        full_horizon_max_paired_gate_difference_m=float(pair.full_horizon_max_paired_gate_difference_m.max()),
        paired_latch_mismatches=int(pair.loc[pair.scope=='full_horizon','full_horizon_paired_latch_mismatch_count'].sum()),
        no_new_ODE_or_sampling_in_analysis=True,original_results_untouched=True,
        interpretation='Three selected saved parameter configurations do not estimate probabilities; event-window starts can differ between models after conditioning.')
    write(TABLES/'analysis_summary.json',summary)
    write(HERE/'preservation_check.json',preservation.verify())
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
