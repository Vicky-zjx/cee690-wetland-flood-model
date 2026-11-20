"""Analyze saved paired parameter results for separately selected historical storms.

No sampling, ODE solving, fitting, optimization or historical-file mutation.
The same 64 parameter vectors are matched by sample_id across fixed events.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import warnings

sys.dont_write_bytecode = True
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/cee690-historical-storm-plots-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import spearmanr, beta, binomtest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from D.d2_analysis.parameter_uncertainty.official_plots import _official_functions, OFFICIAL_REVISION, OFFICIAL_HASHES

SOURCE = HERE.parent/'parameter_uncertainty'
FACTORS = ['m_K_u', 'm_routing', 'm_R_r']
FACTOR_LABELS = ['Infiltration multiplier', 'Routing-rate multiplier', 'River-resistance multiplier']
CASES = ['wetland', 'no_wetland']
METRICS = {
    'peak_stage_m':'m', 'peak_discharge_m3_s':'m3/s',
    'maximum_rise_above_initial_m':'m', 'initial_level_exceedance_h':'h',
    'initial_discharge_excess_m3':'m3', 'total_outlet_volume_m3':'m3',
}
LIMITS = {'H_d1':101.2, 'H_d2':101.0, 'H_w1':98.9, 'H_w2':98.8}
COLORS = ['#242424', '#0072B2', '#D55E00', '#009E73', '#CC79A7']
# Fixed before completing the ensemble; these are CDF evaluation locations,
# never flood thresholds or observations selected to make a result favorable.
CDF_GRIDS = {
    'stage':np.linspace(97.,101.,2001),
    'discharge':np.linspace(0.,3.,601),
    'delta_stage':np.linspace(-1.,1.,2001),
    'delta_discharge':np.linspace(-3.,3.,1201),
}
CP_DOC_URL='https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats._result_classes.BinomTestResult.proportion_ci.html'


def clean(x):
    if isinstance(x, dict):return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x, (list,tuple)):return [clean(v) for v in x]
    if isinstance(x, np.generic):x=x.item()
    if isinstance(x, float) and not np.isfinite(x):return None
    return x


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def boolean(series):
    value=series.astype(str).str.lower()
    if not value.isin(['true','false']).all():raise ValueError('Missing or ambiguous Boolean classification')
    return value.eq('true')


def describe(values):
    x=np.asarray(values,float)
    if not np.isfinite(x).all():raise ValueError('Nonfinite output in a declared valid sample')
    return dict(n=len(x),mean=float(np.mean(x)) if len(x) else None,
        median=float(np.median(x)) if len(x) else None,
        minimum=float(np.min(x)) if len(x) else None,maximum=float(np.max(x)) if len(x) else None,
        sample_sd=float(np.std(x,ddof=1)) if len(x)>1 else None,
        q05=float(np.quantile(x,.05)) if len(x) else None,q95=float(np.quantile(x,.95)) if len(x) else None)


def label(event):
    date=str(event['start_utc'])[:10]
    return date+(' (benchmark)' if event.get('benchmark') else '')


def load_inputs():
    manifest_path=HERE/'events/manifest.json'
    manifest=json.loads(manifest_path.read_text())
    events=manifest['events']
    if len({e['event_id'] for e in events})!=len(events):raise ValueError('Duplicate event IDs')
    benchmarks=[e for e in events if e.get('benchmark')]
    if len(benchmarks)!=1:raise ValueError('Exactly one unchanged benchmark event is required')
    source_plan_path=SOURCE/'results/study_plan.json'
    plan=json.loads(source_plan_path.read_text())
    new_plan_path=HERE/'study_plan.json'
    new_plan=json.loads(new_plan_path.read_text())
    if new_plan['source_parameter_plan_sha256']!=sha(source_plan_path) or new_plan['events_manifest_sha256']!=sha(manifest_path):
        raise ValueError('Event or source-parameter plan differs from the frozen new study plan')
    expected=pd.DataFrame(plan['mc']).set_index('sample_id')
    if len(expected)!=64 or expected.index.duplicated().any():raise ValueError('Original 64-draw plan changed')
    frames={};sources=[]
    for event in events:
        eid=event['event_id'];path=HERE/'results/ensemble'/eid/'paired_results.csv'
        if sha(Path(event['forcing_file']))!=event['forcing_sha256']:raise ValueError(f'{eid}: saved forcing checksum changed')
        frame=pd.read_csv(path)
        if len(frame)!=64 or frame.sample_id.duplicated().any() or set(frame.sample_id)!=set(expected.index):
            raise ValueError(f'{eid}: preserve all 64 unique original sample IDs')
        if 'event_id' not in frame or not frame.event_id.eq(eid).all():raise ValueError(f'{eid}: event IDs mismatch')
        if not np.all(frame.dt_s==15):raise ValueError('Do not mix primary and refinement evaluations')
        frame=frame.set_index('sample_id').loc[expected.index].reset_index()
        for column in FACTORS:
            if not np.allclose(frame[column],expected[column],rtol=2e-14,atol=0):raise ValueError(f'{eid}: changed parameter draw {column}')
        for row in frame.itertuples(index=False):
            p=expected.loc[row.sample_id,'physical_values']
            mapping=dict(K_u=p['K_u'],k_ud_1=p['k_ud'][0],k_ud_2=p['k_ud'][1],k_ur=p['k_ur'],R_r=p['R_r'])
            for column,value in mapping.items():
                if not np.isclose(getattr(row,column),value,rtol=2e-14,atol=0):raise ValueError(f'{eid}/{row.sample_id}: altered {column}')
        valid=boolean(frame.paired_valid)
        if (valid & ~(boolean(frame.wetland_valid)&boolean(frame.no_wetland_valid))).any():raise ValueError('A valid pair includes an invalid model')
        if not boolean(frame.loc[valid,'paired_weather_inputs_exact']).all():raise ValueError(f'{eid}: a valid pair lacks identical saved weather')
        # One unchanged feedback policy need not produce identical actual gate
        # histories in structurally different models. Differences are outcomes.
        boolean(frame.loc[valid,'paired_gate_history_agreement'])
        for metric in METRICS:
            computed=frame['no_wetland_'+metric]-frame['wetland_'+metric]
            name='delta_'+metric
            if name in frame and not np.allclose(frame.loc[valid,name],computed[valid],rtol=1e-10,atol=1e-10):raise ValueError('Incorrect signed paired difference: '+name)
            frame[name]=computed.where(valid,np.nan)
        for case in CASES:
            okay=boolean(frame[case+'_valid'])
            peaks=frame.loc[okay,case+'_peak_stage_m'];initial=frame.loc[okay,case+'_initial_stage_m']
            if not np.isfinite(peaks).all() or np.any(peaks<initial-1e-12):raise ValueError('Initial state omitted from full-window peak')
            rise=frame.loc[okay,case+'_maximum_rise_above_initial_m']
            if not np.allclose(rise,peaks-initial,rtol=0,atol=1e-10):raise ValueError('Maximum rise definition changed')
            q=np.minimum(3.,np.maximum(peaks-97.15,0)/frame.loc[okay,'R_r'])
            if not np.allclose(frame.loc[okay,case+'_peak_discharge_m3_s'],q,rtol=0,atol=1e-10):raise ValueError('Unchanged outlet rating law mismatch')
            for metric in METRICS:
                if not np.isfinite(frame.loc[okay,case+'_'+metric]).all():raise ValueError('Valid model output missing: '+case+'_'+metric)
        if (~valid).any() and frame.loc[~valid,'failure_reason'].fillna('').str.strip().eq('').any():raise ValueError('Invalid pair requires a specific recorded reason')
        uncapped=valid & frame.wetland_peak_discharge_m3_s.lt(3.) & frame.no_wetland_peak_discharge_m3_s.lt(3.)
        residual=frame.delta_peak_stage_m-frame.R_r*frame.delta_peak_discharge_m3_s
        frame['both_peak_outlets_uncapped']=uncapped
        frame['uncapped_peak_rating_identity_residual_m']=residual.where(uncapped,np.nan)
        if np.any(abs(residual[uncapped])>1e-10):raise ValueError('Uncapped delta_H = R_r delta_Q identity failed')
        frames[eid]=frame
        sources.append(dict(event_id=eid,paired_csv=str(path.relative_to(HERE)),paired_csv_sha256=sha(path)))
    return events,frames,dict(manifest_sha256=sha(manifest_path),study_plan_sha256=sha(new_plan_path),original_parameter_plan_sha256=sha(source_plan_path),analysis_source_sha256=sha(__file__),official_revision=OFFICIAL_REVISION,official_plot_source_hashes=OFFICIAL_HASHES,inputs=sources)


def tables_for_events(events,frames,tables):
    stats=[];accounting=[];failures=[];safety=[];associations=[];diagnostics=[];rating=[];storage=[]
    for event in events:
        eid=event['event_id'];f=frames[eid];valid=boolean(f.paired_valid)
        accounting.append(dict(event_id=eid,case='paired',planned=64,valid=int(valid.sum()),invalid=int((~valid).sum()),
            valid_unsafe=int((valid & (f.wetland_unsafe.fillna(False)|f.no_wetland_unsafe.fillna(False))).sum()),
            gate_history_disagreement_count=int((~boolean(f.loc[valid,'paired_gate_history_agreement'])).sum()),
            maximum_paired_gate_opening_difference_m=float(f.loc[valid,'maximum_paired_gate_opening_difference_m'].max()) if valid.any() else None,
            paired_latch_mismatch_count=int(f.loc[valid,'paired_latch_mismatch_count'].sum()),
            denominator='All 64 planned parameter draws; no resampling; actual gate differences retained'))
        for case in CASES:
            okay=boolean(f[case+'_valid']);unsafe=boolean(f.loc[okay,case+'_unsafe'])
            accounting.append(dict(event_id=eid,case=case,planned=64,valid=int(okay.sum()),invalid=int((~okay).sum()),
                valid_unsafe=int(unsafe.sum()),initial_peak_count=int((okay & f[case+'_peak_stage_time_s'].eq(0)).sum()),
                outlet_cap_peak_count=int((okay & f[case+'_peak_discharge_m3_s'].ge(3.)).sum()),
                minimum_safety_margin_m=float(f.loc[okay,case+'_minimum_safety_margin_m'].min()) if okay.any() else None,
                maximum_relative_mass_residual=float(f.loc[okay,case+'_balance_relative'].max()) if okay.any() else None,
                denominator='All model-valid trajectories; unsafe and initial-peak outcomes retained'))
            for row in f.itertuples(index=False):
                record=row._asdict();path=HERE/'results/ensemble'/eid/f'{row.sample_id}__{case}__15s.json'
                if not path.exists() and event.get('benchmark'):path=SOURCE/'results/monte_carlo'/path.name
                if not path.exists():raise ValueError('Missing model diagnosis: '+str(path))
                raw=json.loads(path.read_text())
                diagnostics.append(dict(event_id=eid,sample_id=row.sample_id,case=case,classification=record[case+'_classification'],
                    valid=record[case+'_valid'],unsafe=record[case+'_unsafe'],full_horizon_completed=raw.get('full_integration_horizon_completed'),
                    failure_reason=raw.get('failure_reason',record.get(case+'_failure_reason')),failure_message=raw.get('failure_message'),
                    failure_stage=raw.get('failure_stage'),first_failure_time_s=raw.get('first_failure_time_s'),
                    solver_termination=raw.get('solver_termination'),relative_mass_residual=raw.get('balance_relative'),
                    relative_mass_tolerance=raw.get('balance_relative_tolerance'),source=str(path)))
                storage_keys=['initial_total_storage_m3','final_total_storage_m3','total_storage_change_m3','external_net_integrated_m3','end_storage_balance_residual_m3']
                storage_values={key:record.get(case+'_'+key) for key in storage_keys}
                if pd.isna(storage_values['initial_total_storage_m3']) and raw.get('initial_storage_m3') is not None:
                    initial_values=np.asarray(raw['initial_storage_m3'],float)
                    if np.isfinite(initial_values).all():storage_values['initial_total_storage_m3']=float(initial_values.sum())
                storage.append(dict(event_id=eid,sample_id=row.sample_id,case=case,valid=record[case+'_valid'],
                    **storage_values,total_outlet_volume_m3=record.get(case+'_total_outlet_volume_m3'),
                    initial_storage_source='Saved model diagnosis or saved trajectory; no antecedent-state reconstruction',
                    final_storage_availability='Saved trajectory-derived accounting' if pd.notna(record.get(case+'_final_total_storage_m3')) else 'Not available in preserved benchmark summary or invalid trajectory; not reconstructed',
                    interpretation='All-compartment storage and external net balance; internal flows are not added again'))
                for variable,limit in LIMITS.items():
                    applicable=case=='wetland' or not variable.startswith('H_w')
                    values=[record.get(f'{case}_{variable}_{where}_peak_m',np.nan) for where in ['node','rhs']]
                    finite=[float(x) for x in values if pd.notna(x)]
                    peak=max(finite) if finite else np.nan
                    if applicable and record[case+'_valid'] and len(finite)!=2:raise ValueError('Valid model lacks output-node/RHS safety audit')
                    safety.append(dict(event_id=eid,sample_id=row.sample_id,case=case,variable=variable,applicable=applicable,
                        threshold_m=limit if applicable else None,peak_m=peak if applicable else None,
                        margin_m=limit-peak if applicable else None,exceeded=bool(peak>limit) if applicable and finite else None,
                        valid=record[case+'_valid']))
        for prefix in [*CASES,'delta']:
            for metric,unit in METRICS.items():
                column=prefix+'_'+metric;values=f.loc[valid,column]
                record=dict(event_id=eid,metric=column,unit=unit,planned=64,valid_pairs=int(valid.sum()),
                    conditional_on_valid_pairs=bool((~valid).any()),**describe(values))
                if prefix=='delta':record.update(positive_count=int(values.gt(0).sum()),zero_count=int(values.eq(0).sum()),negative_count=int(values.lt(0).sum()))
                stats.append(record)
        for factor in FACTORS:
            for outcome in ['delta_peak_stage_m','delta_peak_discharge_m3_s']:
                x=f.loc[valid,factor].to_numpy();y=f.loc[valid,outcome].to_numpy()
                defined=len(x)>=3 and np.ptp(x)>0 and np.ptp(y)>0
                rho=float(spearmanr(x,y).statistic) if defined else None
                associations.append(dict(event_id=eid,factor=factor,outcome=outcome,n=len(x),spearman_rho=rho,
                    unavailable_reason=None if defined else 'Fewer than 3 valid pairs or a constant variable',
                    interpretation='Descriptive marginal rank association; not causal, Sobol, a fitted surrogate or an event-probability estimate'))
        invalid=f.loc[~valid,['event_id','sample_id','failure_reason','wetland_classification','no_wetland_classification']]
        failures.extend(invalid.to_dict('records'))
        r=f.loc[f.both_peak_outlets_uncapped,'uncapped_peak_rating_identity_residual_m']
        rating.append(dict(event_id=eid,valid_pairs=int(valid.sum()),uncapped_pairs=len(r),
            maximum_absolute_identity_residual_m=float(abs(r).max()) if len(r) else None,
            identity='delta_H_peak = R_r * delta_Q_peak only when both peaks are on the uncapped linear rating branch'))
    all_pairs=pd.concat(frames.values(),ignore_index=True)
    all_pairs.to_csv(tables/'all_event_paired_metrics.csv',index=False)
    outputs={'outcome_summaries':stats,'event_accounting':accounting,'safety_checks':safety,
        'trajectory_diagnostics':diagnostics,'parameter_rank_associations':associations,'outlet_rating_identity':rating,'water_storage_balance':storage}
    for name,rows in outputs.items():pd.DataFrame(rows).to_csv(tables/f'{name}.csv',index=False)
    pd.DataFrame(failures,columns=['event_id','sample_id','failure_reason','wetland_classification','no_wetland_classification']).to_csv(tables/'invalid_pairs.csv',index=False)
    pd.DataFrame(events).to_csv(tables/'event_weather_summary.csv',index=False)
    return pd.DataFrame(stats),pd.DataFrame(accounting),pd.DataFrame(associations)


def matched_tables(events,frames,tables):
    benchmark=next(e['event_id'] for e in events if e.get('benchmark'))
    base=frames[benchmark].set_index('sample_id');rows=[];summary=[]
    for event in events:
        eid=event['event_id']
        if eid==benchmark:continue
        current=frames[eid].set_index('sample_id').loc[base.index]
        shared=boolean(base.paired_valid)&boolean(current.paired_valid)
        for sid in base.index:
            row=dict(event_id=eid,benchmark_event_id=benchmark,sample_id=sid,shared_valid=bool(shared.loc[sid]),
                     benchmark_valid=bool(base.loc[sid,'paired_valid']),event_valid=bool(current.loc[sid,'paired_valid']))
            for metric in METRICS:
                name='delta_'+metric
                row['benchmark_'+name]=base.loc[sid,name] if shared.loc[sid] else None
                row['event_'+name]=current.loc[sid,name] if shared.loc[sid] else None
                row['change_'+name]=current.loc[sid,name]-base.loc[sid,name] if shared.loc[sid] else None
            rows.append(row)
        for metric in METRICS:
            column='delta_'+metric
            values=current.loc[shared,column]-base.loc[shared,column]
            summary.append(dict(event_id=eid,benchmark_event_id=benchmark,metric=column,planned=64,shared_valid=int(shared.sum()),
                comparison='Event paired effect minus benchmark paired effect, matched by original parameter sample_id',
                **describe(values),positive_count=int(values.gt(0).sum()),zero_count=int(values.eq(0).sum()),negative_count=int(values.lt(0).sum())))
    pd.DataFrame(rows).to_csv(tables/'matched_event_changes.csv',index=False)
    pd.DataFrame(summary).to_csv(tables/'matched_event_change_summaries.csv',index=False)
    return pd.DataFrame(rows)


def weather_figure(events,plots):
    fields=['rainfall_total_mm','max_hourly_mm','max_24h_mm']
    titles=['Total window rainfall','Maximum hourly rainfall','Maximum rolling 24 h rainfall']
    fig,axes=plt.subplots(1,3,figsize=(13,4.2),sharey=True,layout='constrained')
    labels=[label(e) for e in events];colors=COLORS[:len(events)]
    for ax,field,title in zip(axes,fields,titles):
        values=[e.get(field,np.nan) for e in events]
        ax.barh(labels,values,color=colors)
        for i,value in enumerate(values):
            if pd.notna(value):ax.text(float(value),i,f' {value:.1f}',va='center',fontsize=8)
        ax.set(title=title,xlabel='Rainfall (mm)');ax.grid(axis='x',alpha=.2)
        ax.set_xlim(0,max([float(v) for v in values if pd.notna(v)]+[1])*1.25)
    axes[0].invert_yaxis()
    for ax in axes[1:]:ax.tick_params(axis='y',labelleft=False)
    fig.suptitle('Selected fixed historical storm inputs — no occurrence-probability weights')
    fig.savefig(plots/'01_event_rainfall_summary.png',dpi=180);plt.close(fig)


def clopper_pearson(k,n):
    k=np.asarray(k,dtype=int)
    if n<1 or np.any(k<0) or np.any(k>n):raise ValueError('Invalid binomial counts')
    lower=np.zeros(k.shape,dtype=float);upper=np.ones(k.shape,dtype=float)
    low=k>0;high=k<n
    lower[low]=beta.ppf(.025,k[low],n-k[low]+1)
    upper[high]=beta.ppf(.975,k[high]+1,n-k[high])
    return lower,upper


def ecdf_figure(events,frames,plots,tables):
    official,_=_official_functions();records=[];cp_records=[];warnings_seen=[];checks=[]
    sample_sizes=sorted({int(boolean(f.paired_valid).sum()) for f in frames.values()}-{0})
    for n in sample_sizes:
        k=np.arange(n+1);low,high=clopper_pearson(k,n)
        reference=np.array([[binomtest(int(j),n).proportion_ci(.95,method='exact').low,binomtest(int(j),n).proportion_ci(.95,method='exact').high] for j in k])
        error=float(np.max(abs(np.column_stack([low,high])-reference)))
        if error>1e-10:raise ValueError('Beta-quantile intervals disagree with official SciPy Clopper-Pearson implementation')
        checks.append(dict(n=n,all_possible_success_counts_checked=n+1,maximum_absolute_interval_difference=error,scipy_version=scipy.__version__,passed=True))
    fig,axes=plt.subplots(2,3,figsize=(14,8.7))
    for row,(metric,unit,title) in enumerate([('peak_stage_m','m','River peak stage'),('peak_discharge_m3_s','m³/s','Outlet peak discharge')]):
        for col,prefix in enumerate(['wetland','no_wetland','delta']):
            ax=axes[row,col];column=prefix+'_'+metric;observed=[]
            for color,event in zip(COLORS,events):
                eid=event['event_id'];f=frames[eid];x=f.loc[boolean(f.paired_valid),column].to_numpy(float);n=len(x)
                if not n:continue
                observed.extend(x.tolist())
                grid_key=('delta_' if prefix=='delta' else '')+('stage' if metric=='peak_stage_m' else 'discharge')
                grid=CDF_GRIDS[grid_key];k=np.searchsorted(np.sort(x),grid,side='right')
                cp_lower,cp_upper=clopper_pearson(k,n)
                ax.fill_between(grid,cp_lower,cp_upper,color=color,alpha=.10,step='post')
                for j in range(len(grid)):
                    cp_records.append(dict(event_id=eid,metric=column,grid_name=grid_key,threshold=float(grid[j]),count_leq=int(k[j]),n=n,planned=64,
                        empirical_cdf=float(k[j]/n),pointwise_lower=float(cp_lower[j]),pointwise_upper=float(cp_upper[j]),
                        method='Two-sided Clopper-Pearson 95% at a fixed grid location',conditional_on_valid_pairs=n<64))
                unique,counts=np.unique(x,return_counts=True);exact=np.cumsum(counts)/n
                # Duplicate the minimum at F=0 so an atom at the first value
                # has a visible vertical jump instead of starting above zero.
                plot_x=np.r_[unique[0],unique];plot_F=np.r_[0.,exact]
                ax.step(plot_x,plot_F,where='post',color=color,lw=1.5,label=f'{label(event)}: {n}/64')
                if n>=2:
                    temporary=plt.figure();number=temporary.number
                    with plt.rc_context({'text.usetex':False}), warnings.catch_warnings(record=True) as caught:
                        warnings.simplefilter('always')
                        sx,F,lower,upper=official(x,95.,number,x_label=column,save_plots=False)
                    warnings_seen.extend(str(w.message) for w in caught);plt.close(number)
                    ax.plot(sx,F,':',color=color,lw=.9,alpha=.9)
                    for i in range(n):records.append(dict(event_id=eid,metric=column,rank=i+1,n=n,planned=64,
                        sorted_value=float(sx[i]),official_rank_position=float(F[i]),pointwise_lower=float(lower[i]),pointwise_upper=float(upper[i]),
                        ties_present=len(unique)<n,continuous_exact_coverage_asserted=False,displayed_as_shading=False))
                if prefix=='delta':ax.axvline(0,color='0.5',lw=.6,zorder=0)
            if observed:
                lo=min(observed);hi=max(observed);span=max(hi-lo,.01 if metric=='peak_stage_m' else .02)
                if prefix=='delta':lo=min(lo,0.);hi=max(hi,0.)
                ax.set_xlim(lo-.04*span,hi+.04*span)
            ax.set(title={'wetland':'Wetland','no_wetland':'No wetland','delta':'No wetland − wetland'}[prefix],
                   xlabel=(f'Peak-stage difference ({unit})' if metric=='peak_stage_m' else f'Peak-discharge difference ({unit})') if prefix=='delta' else f'{title} ({unit})',ylabel='Cumulative fraction',ylim=(-.02,1.02));ax.grid(alpha=.18)
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.065),ncol=3,fontsize=8)
    fig.suptitle('Separate parameter distributions conditional on each fixed storm',y=.985)
    fig.text(.5,.013,'Solid: ties-aware count/N ECDF. Dotted: official rank/(N+1). Shading: 95% pointwise Clopper–Pearson intervals on fixed grids.\nCoverage is pointwise at fixed grid locations (connecting shade is a visual guide), not simultaneous. Official Beta limits are archived only.',ha='center',fontsize=8)
    fig.subplots_adjust(top=.92,bottom=.21,left=.065,right=.99,hspace=.34,wspace=.25)
    fig.savefig(plots/'02_per_event_peak_distributions.png',dpi=180);plt.close(fig)
    pd.DataFrame(records).to_csv(tables/'official_pointwise_ecdf_arrays.csv',index=False)
    pd.DataFrame(cp_records).to_csv(tables/'fixed_grid_pointwise_cdf_intervals.csv.gz',index=False,compression='gzip')
    pd.DataFrame(checks).to_csv(tables/'clopper_pearson_implementation_checks.csv',index=False)
    return warnings_seen


def matched_figure(events,matched,plots):
    others=[e for e in events if not e.get('benchmark')]
    fig,axes=plt.subplots(2,len(others),figsize=(3.45*len(others),7),squeeze=False,layout='constrained')
    for col,event in enumerate(others):
        part=matched.loc[matched.event_id.eq(event['event_id']) & matched.shared_valid.eq(True)]
        for row,(metric,unit) in enumerate([('delta_peak_stage_m','m'),('delta_peak_discharge_m3_s','m³/s')]):
            ax=axes[row,col];x=part['benchmark_'+metric].to_numpy(float);y=part['event_'+metric].to_numpy(float)
            ax.scatter(x,y,s=18,color=COLORS[col+1],alpha=.75)
            if len(x):
                lo=min(0,float(x.min()),float(y.min()));hi=max(float(x.max()),float(y.max()),1e-8)
                ax.plot([lo,hi],[lo,hi],color='0.5',ls='--',lw=.8)
            ax.axhline(0,color='0.7',lw=.5);ax.axvline(0,color='0.7',lw=.5);ax.grid(alpha=.18)
            ax.set(title=f'{label(event)}\n{len(part)}/64 shared valid draws',xlabel=f'Benchmark paired reduction ({unit})',ylabel=f'Event paired reduction ({unit})')
    fig.suptitle('Same parameter vectors matched across events — no independent-event sampling claim')
    fig.savefig(plots/'03_matched_event_peak_effects.png',dpi=180);plt.close(fig)


def association_figure(events,assoc,plots):
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),sharey=True,layout='constrained')
    for ax,outcome,title in zip(axes,['delta_peak_stage_m','delta_peak_discharge_m3_s'],['Peak-stage reduction','Peak-discharge reduction']):
        for color,event in zip(COLORS,events):
            part=assoc.loc[assoc.event_id.eq(event['event_id']) & assoc.outcome.eq(outcome)].set_index('factor').reindex(FACTORS)
            ax.plot(range(3),part.spearman_rho,'o-',color=color,label=label(event),ms=5,lw=1)
        ax.set(title=title,xticks=range(3),xticklabels=['Infiltration','Routing','Resistance'],ylim=(-1.05,1.05),xlabel='Parameter multiplier')
        ax.axhline(0,color='0.5',lw=.8);ax.grid(alpha=.18)
    axes[0].set_ylabel('Descriptive Spearman rank association');axes[1].legend(fontsize=8)
    fig.suptitle('Association, not causal sensitivity: ΔH = Rᵣ ΔQ on the uncapped outlet branch')
    fig.savefig(plots/'04_parameter_rank_associations.png',dpi=180);plt.close(fig)


def hydrograph_figure(events,frames,plots,tables):
    """Show the first original draw, fixed before inspecting new event outcomes."""
    sample_id='mc_00000'
    fig,axes=plt.subplots(len(events),2,figsize=(13,2.2*len(events)),sharex=True,squeeze=False)
    records=[];stage_values=[];rain_values=[]
    for i,event in enumerate(events):
        eid=event['event_id'];frame=frames[eid].set_index('sample_id');record=frame.loc[sample_id]
        stage_ax,rain_ax=axes[i]
        for case,color,style in [('wetland','#0072B2','-'),('no_wetland','#D55E00','--')]:
            if event.get('benchmark'):
                path=SOURCE/'results/pilot'/f'pilot_mc_00000__{case}__15s.npz'
            else:path=HERE/'results/ensemble'/eid/f'{sample_id}__{case}__15s.npz'
            available=bool(record[case+'_valid'])
            info=dict(event_id=eid,sample_id=sample_id,case=case,valid=available,
                selection='First original parameter draw; declared before inspecting historical-event hydrographs',
                trajectory_source=str(path),forcing_sha256=event['forcing_sha256'])
            if available:
                if not path.exists():raise ValueError('Valid illustrative trajectory missing: '+str(path))
                with np.load(path) as z:
                    time=z['time'].copy();stage=z['downstream_stage'].copy()
                if time.ndim!=1 or stage.shape!=time.shape or not np.isfinite(stage).all():raise ValueError('Invalid saved hydrograph shape/values')
                if time[0]!=0 or time[-1]!=144*3600 or np.any(np.diff(time)<=0):raise ValueError('Illustrative hydrograph does not cover the original full 144 h')
                if not np.isclose(stage.max(),record[case+'_peak_stage_m'],rtol=0,atol=1e-10):raise ValueError('Illustrative trajectory disagrees with its saved MC peak')
                stage_ax.plot(time/3600,stage,color=color,ls=style,lw=1.2,label=case.replace('_',' ').title())
                stage_values.extend([float(stage.min()),float(stage.max())])
                info.update(trajectory_sha256=sha(path),time_start_s=float(time[0]),time_end_s=float(time[-1]),
                    initial_stage_m=float(stage[0]),full_window_peak_stage_m=float(stage.max()),peak_time_h=float(time[np.argmax(stage)]/3600))
            else:
                stage_ax.text(.98,.93 if case=='wetland' else .82,case.replace('_',' ')+': no valid complete trajectory',
                    ha='right',va='top',transform=stage_ax.transAxes,fontsize=8,color=color)
                info['failure_reason']=record.get(case+'_failure_reason')
            records.append(info)
        initial=float(record.wetland_initial_stage_m)
        stage_ax.axhline(initial,color='0.45',ls=':',lw=.9,label='Initial stage (not flood threshold)')
        stage_ax.set_ylabel('River stage (m)');stage_ax.grid(alpha=.2)
        stage_ax.set_title(f"{label(event)} — window starts {str(event['start_utc'])[11:16]} UTC",fontsize=10)
        forcing=pd.read_csv(Path(event['forcing_file']))
        if event.get('benchmark'):
            t=forcing.time_s.to_numpy(float)/3600;rain=forcing.P.to_numpy(float)*3.6e6
            rain_ax.plot(t,rain,color='#555555',lw=1);rain_ax.fill_between(t,0,rain,color='#555555',alpha=.2)
        else:
            ts=pd.to_datetime(forcing.timestamp_utc,utc=True);start=pd.Timestamp(event['start_utc'])
            t=(ts-start).dt.total_seconds().to_numpy()/3600;rain=forcing.precipitation_mm_per_h.to_numpy(float)
            rain_ax.stairs(rain,np.r_[t,t[-1]+1],color='#555555',fill=True,alpha=.45)
        rain_values.extend(rain.tolist());rain_ax.set_ylabel('Rainfall (mm/h)');rain_ax.grid(alpha=.2)
        rain_ax.set_title('Preserved event rainfall; no peak-time alignment',fontsize=10)
    if stage_values:
        lo=min(stage_values);hi=max(stage_values);pad=max(.03,(hi-lo)*.06)
        for ax in axes[:,0]:ax.set_ylim(lo-pad,hi+pad)
    for ax in axes[:,1]:ax.set_ylim(0,max(rain_values+[1.])*1.08)
    for ax in axes.flat:ax.set_xlim(0,144);ax.set_xticks(np.arange(0,145,24))
    for ax in axes[-1]:ax.set_xlabel('Hours from the saved event-window start')
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.035),ncol=3,fontsize=9)
    fig.suptitle('One illustrative parameter draw: mc_00000 — selected in advance, not representative of all 64',y=.995,fontsize=13)
    fig.text(.5,.012,'Full 144 h shown including initial conditions. The benchmark retains its original 120 s rain ramps; other events retain hourly constant rain.\nAll panels use the same elapsed-window axes; storms and hydrographs are not shifted to their peaks.',ha='center',fontsize=8)
    fig.subplots_adjust(left=.075,right=.99,top=.95,bottom=.105,hspace=.4,wspace=.22)
    fig.savefig(plots/'05_fixed_draw_hydrographs.png',dpi=180);plt.close(fig)
    pd.DataFrame(records).to_csv(tables/'illustrative_hydrograph_sources.csv',index=False)


def main():
    events,frames,provenance=load_inputs()
    tables,plots=HERE/'tables',HERE/'plots';tables.mkdir(exist_ok=True);plots.mkdir(exist_ok=True)
    stats,accounting,assoc=tables_for_events(events,frames,tables)
    matched=matched_tables(events,frames,tables)
    weather_figure(events,plots)
    plot_warnings=ecdf_figure(events,frames,plots,tables)
    matched_figure(events,matched,plots);association_figure(events,assoc,plots)
    hydrograph_figure(events,frames,plots,tables)
    pilots=[]
    for event in events:
        path=HERE/'results/pilot'/event['event_id']/'numerical_comparison.csv'
        if path.exists():
            f=pd.read_csv(path);f['event_id']=event['event_id'];pilots.append(f)
    if pilots:pd.concat(pilots,ignore_index=True).to_csv(tables/'numerical_pilot_comparisons.csv',index=False)
    pilot_summary_path=HERE/'results/pilot/pilot_summary.json'
    pilot_summary=json.loads(pilot_summary_path.read_text()) if pilot_summary_path.exists() else None
    audit=dict(**provenance,event_count=len(events),planned_parameter_draws_per_event=64,
        distinct_parameter_vectors=64,event_occurrence_probabilities_assigned=False,
        illustrative_hydrograph_sample_id='mc_00000',illustrative_sample_selected_before_new_hydrograph_inspection=True,parameter_law='Original independent exploratory log-uniform multipliers on [0.5,2]',
        common_initial_conditions_reconstruct_antecedent_weather=False,weather_or_parameter_resampling_performed=False,
        physical_model_runs_performed_by_analysis=0,surrogate_fit_performed=False,unsafe_or_negative_effect_filter_applied=False,
        conditional_valid_subsets_used={e['event_id']:int(boolean(frames[e['event_id']].paired_valid).sum())<64 for e in events},
        numerical_pilot_event_ids=[str(f.event_id.iloc[0]) for f in pilots],numerical_pilot_summary=pilot_summary,
        numerical_pilot_summary_sha256=sha(pilot_summary_path) if pilot_summary_path.exists() else None,plotting_warnings=plot_warnings,
        cross_event_caveat='Benchmark preserves original B trace-zero lower bounds and 120 s ramps; new complete C-derived windows use their declared hourly adapter. Processing differs and is not normalized after seeing results.',
        ecdf_scope='Separate fixed-event distributions; exact count/N curves and official rank/(N+1) dotted lines. Shading uses atom-safe 95% Clopper-Pearson CDF intervals at fixed grid locations; pointwise, not simultaneous. Official order-statistic Beta intervals retained only as diagnostics.',
        cdf_interval_method=dict(name='Clopper-Pearson two-sided exact binomial intervals',confidence=.95,scipy_version=scipy.__version__,documentation=CP_DOC_URL,
            interpretation='Conservative pointwise coverage conditional on exploratory iid parameter law; supports atoms; no weather-occurrence probability',
            grids={name:dict(start=float(grid[0]),end=float(grid[-1]),step=float(grid[1]-grid[0]),count=len(grid)) for name,grid in CDF_GRIDS.items()},
            fixed_grids_selected_before_ensemble_completion=True,plot_zoom_changes_computation_grid=False),
        output_files=[str(p.relative_to(HERE)) for p in sorted(tables.glob('*'))]+[str(p.relative_to(HERE)) for p in sorted(plots.glob('*.png'))])
    (tables/'analysis_audit.json').write_text(json.dumps(clean(audit),indent=2,allow_nan=False)+'\n')
    captions='''# Historical-storm robustness figures

1. **Event rainfall summary.** Each bar describes one selected fixed historical storm window. No frequency or probability weights are assigned. The preserved 2022 benchmark retains its original processing, including trace-zero lower bounds and 120 s ramps; new complete C-derived windows use their separately documented hourly adapter. The selected additional windows have weaker recorded 24-hour rainfall than the benchmark, so they do not establish robustness to more severe storms.
2. **Separate peak distributions.** The same original 64 exploratory parameter vectors are propagated for each event. Solid curves are ties-aware count/N empirical distributions and dotted curves are unchanged official rank/(N+1) plotting positions. Because exact zero effects and initial-peak atoms occur, shading instead uses two-sided 95% Clopper–Pearson intervals for the CDF at predeclared fixed grid locations. These binomial intervals support atoms and have conservative pointwise coverage at fixed grid locations; connecting shading is a visual guide, not a simultaneous or between-grid confidence claim. This is an explicit statistical adaptation, not the official utility's order-statistic interval. The unchanged official Beta output is preserved only in a diagnostic CSV, without a continuous-coverage claim for tied outcomes. See [CDF interval method](cdf_interval_method.md) for grids, formulas and the official SciPy reference. Invalid-pair counts remain in the tables; any reduced denominator is conditional on complete valid pairs. Unsafe valid outcomes, negative effects and initial-state peaks are retained. No downstream flood threshold is defined.
3. **Matched event effects.** Every point matches the same parameter sample_id across the benchmark and another event. Positive reductions mean no-wetland minus wetland. The dashed identity line denotes equal paired effect. These are controlled forcing comparisons with fixed initial states, not independent sampled events or reconstructions of antecedent conditions. No pooled weather probability is calculated.
4. **Parameter rank associations.** Spearman coefficients describe marginal associations within each event's retained parameter draws; they are not causal effects, Sobol indices, a new fitted model or a new one-at-a-time experiment. On the uncapped outlet branch, ΔH = Rᵣ ΔQ, so stage associations with resistance are partly inherent in the rating conversion. Constant-output associations are unavailable, not forced to zero.
5. **Predeclared illustrative hydrographs.** The first original parameter draw, `mc_00000`, was selected before inspecting the new historical-event hydrographs. Both model trajectories and the initial-stage reference are shown for each complete 144-hour window, alongside its preserved rainfall input. Window elapsed times are not shifted to rainfall or river peaks. This one draw illustrates dynamics and initial-peak floors; it is not claimed representative of all 64 draws. The benchmark trajectory is read from its preserved original numerical-pilot NPZ, checked against the benchmark MC summary, and is not rerun. Missing invalid trajectories are reported, not replaced.

Peak metrics cover the complete unchanged event window including t=0. Initial-level duration and initial-discharge excess volume are response benchmarks, not flood-exceedance measures. The wetland and no-wetland configurations retain their different accepted initial water volumes; their structural contrast is conditional on that choice. The same feedback policy is retained, but actual gate histories may differ with the modeled states. Such differences are reported without invalidating or discarding trajectories; the structural contrast is not automatically a pure wetland-storage effect under identical releases. The 64 reused parameter vectors create cross-event dependence, which is preserved through matched comparisons. No new weather, parameter draw, physical run, optimization or surrogate fit is created by this analysis script.
'''
    (HERE/'figure_captions.md').write_text(captions)
    print(json.dumps(clean(accounting.to_dict('records')),indent=2,allow_nan=False))


if __name__=='__main__':main()
