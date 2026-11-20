"""Report saved Stage D execution artifacts; never invoke a model or solver.

Require counts.json for nominal, low and high before assembling final outputs.
The report preserves cold exclusions and invalid physical peaks as missing.
"""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import sys
from collections import Counter
sys.dont_write_bytecode=True
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/cee690-d-revision-report-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
RESULTS=HERE/'results'
TABLES=HERE/'tables'
PLOTS=HERE/'plots'
GROUPS=('nominal','low','high')


def read(path):return json.loads(Path(path).read_text())
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,data):Path(path).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
def csv(name,rows):pd.DataFrame(rows).to_csv(TABLES/name,index=False)
def valid(d):return bool(d.get('full_valid_trajectory_available',False))
def safe(d):return valid(d) and not (d.get('safety_violation') or d.get('safety_violation_any_successful_rhs'))
def deterministic(d):return d['scenario_id']=='deterministic_september_2022'
def fault(d):return d.get('first_violation') or {}
def ftime(d):return d.get('first_failure_time_s',fault(d).get('time_s'))
def case_label(d):return d['_case']+(' '+d['_design'] if d['_case']=='wetland' else '')


def load_completed():
    loaded=[];plans={}
    for group in GROUPS:
        directory=RESULTS/group
        if not (directory/'counts.json').is_file():
            raise RuntimeError('Runs not complete; missing counts.json: '+group)
        plan=read(directory/'execution_plan.json');plans[group]=plan
        local=[]
        for item in plan['plans']:
            sid=item['scenario_id'] or 'deterministic_september_2022'
            name=f'{sid}__{item["case"]}__{item["design"]}__{item["dt_s"]}s'
            path=directory/(name+'.json')
            if not path.is_file():raise RuntimeError('Missing completed case record: '+str(path))
            d=read(path)
            d.update(_group=group,_case=item['case'],_design=item['design'],_path=str(path),_npz=str(path.with_suffix('.npz')))
            if d.get('scenario_id')!=sid or d.get('dt_s')!=item['dt_s']:
                raise ValueError('Case identity mismatch: '+str(path))
            if valid(d) and not path.with_suffix('.npz').exists():raise ValueError('Missing complete trajectory NPZ: '+str(path))
            if not valid(d) and any(d.get(k) is not None for k in ('physical_river_peak_m','physical_river_peak_discharge_m3_s')):
                raise ValueError('Invalid trajectory was assigned a physical peak: '+str(path))
            local.append(d)
        counts=Counter(f'{d["_case"]} / {d["classification"]}' for d in local)
        if dict(counts)!=read(directory/'counts.json'):
            raise ValueError('Case counts do not agree with completion record: '+group)
        loaded.extend(local)
    return loaded,plans


def original_for(d):
    if deterministic(d):
        ref=d.get('saved_B_baseline_comparison')
        if not ref:raise ValueError('Missing revised-versus-original comparison')
        return ref
    return read(ROOT/f'D/g5_diagnostics/results/original_20_120s/{d["scenario_id"]}__{d["_case"]}__best.json')


def safety_rows(d):
    rows=[]
    for item in (d.get('all_rhs_observations') or {}).get('safety_successful_rhs',[]):
        rows.append(dict(group=d['_group'],dt_s=d['dt_s'],scenario_id=d['scenario_id'],case=d['_case'],design=d['_design'],
            classification=d['classification'],full_trajectory_valid=valid(d),full_horizon_safety_feasible=safe(d) if valid(d) else None,
            peak_scope='Every successful official RHS/RK evaluation over complete event' if valid(d) else 'Successful RHS/RK evaluations before interruption only; full-event safety unknown',
            **item))
    return rows


def deterministic_rows(data):
    rows=[]
    for d in data:
        if not deterministic(d):continue
        baseline=original_for(d)
        row=dict(group=d['_group'],Cw_US=d['hydraulic_variant_settings']['Cw_US'],dt_s=d['dt_s'],case=d['_case'],design=d['_design'],
            classification=d['classification'],valid=valid(d),safety_feasible=safe(d) if valid(d) else None,
            original_B_peak_Hr_m=baseline['saved_physical_river_peak_m'],revised_D_peak_Hr_m=d.get('physical_river_peak_m'),
            revised_minus_original_peak_m=baseline.get('river_peak_change_m'),
            original_B_peak_Qout_m3_s=baseline['saved_physical_peak_outlet_m3_s'],revised_D_peak_Qout_m3_s=d.get('physical_river_peak_discharge_m3_s'),
            balance_relative=d.get('balance_relative'),balance_max_m3=d.get('balance_max_m3'),
            derivative_mass_identity_max_m3_s=d.get('all_rhs_observations',{}).get('derivative_physical_mass_identity_max_abs_m3_s'),
            identical_initial_states=baseline.get('identical_initial_states'),identical_inputs_at_original_nodes=baseline.get('identical_inputs_at_original_nodes'),
            comparison_kind=baseline['kind'],original_saved_file=baseline['saved_trajectory'],original_saved_sha256=baseline['saved_sha256'])
        for s in safety_rows(d):row[s['variable']+'_peak_m']=s['peak_m'];row[s['variable']+'_margin_m']=s['margin_m']
        rows.append(row)
    return rows


def scenario_rows(data):
    rows=[]
    for d in data:
        if deterministic(d):continue
        old=original_for(d);v=fault(d);o=fault(old);weather=Path(d['forcing_file'])
        if sha(weather)!=d['forcing_file_sha256'] or d['forcing_file_sha256']!=old['forcing_file_sha256']:
            raise ValueError('Scenario weather is not byte-identical to original: '+d['scenario_id'])
        if d['original_seed']!=old['original_seed']:raise ValueError('Seed provenance mismatch')
        changed=None if ftime(d) is None else ftime(d)-o['time_s']
        rows.append(dict(group=d['_group'],Cw_US=d['hydraulic_variant_settings']['Cw_US'],dt_s=d['dt_s'],scenario_id=d['scenario_id'],case=d['_case'],
            seed=json.dumps(d['original_seed']),weather_sha256=d['forcing_file_sha256'],original_weather_hash_verified=True,
            classification=d['classification'],failure_reason=d.get('failure_reason'),failure_stage=d.get('failure_stage'),
            ode_started=d.get('ode_started',False),valid=valid(d),full_horizon_safety_feasible=safe(d) if valid(d) else None,
            first_failure_time_s=ftime(d),first_failure_utc=d.get('first_failure_timestamp_utc',v.get('timestamp_utc')),
            branch=v.get('branch'),H_m=v.get('water_depth_m'),a_m=v.get('opening_height_m'),H_over_a=v.get('H_over_a'),
            upstream_m=v.get('upstream_stage_m'),tailwater_m=v.get('tailwater_stage_m'),head_length_ratio=v.get('head_length_ratio'),
            revised_applicability_residual=v.get('applicability_residual'),revised_residual_units=v.get('applicability_residual_units'),
            legacy_G5_residual_m=v.get('legacy_G5_residual_m'),legacy_G5_role='Descriptive only for revised law',
            old_120s_classification=old['classification'],old_120s_first_failure_time_s=o['time_s'],
            detection_time_shift_from_old_120s_s=changed,
            shift_interpretation='Later interruption is not a complete or validated trajectory',
            physical_peak_Hr_m=d.get('physical_river_peak_m'),physical_peak_Qout_m3_s=d.get('physical_river_peak_discharge_m3_s'),
            nonpositive_temperature_hours=d.get('nonpositive_temperature_hours'),
            first_failure_detection=v.get('time_s_meaning',v.get('detection'))))
    return rows


def contrast_rows(data):
    rows=[]
    for group in GROUPS:
        for dt in sorted({d['dt_s'] for d in data if d['_group']==group}):
            selected=[d for d in data if d['_group']==group and d['dt_s']==dt and deterministic(d)]
            no=next(d for d in selected if d['_case']=='no_wetland')
            for wet in (d for d in selected if d['_case']=='wetland'):
                both=valid(wet) and valid(no);feasible=safe(wet) and safe(no)
                rows.append(dict(group=group,dt_s=dt,wetland_design=wet['_design'],both_trajectories_valid=both,
                    both_safety_feasible=feasible if both else None,no_wetland_peak_m=no.get('physical_river_peak_m'),wetland_peak_m=wet.get('physical_river_peak_m'),
                    no_wetland_minus_wetland_peak_m=no['physical_river_peak_m']-wet['physical_river_peak_m'] if both else None,
                    interpretation=('Fixed-control structural comparison; unsafe D trajectories, not an admissible design benefit' if both and not feasible else
                                    'Fixed-control structural comparison; no new D optimization was performed'),
                    saved_B_design_labels_are_not_D_optima=True,initial_extra_wetland_water_m3=5500))
    return rows


def resolution_rows(data):
    rows=[];nom=[d for d in data if d['_group']=='nominal']
    keys=sorted({(d['scenario_id'],d['_case'],d['_design']) for d in nom})
    for sid,case,design in keys:
        values={d['dt_s']:d for d in nom if (d['scenario_id'],d['_case'],d['_design'])==(sid,case,design)}
        for coarse,fine in ((120,60),(60,30)):
            a,b=values[coarse],values[fine]
            row=dict(scenario_id=sid,case=case,design=design,coarse_dt_s=coarse,fine_dt_s=fine,
                coarse_classification=a['classification'],fine_classification=b['classification'],same_classification=a['classification']==b['classification'],
                both_complete=valid(a) and valid(b),coarse_first_failure_s=ftime(a),fine_first_failure_s=ftime(b),
                first_failure_time_difference_s=None if ftime(a) is None or ftime(b) is None else ftime(b)-ftime(a))
            if valid(a) and valid(b):
                with np.load(a['_npz'],allow_pickle=False) as ca,np.load(b['_npz'],allow_pickle=False) as fi:
                    index=np.searchsorted(fi['time'],ca['time'])
                    if np.any(index>=len(fi['time'])) or not np.array_equal(fi['time'][index],ca['time']):raise ValueError('Unmatched resolution grids')
                    names=ca['output_names'].tolist();stages=[i for i,n in enumerate(names) if n.startswith('H_')]
                    row.update(max_common_grid_river_difference_m=float(np.max(abs(ca['downstream_stage']-fi['downstream_stage'][index]))),
                        max_common_grid_all_stage_difference_m=float(np.max(abs(ca['y'][stages]-fi['y'][stages][:,index]))),
                        max_common_grid_storage_difference_m3=float(np.max(abs(ca['x']-fi['x'][:,index]))),
                        native_grid_river_peak_difference_m=b['physical_river_peak_m']-a['physical_river_peak_m'],
                        note='Finite-grid consistency on complete trajectories; not physical validation or a general error bound')
            else:row['note']='No complete trajectory error comparison; diagnostic first-failure timing only'
            rows.append(row)
    return rows


def coefficient_rows(data,plans):
    rows=[]
    for group in ('low','nominal','high'):
        vals=[d for d in data if d['_group']==group and d['dt_s']==120]
        warm=[d for d in vals if not deterministic(d) and d['classification']!='meteorological_liquid_water_applicability']
        cold=[d for d in vals if not deterministic(d) and d['classification']=='meteorological_liquid_water_applicability']
        det=[d for d in vals if deterministic(d)]
        cw=plans[group]['variant']['Cw_US']
        row=dict(group=group,Cw_US=cw,Cw_SI=cw*np.sqrt(.3048),warm_weather_IDs=len({d['scenario_id'] for d in warm}),cold_weather_IDs=len({d['scenario_id'] for d in cold}),
            warm_case_evaluations=len(warm),warm_case_completions=sum(valid(d) for d in warm),cold_case_preflight_exclusions=len(cold),
            deterministic_case_completions=sum(valid(d) for d in det),deterministic_safe_cases=sum(safe(d) for d in det),
            selection_status='Predetermined sensitivity member, not calibrated or selected for feasibility')
        for d in det:
            key=d['_case']+'_'+d['_design'];row[key+'_Hr_peak_m']=d.get('physical_river_peak_m')
            for s in safety_rows(d):row[key+'_'+s['variable']+'_margin_m']=s['margin_m']
        rows.append(row)
    return rows


def source_rows(data):
    groups={}
    for d in data:
        manifest=d['fixed_configuration']
        merged=dict(manifest['source_sha256']);merged.update(manifest.get('stage_d_source_sha256',{}))
        for rel,digest in merged.items():
            groups.setdefault((rel,digest),set()).add((d['_group'],d['scenario_id'],d['_case'],d['_design'],d['dt_s']))
    rows=[]
    for (rel,digest),uses in sorted(groups.items()):
        path=ROOT/rel;current=sha(path) if path.exists() else None
        rows.append(dict(source=rel,recorded_sha256=digest,current_sha256=current,hash_matches_current=digest==current,recorded_case_uses=len(uses)))
    return rows


def plot_deterministic(data):
    d=next(d for d in data if d['_group']=='nominal' and deterministic(d) and d['_case']=='wetland' and d['_design']=='best' and d['dt_s']==30)
    no=next(d for d in data if d['_group']=='nominal' and deterministic(d) and d['_case']=='no_wetland' and d['dt_s']==30)
    with np.load(d['_npz'],allow_pickle=False) as f:new={k:f[k] for k in f.files}
    with np.load(ROOT/d['saved_B_baseline_comparison']['saved_trajectory'],allow_pickle=False) as f:old={k:f[k] for k in f.files}
    with np.load(no['_npz'],allow_pickle=False) as f:nnew={k:f[k] for k in f.files}
    with np.load(ROOT/no['saved_B_baseline_comparison']['saved_trajectory'],allow_pickle=False) as f:nold={k:f[k] for k in f.files}
    names=new['output_names'].tolist();nnames=nnew['output_names'].tolist()
    fig,axes=plt.subplots(3,2,figsize=(11.5,9),sharex=True,constrained_layout=True)
    limits={'H_d1':101.2,'H_d2':101.0,'H_w1':98.9,'H_w2':98.8}
    for ax,name in zip(axes.flat,['H_d1','H_d2','H_w1','H_w2','H_r','Q_out']):
        j=names.index(name)
        ax.plot(old['time']/3600,old['y'][j],color='#286d8d',lw=1.5,label='Original B, historical B-selected design')
        ax.plot(new['time']/3600,new['y'][j],color='#c2592b',lw=1.5,label='Revised D, historical B-selected design')
        if name in limits:ax.axhline(limits[name],color='#a22d3b',ls='--',lw=1,label='Unchanged safety limit')
        else:
            jj=nnames.index(name)
            ax.plot(nold['time']/3600,nold['y'][jj],color='#585858',ls=':',lw=1.2,label='Original B, no wetland')
            ax.plot(nnew['time']/3600,nnew['y'][jj],color='#8b6caf',ls='-.',lw=1.2,label='Revised D, no wetland')
        ax.set(title=name.replace('_',' '),ylabel='Discharge (m³/s)' if name=='Q_out' else 'Elevation (m)')
        ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=7,loc='lower right');axes[2,0].legend(fontsize=7,loc='upper left')
    for ax in axes[-1]:ax.set_xlabel('Hours since 2022-09-09 00:00 UTC')
    fig.suptitle('Same event and controls; revised nominal hydraulics changes routing\nDam 1 safety is exceeded: lower river peaks are not an admissible optimization gain',fontsize=12)
    fig.savefig(PLOTS/'deterministic_comparison.png',dpi=180);plt.close(fig)


def plot_failures(rows):
    chosen=[r for r in rows if r['group']=='nominal' and r['dt_s']==30 and r['case']=='wetland']
    warm=sorted([r for r in chosen if r['ode_started']],key=lambda r:r['scenario_id'])
    cold=sorted([r for r in chosen if not r['ode_started']],key=lambda r:r['scenario_id'])
    fig,(ax,bx)=plt.subplots(1,2,figsize=(12,7.5),gridspec_kw={'width_ratios':[3,2]},constrained_layout=True)
    y=np.arange(len(warm))
    old=np.array([r['old_120s_first_failure_time_s']/3600 for r in warm]);new=np.array([r['first_failure_time_s']/3600 for r in warm])
    for k in y:ax.plot([old[k],new[k]],[k,k],color='#bdbdbd',lw=1)
    ax.scatter(old,y,s=25,color='#286d8d',label='Original B: unsupported H/a (120 s)')
    ax.scatter(new,y,s=35,marker='x',color='#c2592b',label='Revised D: crest H/L domain (30 s)')
    ax.set_yticks(y,[r['scenario_id'] for r in warm],fontsize=8);ax.invert_yaxis()
    ax.set(xlabel='First detected failure time (h)',title='All 16 warm scenarios still stop')
    ax.grid(axis='x',alpha=.2);ax.set_ylim(len(warm)+2.5,-.5);ax.legend(fontsize=8,loc='lower right')
    cy=np.arange(len(cold));hours=[r['first_failure_time_s']/3600 for r in cold]
    bx.barh(cy,hours,color='#858585',height=.5)
    bx.set_yticks(cy,[r['scenario_id'] for r in cold],fontsize=8);bx.invert_yaxis()
    for k,r in enumerate(cold):bx.text(hours[k]+1,k,f'{r["nonpositive_temperature_hours"]} cold hours',va='center',fontsize=8)
    bx.set_xlim(0,max(hours)*1.5);bx.set(xlabel='First excluded weather interval (h)',title='Four cold inputs: no ODE was started')
    bx.grid(axis='x',alpha=.2)
    fig.suptitle('Original weather IDs retained; later interruption is not successful propagation\nBoth physical cases have the same first-failure times in this nominal comparison',fontsize=12)
    fig.savefig(PLOTS/'scenario_first_failures.png',dpi=180);plt.close(fig)


def plot_resolution_sensitivity(data,coeff):
    fig,axes=plt.subplots(1,3,figsize=(13,4.4),constrained_layout=True)
    labels=[('wetland','reference','A/A reference','#286d8d'),('wetland','best','Historical B-selected design','#c2592b'),('no_wetland','best','No wetland','#675887')]
    for case,design,label,color in labels:
        vals=sorted([d for d in data if d['_group']=='nominal' and deterministic(d) and d['_case']==case and d['_design']==design],key=lambda d:d['dt_s'])
        fine=vals[0]['physical_river_peak_m']
        axes[0].plot([d['dt_s'] for d in vals],[1000*(d['physical_river_peak_m']-fine) for d in vals],'-o',label=label,color=color,ms=4)
        ordered=sorted(coeff,key=lambda r:r['Cw_US']);key=case+'_'+design+'_Hr_peak_m'
        axes[2].plot([r['Cw_US'] for r in ordered],[r[key] for r in ordered],'-o',label=label,color=color,ms=4)
    ordered=sorted(coeff,key=lambda r:r['Cw_US'])
    axes[1].plot([r['Cw_US'] for r in ordered],[100*r['wetland_best_H_d1_margin_m'] for r in ordered],'-o',color='#c2592b')
    axes[1].axhline(0,color='#a22d3b',ls='--',lw=1)
    axes[1].set(xlabel='Cw (US units: ft½/s)',ylabel='Dam 1 safety margin (cm)',title='Fixed-limit safety changes with assumed Cw')
    axes[0].set(xlabel='Timestep (s)',ylabel='Peak difference from 30 s (mm)',title='Nominal finite-grid consistency');axes[0].set_xticks([30,60,120])
    axes[2].set(xlabel='Cw (US units: ft½/s)',ylabel='Physical river peak (m)',title='Predetermined coefficient sensitivity')
    axes[2].legend(fontsize=8)
    for ax in axes:ax.grid(alpha=.2);ax.ticklabel_format(axis='y',useOffset=False)
    fig.suptitle('No coefficient was selected to obtain feasibility; all 16 warm weather inputs still fail',fontsize=12)
    fig.savefig(PLOTS/'resolution_and_coefficient_sensitivity.png',dpi=180);plt.close(fig)


def markdown_table(rows,columns):
    def fmt(v):
        if v is None:return 'Not available'
        if isinstance(v,(float,np.floating)):return f'{v:.8g}'
        return str(v)
    head='| '+' | '.join(label for key,label in columns)+' |\n| '+' | '.join('---' for _ in columns)+' |\n'
    return head+'\n'.join('| '+' | '.join(fmt(row.get(key)) for key,label in columns)+' |' for row in rows)


def main():
    data,plans=load_completed()
    TABLES.mkdir(exist_ok=True);PLOTS.mkdir(exist_ok=True)
    det=deterministic_rows(data);scenarios=scenario_rows(data);contrasts=contrast_rows(data)
    safety=[r for d in data for r in safety_rows(d)];resolution=resolution_rows(data);coeff=coefficient_rows(data,plans);sources=source_rows(data)
    for name,rows in [('deterministic_comparison.csv',det),('safety_constraints.csv',safety),('wetland_contrasts.csv',contrasts),
        ('all20_first_failures.csv',scenarios),('numerical_resolution.csv',resolution),('coefficient_sensitivity_summary.csv',coeff),('source_hashes.csv',sources)]:csv(name,rows)
    nominal=[d for d in data if d['_group']=='nominal'];warm=[d for d in nominal if not deterministic(d) and d['dt_s']==120 and d.get('ode_started')]
    cold=[d for d in nominal if not deterministic(d) and d['dt_s']==120 and not d.get('ode_started')]
    worst=max((abs(r['first_failure_time_difference_s']) for r in resolution if r['fine_dt_s']==30 and r['first_failure_time_difference_s'] is not None),default=0.)
    worst_coarse=max((abs(r['first_failure_time_difference_s']) for r in resolution if r['fine_dt_s']==60 and r['first_failure_time_difference_s'] is not None),default=0.)
    completed=[d for d in data if valid(d)];fine=[r for r in resolution if r['fine_dt_s']==30 and r['both_complete']]
    reasoncounts=dict(Counter(d.get('failure_reason') for d in warm))
    summary=dict(report_scope='Saved-result reporting only; no ODE, fitting, optimization, data generation or parameter selection',
        distinct_original_weather_scenarios=20,warm_original_weather_scenarios=16,cold_original_weather_scenarios=4,
        revised_case_evaluations=len(data),revised_completed_trajectories=len(completed),
        evaluation_count_scope='Formal nominal/low/high revised matrix only; auxiliary initial smoke cases and original-B control runs are separate',
        auxiliary_initial_smoke_case_evaluations=len(list((RESULTS/'first_pair').glob('*.json'))),
        separate_original_B_control_case_evaluations=len(list((RESULTS/'original_regression').glob('*__*.json'))),
        nominal_120s_attempted_warm_scenario_pairs=len({d['scenario_id'] for d in warm}),
        nominal_120s_cold_scenario_ids_accounted=len({d['scenario_id'] for d in cold}),
        nominal_120s_warm_case_evaluations=len(warm),nominal_120s_warm_completed=sum(valid(d) for d in warm),
        nominal_120s_cold_case_exclusions=len(cold),nominal_warm_failure_reasons=reasoncounts,
        all_nominal_warm_first_failures_are_crest_domain=bool(warm) and set(reasoncounts)=={'broad_crest_head_length_outside_supported_range'},
        complete_nominal_weather_pairs_120s=sum(all(valid(d) for d in nominal if d['scenario_id']==sid and d['dt_s']==120) for sid in {d['scenario_id'] for d in warm+cold}),
        maximum_complete_balance_relative=max(d['balance_relative'] for d in completed),
        maximum_successful_RHS_mass_identity_m3_s=max(d['all_rhs_observations']['derivative_physical_mass_identity_max_abs_m3_s'] for d in data if d['all_rhs_observations']['derivative_physical_mass_identity_max_abs_m3_s'] is not None),
        max_nominal_60_to_30_common_grid_river_difference_m=max(r['max_common_grid_river_difference_m'] for r in fine),
        max_nominal_60_to_30_common_grid_all_stage_difference_m=max(r['max_common_grid_all_stage_difference_m'] for r in fine),
        max_nominal_60_to_30_native_grid_peak_difference_m=max(abs(r['native_grid_river_peak_difference_m']) for r in fine),
        max_nominal_120_to_60_failure_detection_shift_s=worst_coarse,
        max_nominal_60_to_30_failure_detection_shift_s=worst,
        all_resolution_classifications_match=all(r['same_classification'] for r in resolution),
        original_source_hashes_current=all(r['hash_matches_current'] for r in sources if not r['source'].startswith('D/gate_hydraulics_revision/')),
        revised_source_hashes_current=all(r['hash_matches_current'] for r in sources if r['source'].startswith('D/gate_hydraulics_revision/')),
        hydraulic_coverage_resolved=False,ready_for_bounded_D2=False,
        reason='Warm original inputs enter an unsupported crest head/length range; nominal deterministic dam 1 also exceeds its unchanged safety limit. No extreme/cold scenario was dropped or modified.')
    write(RESULTS/'summary.json',summary)
    plot_deterministic(data);plot_failures(scenarios);plot_resolution_sensitivity(data,coeff)
    n30=[dict(r,design_display='Historical B-selected wetland design' if r['case']=='wetland' and r['design']=='best' else 'A/A reference wetlands' if r['case']=='wetland' else 'No wetland') for r in det if r['group']=='nominal' and r['dt_s']==30]
    first=next(r for r in scenarios if r['group']=='nominal' and r['dt_s']==30 and r['scenario_id']=='copula_var_00000' and r['case']=='wetland')
    worstmargin=min(r['H_d1_margin_m'] for r in n30)
    readme=f'''# Stage D revised-hydraulic execution results

This report reads saved execution artifacts only. The revised nominal model is **not ready for D.2 uncertainty propagation**. It completes the original deterministic event, but dam 1 exceeds its unchanged safety limit by approximately {-worstmargin:.6f} m. All 16 original warm weather scenarios still stop in both models, now at the explicit broad-crest head/length applicability boundary. The four cold scenarios remain preflight exclusions.

The implementation and its physical assumptions are described separately in the Stage D hydraulic-method documentation. This report does not validate those assumptions merely because a numerical trajectory completes. No gates, weather histories, initial states, safety limits or other fixed parameters were changed during reporting.

## Original B versus revised D

Nominal coefficient Cw = 2.8 ft½/s ({2.8*np.sqrt(.3048):.8f} m½/s). The table uses completed 30 s revised runs and the preserved 120 s B reference. Source arrays are compared on exactly aligned original 120 s nodes. Revised native-grid peaks also include finer samples; these are structural-model comparisons, not unchanged-model regression tests.

{markdown_table(n30,[('case','Case'),('design_display','Design'),('original_B_peak_Hr_m','Original B river peak (m)'),('revised_D_peak_Hr_m','Revised D river peak (m)'),('H_d1_margin_m','Dam 1 margin (m)'),('safety_feasible','Safety feasible')])}

The unchanged original implementation was separately rerun in [original_regression](results/original_regression/); its exact regression results must be distinguished from the revised trajectories. A reduced river peak under a dam-safety violation is not an admissible optimization gain. No Stage D optimization was performed, and the saved B/B label means the design selected under Stage B physics, not an optimum under the revised hydraulics.

![Original and revised deterministic hydrographs](plots/deterministic_comparison.png)

*Same event, gates, physical initial states and saved wetland design. The revised nominal curves complete numerically but exceed dam 1's unchanged safety limit. The figure uses 30 s D output and preserved 120 s B output; it does not imply calibration.*

[Deterministic comparison](tables/deterministic_comparison.csv) records all coefficients/resolutions, source hashes, river/outlet peaks and safety margins. [Wetland contrasts](tables/wetland_contrasts.csv) reports no-wetland minus wetland peaks under identical controls. Those structural contrasts retain the extra 5,500 m³ initial wetland water and flag unsafe pairs explicitly; they are not new optimization benefits.

## Every original weather scenario

Each resolution evaluates all 20 original IDs in each physical model: 16 warm trajectories are attempted and four cold inputs are excluded before ODE entry. At each nominal resolution, **0/16 attempted warm scenario pairs complete**. The four cold scenarios are separately retained as preflight exclusions, accounting for all 20 original IDs. The low and high coefficient sensitivity runs at 120 s also have zero complete warm pairs. No invalid trajectory receives a complete-event peak, discharge or safety result.

For `copula_var_00000`, the nominal 30 s wetland run first fails at {first['first_failure_time_s']:.0f} s ({first['first_failure_time_s']/3600:.6f} h), with H = {first['H_m']:.9f} m, a = {first['a_m']:.9f} m, and H/a = {first['H_over_a']:.6f}. The actual reason is `{first['failure_reason']}`. With the declared 1 m crest, H/L falls below 0.05. This is a bounded physical applicability exclusion, not an ODE numerical failure and not the old H/a >= 1.25 rejection. The old G5 residual remains descriptive only.

[All-scenario first failures](tables/all20_first_failures.csv) retains every case, coefficient and timestep; original seeds/weather hashes; UTC and elapsed first-failure time; H, a, H/a, tailwater; the new physical residual; and the historical 120 s G5 failure time. Later interruption does not count as successful uncertainty propagation. Cold rows identify the first excluded weather interval rather than an ODE failure time. No precipitation/temperature was clipped, no weather was regenerated, and no scenario was excluded for flood magnitude or negative wetland benefit.

![First failures and cold exclusions](plots/scenario_first_failures.png)

*All 20 original IDs remain visible. Warm-scenario points show first evaluated failing RK states, not exact continuous-time crossings. The nominal 30 s and original 120 s failures correspond to different physical applicability rules. Cold inputs are shown separately because their ODEs never start.*

## Resolution and coefficient sensitivity

The original weather and forcing conventions are unchanged: stochastic rain/radiation remain constant within each saved hour with official solver restarts at hourly boundaries; deterministic weather retains its original linear knots. All nominal outcomes agree across 120, 60 and 30 s. The largest complete deterministic 60-to-30 s common-grid river difference is {summary['max_nominal_60_to_30_common_grid_river_difference_m']:.6g} m; the corresponding largest difference over all reported stages is {summary['max_nominal_60_to_30_common_grid_all_stage_difference_m']:.6g} m. The largest native-grid river-peak difference is {summary['max_nominal_60_to_30_native_grid_peak_difference_m']:.6g} m ({1000*summary['max_nominal_60_to_30_native_grid_peak_difference_m']:.6g} mm): the finer output grid samples a peak between coarser nodes, so this need not equal the common-grid error. The maximum change in first-failure detection time is {worst_coarse:.0f} s from 120 to 60 s, and {worst:.0f} s from 60 to 30 s. No full-trajectory convergence claim is made for interrupted weather runs.

{markdown_table(coeff,[('group','Sensitivity member'),('Cw_US','Cw (ft½/s)'),('deterministic_safe_cases','Safe deterministic cases / 3'),('warm_case_completions','Complete warm cases / 32'),('cold_weather_IDs','Cold IDs retained')])}

The coefficient values were declared before these runs; the high value is not selected or recalibrated to obtain safety. Even though its deterministic cases meet the modeled limits, all 16 warm weather scenarios still fail the crest-domain check. [Resolution comparisons](tables/numerical_resolution.csv) separate complete common-grid errors from failure-time changes. [Coefficient summary](tables/coefficient_sensitivity_summary.csv) records actual outcomes without inventing peaks for invalid runs.

![Resolution and coefficient comparison](plots/resolution_and_coefficient_sensitivity.png)

*Left: finite-grid changes in deterministic physical river peaks. Center: dam 1 margin relative to its unchanged limit. Right: actual deterministic river peaks for the predetermined coefficient set. Nominal and low cases are unsafe; the high case is only a sensitivity member, and it does not repair stochastic applicability.*

## Conservation, safety and provenance

The largest relative cumulative water-balance residual among the {len(completed)} completed revised trajectories is {summary['maximum_complete_balance_relative']:.6g}, below the unchanged 1e-4 diagnostic tolerance. The largest successful-RHS mass-identity error is {summary['maximum_successful_RHS_mass_identity_m3_s']:.6g} m³/s. These checks concern numerical accounting, not empirical validation. Internal transfers cancel once; wetland drainage/overflow remain, and no independent dam spillway or state cap is introduced.

[Safety constraints](tables/safety_constraints.csv) include both dam limits and, where present, both wetland limits at every successful RHS/RK evaluation. Incomplete trajectories report observed partial-window maxima only; their full-event safety remains unknown. [Source hashes](tables/source_hashes.csv) compare recorded dependencies with current files. Each execution plan retains the original settings, seed sequence, time window and weather hashes. The preservation manifest protects A/B/C and previous D studies.

The **formal revised matrix contains {len(data)} case evaluations**: nominal 3 resolutions × 43 cases plus two coefficient sensitivities × 43 cases. This denominator excludes {summary['auxiliary_initial_smoke_case_evaluations']} initial smoke-case evaluations and {summary['separate_original_B_control_case_evaluations']} original-B deterministic controls, which are retained separately. The repeated resolutions, coefficients and controls are not independent weather draws. There are exactly 20 original weather scenarios. No Monte Carlo probability, failure probability, polynomial response sensitivity fit or globally feasible design is inferred. A physical treatment covering the encountered crest-domain states, with its geometry and validity justified, remains necessary before bounded D.2 propagation can proceed.
'''
    (HERE/'README_results.md').write_text(readme)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
