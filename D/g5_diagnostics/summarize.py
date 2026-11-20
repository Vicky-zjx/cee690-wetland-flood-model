"""Summarize only the saved, bounded Stage D diagnostic runs; never simulate."""
from pathlib import Path
import hashlib
import json
import os
import sys
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/cee690-stage-d-matplotlib')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import runner

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def write_json(path,value):
    path.write_text(json.dumps(runner.json_clean(value),indent=2,allow_nan=False)+'\n')

def main():
    folders={f'original_reference_{dt}s':3 for dt in (120,60,30)}
    folders.update({'original_20_120s':40,'original_pilot_60s':20,'original_pilot_30s':20,
                    'schedule_scale_025_reference_120s':3,'schedule_scale_050_reference_120s':3,'schedule_scale_025_pilot_120s':20})
    folders.update({f'constant_005_{kind}_{dt}s':3 if kind=='reference' else 20 for kind in ('reference','pilot') for dt in (120,60,30)})
    frames=[]; reports={}; safety=[]
    for name,count in folders.items():
        directory=HERE/'results'/name
        frame=pd.read_csv(directory/'scenario_diagnostics.csv');assert len(frame)==count,(name,len(frame),count)
        assert (directory/'counts.json').is_file(),name
        frame.insert(0,'experiment',name)
        frame['forcing_conservation_pass']=frame['forcing_conservation_pass'].astype(object)
        if 'schedule_kind' not in frame:frame['schedule_kind']=np.where(frame.gate_scale==1,'original','uniform_scale')
        if 'constant_opening' not in frame:frame['constant_opening']=np.nan
        for i,r in frame.iterrows():
            basename=f'{r.scenario_id}__{r.case_name}__{r.design_source_label}'
            report=json.loads((directory/(basename+'.json')).read_text())
            reports[(name,r.scenario_id,r.case_name,r.design_source_label)]=report
            v=report.get('first_violation') or {}; obs=report['all_rhs_observations']
            frame.loc[i,'safety_peak_scope']='complete trajectory, all successful RK stages' if r.valid else 'partial successful RK stages only; full safety unknown'
            frame.loc[i,'time_s_meaning']='preflight offending weather interval, no ODE started' if not r.ode_started else ('first detected RHS exception, not exact continuous crossing' if not r.valid else 'no failure')
            fc=report.get('forcing_conservation',{})
            if 'mapped_precipitation_m' in fc:
                pe=fc['mapped_precipitation_m']-fc['source_precipitation_m'];se=fc['mapped_solar_Wh_m2']-fc['source_solar_Wh_m2']
                frame.loc[i,'rain_mapping_error_m']=pe;frame.loc[i,'solar_mapping_error_Wh_m2']=se
                frame.loc[i,'forcing_conservation_pass']=abs(pe)<=1e-12 and abs(se)<=1e-7
            loc=report.get('failure_location') or {}
            frame.loc[i,'failure_file']=loc.get('file');frame.loc[i,'failure_line']=loc.get('line')
            frame.loc[i,'failure_function']=loc.get('function')
            frame.loc[i,'first_temperature_C']=v.get('temperature_C')
            frame.loc[i,'first_failed_RK_evaluation']=(obs.get('first_rhs_exception') or {}).get('rk_evaluation')
            for record in obs['safety_successful_rhs']:
                safety.append(dict(experiment=name,scenario_id=r.scenario_id,case_name=r.case_name,design=r.design_source_label,
                                   complete=bool(r.valid),scope='All successful official RHS/RK evaluations; partial if run incomplete',**record))
        frames.append(frame)
    allrows=pd.DataFrame([record for frame in frames for record in frame.to_dict('records')])
    # Older per-run CSVs omit schedule_kind only; the immutable manifest declares the actual schedule.
    allrows.to_csv(HERE/'scenario_diagnostics.csv',index=False)
    pd.DataFrame(safety).to_csv(HERE/'safety_constraints.csv',index=False)
    originals=allrows[allrows.experiment=='original_20_120s']
    historical=pd.read_csv(HERE/'recovered/original_all_case_diagnostics.csv')
    comparison=originals.merge(historical,on='scenario_id',suffixes=('_new','_old'))
    comparison=comparison[comparison.case_name==comparison['case']]
    assert len(comparison)==40
    mapping={'time_s':'failure_time_s','branch':'failure_branch','water_depth_m':'upstream_depth_at_failure_m',
             'opening_height_m':'gate_opening_at_failure_m','residual_m':'G5_residual_at_failure_m'}
    repro=[]
    for _,r in comparison.iterrows():
        row=dict(scenario_id=r.scenario_id,case_name=r.case_name,forcing_hash_match=r.forcing_sha256_new==r.forcing_sha256_old,
                 reason_match=r.failure_reason_new==r.failure_reason_old)
        for new,old in mapping.items():row[new+'_difference']=None if pd.isna(r[new]) or pd.isna(r[old]) else float(r[new]-r[old])
        repro.append(row)
    pd.DataFrame(repro).to_csv(HERE/'original_reproduction_check.csv',index=False)
    # Input availability, identities, units, conservation and retained cold windows.
    inputrows=[]
    for path in sorted((HERE/'recovered/scenarios').glob('*.csv')):
        frame=pd.read_csv(path);meta=json.loads(path.with_suffix('.json').read_text());times=pd.to_datetime(frame.timestamp_utc,utc=True)
        d=reports[('original_20_120s',path.stem,'wetland','best')]
        inputrows.append(dict(scenario_id=path.stem,seed_sequence=json.dumps(meta['seed_sequence_entropy']),
            forcing_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),saved_hash_match=hashlib.sha256(path.read_bytes()).hexdigest()==meta['forcing_sha256'],
            rows=len(frame),start_utc=times.iloc[0].isoformat(),end_exclusive_utc=(times.iloc[-1]+pd.Timedelta(hours=1)).isoformat(),
            hourly_utc=bool(np.all(np.diff(times.astype('int64'))==3600*10**9)),rainfall_mm=frame.precipitation_mm_per_h.sum(),
            minimum_temperature_C=frame.temperature_C.min(),nonpositive_temperature_hours=int((frame.temperature_C<=0).sum()),
            missing_values=int(frame.isna().sum().sum()),negative_rain_solar_wind=bool((frame[['precipitation_mm_per_h','solar_W_per_m2','wind_m_per_s']]<0).any().any()),
            **d['forcing_conservation']))
    pd.DataFrame(inputrows).to_csv(HERE/'input_audit.csv',index=False)
    # Every pilot row retained, including the cold preflight exclusion.
    pairs=[]
    for name in ['original_20_120s','original_pilot_60s','original_pilot_30s','schedule_scale_025_pilot_120s']+[f'constant_005_pilot_{dt}s' for dt in (120,60,30)]:
        f=allrows[allrows.experiment==name]
        if name=='original_20_120s':f=f[f.scenario_id.str[-5:].astype(int)<5]
        for sid,g in f.groupby('scenario_id'):
            assert len(g)==2
            pairs.append(dict(experiment=name,scenario_id=sid,dt_s=int(g.dt_s.iloc[0]),
                wetland_status=g[g.case_name=='wetland'].classification.iloc[0],no_wetland_status=g[g.case_name=='no_wetland'].classification.iloc[0],
                complete_pair=bool(g.valid.all()),safety_feasible_pair=bool(g.valid.all() and all(v == False for v in g.rhs_safety_violation)),
                weather_hash_match=g.forcing_sha256.nunique()==1,seed_match=g.original_seed.nunique()==1))
    pd.DataFrame(pairs).to_csv(HERE/'pilot_summary.csv',index=False)
    # All original/constant trajectories and failures at 120,60,30 seconds.
    resolutions=[]
    for kind,scenario_ids in [('reference',['deterministic_september_2022']),('pilot',sorted(set(p['scenario_id'] for p in pairs)))]:
        for schedule in ('original','constant_005'):
            for sid in scenario_ids:
                for case,design in ([('wetland','reference'),('wetland','best'),('no_wetland','best')] if kind=='reference' else [('wetland','best'),('no_wetland','best')]):
                    for coarse,fine in ((120,60),(60,30)):
                        def key(dt):
                            name=f'{schedule}_{kind}_{dt}s'
                            if schedule=='original' and kind=='pilot' and dt==120:name='original_20_120s'
                            return name,sid,case,design
                        kc,kf=key(coarse),key(fine);dc,df=reports[kc],reports[kf]
                        row=dict(schedule=schedule,scope=kind,scenario_id=sid,case=case,design=design,coarse_dt_s=coarse,fine_dt_s=fine,
                                 coarse_status=dc['classification'],fine_status=df['classification'])
                        vc,vf=dc.get('first_violation') or {},df.get('first_violation') or {}
                        if vc.get('time_s') is not None and vf.get('time_s') is not None:
                            row['first_detected_failure_time_difference_s']=vc['time_s']-vf['time_s']
                            row['failure_branch_same']=vc.get('branch')==vf.get('branch')
                        if dc['full_valid_trajectory_available'] and df['full_valid_trajectory_available']:
                            def load(k):return np.load(HERE/'results'/k[0]/f'{sid}__{case}__{design}.npz')
                            with load(kc) as ac,load(kf) as af:
                                ratio=coarse//fine
                                assert np.array_equal(ac['time'],af['time'][::ratio])
                                hr_c=ac['downstream_stage'];hr_f=af['downstream_stage']
                                row.update(peak_river_change_m=float(hr_c.max()-hr_f.max()),
                                    max_common_grid_river_difference_m=float(np.max(np.abs(hr_c-hr_f[::ratio]))),
                                    max_common_grid_state_difference_m3=float(np.max(np.abs(ac['x']-af['x'][:,::ratio]))),
                                    coarse_peak_time_s=float(ac['time'][np.argmax(hr_c)]),fine_peak_time_s=float(af['time'][np.argmax(hr_f)]),
                                    coarse_balance_relative=dc['balance_relative'],fine_balance_relative=df['balance_relative'])
                                names=list(ac['output_names']);levels=[names.index(n) for n in ('H_d1','H_d2')+ (('H_w1','H_w2') if case=='wetland' else ())]
                                row['max_common_grid_safety_level_difference_m']=float(np.max(np.abs(ac['y'][levels]-af['y'][levels,::ratio])))
                        resolutions.append(row)
    pd.DataFrame(resolutions).to_csv(HERE/'numerical_resolution.csv',index=False)
    regression=[]
    for (name,sid,case,design),d in reports.items():
        if name=='original_reference_120s':regression.append(dict(case=case,design=design,**d['regression']))
    write_json(HERE/'regression_results.json',dict(original_120s_regression=regression,
        gate_and_full_rhs_tests=json.loads((HERE/'results/gate_tests/gate_regression_results.json').read_text())))
    before=json.loads((HERE/'preservation_before.json').read_text())
    changed=[rel for rel,h in before.items() if not (ROOT/rel).is_file() or hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()!=h]
    write_json(HERE/'preservation_check.json',dict(files_checked=len(before),changed_or_missing=changed,passed=not changed))
    original_summary=originals.groupby(['case_name','classification']).size()
    complete=allrows[allrows.valid==True]
    summary=dict(total_case_evaluations=len(allrows),original_counts={str(k):int(v) for k,v in original_summary.items()},
        max_rain_mapping_error_m=max(abs(r['precipitation_error_m']) for r in inputrows),
        max_solar_mapping_error_Wh_m2=max(abs(r['solar_error_Wh_m2']) for r in inputrows),
        complete_simulations=len(complete),max_complete_balance_relative=float(complete.balance_relative.max()),
        max_derivative_mass_identity_m3_s=float(allrows.derivative_balance_max_m3s.max()),
        protected_files_unchanged=not changed,hydraulic_G5_resolved=False,
        physical_extension_implemented=False,operating_alternative_accepted=False,
        reason='Missing crest/coefficient/transition specification. Tested common schedules do not meet both reference-event safety and complete pilot applicability.')
    write_json(HERE/'summary.json',summary)
    plot_comparison()
    print(json.dumps(summary,indent=2))

def plot_comparison():
    colors=['#2456A6','#858585','#D28E00','#B63A4E']
    tags=['original_reference_120s','schedule_scale_050_reference_120s','schedule_scale_025_reference_120s','constant_005_reference_120s']
    labels=['Original schedule','0.50 x original','0.25 x original','Constant g = 0.005']
    fig,axs=plt.subplots(3,2,figsize=(11,9),sharex=True,layout='constrained')
    for tag,label,color in zip(tags,labels,colors):
        with np.load(HERE/'results'/tag/'deterministic_september_2022__wetland__reference.npz') as d:
            t=d['time']/3600;names=list(d['output_names'])
            for i in range(2):
                axs[0,i].plot(t,d['u'][i],color=color,label=label,lw=1.6)
                axs[1,i].plot(t,d['y'][names.index(f'H_d{i+1}')],color=color,lw=1.6)
                axs[2,i].plot(t,d['y'][names.index(f'H_w{i+1}')],color=color,lw=1.6)
    for i in range(2):
        axs[0,i].set(title=f'Gate {i+1}',ylabel='Normalized opening g')
        axs[1,i].axhline([101.2,101.0][i],color='black',ls='--',lw=1,label='Unchanged safety limit')
        axs[1,i].set(ylabel=f'Dam pool {i+1} level (m)')
        axs[2,i].axhline([98.9,98.8][i],color='black',ls='--',lw=1)
        axs[2,i].set(ylabel=f'Wetland {i+1} level (m)',xlabel='Hours since 2022-09-09 00:00 UTC',xlim=(0,144))
    for ax in axs.ravel():ax.grid(alpha=.2)
    axs[0,0].legend(fontsize=8);axs[1,0].legend(fontsize=8)
    fig.suptitle('Operating alternatives retain the original hydraulic law\nSeptember 2022 A/A reference; smaller openings cause dam safety exceedance',fontsize=13)
    fig.savefig(HERE/'plots/fixed_schedule_safety.png',dpi=180);plt.close(fig)
    fixed=runner.load_fixed()
    from B.gate_hydraulics import sluice_gate_flow
    gates=runner.adapter().wet.gate_parameters(fixed['model_config'])
    fig,axs=plt.subplots(1,2,figsize=(10.5,4),sharey=True,layout='constrained')
    for i,(ax,g) in enumerate(zip(axs,gates)):
        a=.1;ratio=np.linspace(1.25,3,220)
        q=[sluice_gate_flow(g.gate_sill_z+r*a,g.gate_sill_z-1,a/g.gate_max_opening,g).discharge for r in ratio]
        ax.axvspan(0,1,color='#E4E4E4',label='Weir: unsupported')
        ax.axvspan(1,1.25,color='#F8DDAA',label='Transition: unsupported')
        ax.plot(ratio,q,color='#2456A6',label='Original free sluice law')
        ax.scatter([1.25],[q[0]],color='#2456A6',s=20,zorder=3)
        ax.axvline(1,ls=':',color='#555555');ax.axvline(1.25,ls=':',color='#555555')
        ax.set(title=f'Gate {i+1}: width {g.gate_width:.2f} m',xlabel='Upstream head / opening, H/a',xlim=(0,3),ylim=(0,None))
        ax.grid(alpha=.2)
    axs[0].set_ylabel('Forward discharge (m³/s)');axs[0].legend(fontsize=8,loc='upper left')
    fig.suptitle('Hydraulic coverage remains unchanged\nArithmetic diagnostic slice: a = 0.10 m, tailwater below sill; no low-head curve invented',fontsize=12)
    fig.savefig(HERE/'plots/hydraulic_coverage.png',dpi=180);plt.close(fig)

if __name__=='__main__':main()
