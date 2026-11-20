"""Match reported CSV values to saved per-run simulation JSONs; no ODE/NPZ."""
from pathlib import Path
from collections import Counter
import hashlib
import json
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];D2=ROOT/'D/d2_analysis'
inputs={};counts=Counter();errors=[];records=[]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    inputs[str(p.relative_to(ROOT))]=sha(p)
    return json.loads(p.read_text()) if p.suffix=='.json' else pd.read_csv(p,float_precision='round_trip')
def ck(group,condition,source,field):
    counts[group]+=1
    if not condition:errors.append(dict(group=group,source=str(source),field=field))
def near(group,a,b,source,field,tol=1e-9):
    ck(group,bool(np.allclose(a,b,rtol=0,atol=tol)),source,field)

weather=read(D2/'results/metrics.csv')
for _,r in weather.iterrows():
    path=ROOT/r.source;d=read(path)
    ck('weather_source',sha(path)==r.source_sha256,r.source,'source_sha256')
    ck('weather_classification',bool(d.get('full_valid_trajectory_available',False))==bool(r.valid),r.source,'valid')
    if r.valid:
        near('weather_metrics',r.river_peak_m,d['physical_river_peak_m'],r.source,'peak_H')
        near('weather_metrics',r.outlet_peak_m3_s,d['physical_river_peak_discharge_m3_s'],r.source,'peak_Q')
        near('weather_metrics',r.river_peak_time_h*3600,d['physical_river_peak_time_s'],r.source,'peak_time',1e-7)
        near('weather_metrics',r.balance_relative,d['balance_relative'],r.source,'mass_relative',1e-15)
        ck('weather_safety',bool(r.safety_feasible)==(not d['safety_violation']),r.source,'safety')
    records.append(dict(study='weather_all_resolutions_and_controls',source=r.source,valid=bool(r.valid)))

folders=[D2/'parameter_uncertainty/results'/p for p in ('sensitivity','monte_carlo')]
folders += sorted((D2/'historical_storm_robustness/results/ensemble').iterdir())
for folder in folders:
    csvpath=folder/'paired_results.csv'
    if not csvpath.is_file():continue
    frame=read(csvpath)
    for _,r in frame.iterrows():
        for case in ('wetland','no_wetland'):
            path=folder/f'{r.sample_id}__{case}__15s.json';d=read(path)
            ck('run_validity',bool(r[case+'_valid'])==bool(d['full_valid_trajectory_available']),path,'valid')
            ck('run_safety',bool(r[case+'_unsafe'])==bool(d['safety_violation']),path,'unsafe')
            ck('run_completion',bool(d['full_integration_horizon_completed']),path,'completed')
            for column,field in [('peak_stage_m','physical_river_peak_m'),('peak_discharge_m3_s','physical_river_peak_discharge_m3_s'),('peak_stage_time_s','physical_river_peak_time_s'),('balance_relative','balance_relative')]:
                near('run_physical_metrics',r[case+'_'+column],d[field],path,field)
            safety=d['all_rhs_observations']['safety_successful_rhs']
            near('run_safety',r[case+'_minimum_safety_margin_m'],min(v['margin_m'] for v in safety),path,'minimum_safety_margin_m')
            rec=d['parameter_evaluation']['effective_physical_values']
            near('run_parameters',[r.K_u,r.k_ud_1,r.k_ud_2,r.k_ur,r.R_r],[rec['K_u'],*rec['k_ud'],rec['k_ur'],rec['R_r']],path,'actual_parameters',0)
            initial=sum(d['initial_storage_m3']);final=sum(d['last_completed_state_m3'])
            for column,val in [('initial_total_storage_m3',initial),('final_total_storage_m3',final)]:
                k=case+'_'+column
                if k in r and pd.notna(r[k]):near('run_storage',r[k],val,path,column,1e-8)
            records.append(dict(study=str(folder.relative_to(D2)),source=str(path.relative_to(ROOT)),valid=bool(d['full_valid_trajectory_available'])))

frame=read(D2/'antecedent_conditions/tables/model_window_metrics.csv')
for _,r in frame.iterrows():
    path=ROOT/r.source;d=read(path)
    ck('antecedent_source',sha(path)==r.source_sha256,path,'source_hash')
    ck('antecedent_validity',bool(r.valid)==bool(d['full_valid_trajectory_available']),path,'valid')
    if r.scope=='full_horizon':
        for col,field in [('peak_stage_m','physical_river_peak_m'),('peak_discharge_m3_s','physical_river_peak_discharge_m3_s')]:near('antecedent_full_peaks',r[col],d[field],path,field)
        near('antecedent_full_peaks',r.peak_stage_time_from_window_start_h*3600,d['physical_river_peak_time_s'],path,'peak_time',1e-7)
        near('antecedent_storage',r.initial_total_storage_m3,sum(d['initial_storage_m3']),path,'initial_total',1e-8)
        near('antecedent_storage',r.final_total_storage_m3,sum(d['last_completed_state_m3']),path,'final_total',1e-8)
        records.append(dict(study='antecedent_primary_and_reused_baselines',source=str(path.relative_to(ROOT)),valid=bool(r.valid)))
    elif r.antecedent_hours>0:
        phase=d['phase_diagnostics']['event' if r.scope=='event_window' else 'antecedent']
        for col,field in [('outlet_volume_m3','outlet_volume_m3'),('external_net_volume_m3','external_net_volume_m3'),('end_mass_residual_m3','final_balance_residual_m3')]:near('antecedent_direct_phase_flux',r[col],phase[field],path,col,1e-8)
        near('antecedent_phase_storage',r.initial_total_storage_m3,sum(phase['initial_storage_m3']),path,'initial',1e-8)
        near('antecedent_phase_storage',r.final_total_storage_m3,sum(phase['final_storage_m3']),path,'final',1e-8)
        ck('antecedent_phase_status',bool(phase['complete']) and bool(phase['balance_passed']) and bool(phase['forcing_conservation_passed']),path,'phase_status')
        ck('antecedent_phase_status',bool(r.phase_rhs_unsafe)==bool(phase['safety_exceeded_at_successful_rhs']),path,'phase_safety')

report=dict(status='pass' if not errors else 'fail',checks=sum(counts.values()),checks_by_group=dict(counts),errors=errors,
    unique_saved_input_files=len(inputs),reported_record_references=len(records),unique_run_JSONs=len({r['source'] for r in records}),
    record_roles=dict(Counter(r['study'] for r in records)),source_sha256=inputs,script_sha256=sha(Path(__file__)),
    limitations='Repeated JSONs/roles are audit evidence, not new or independent simulations. Event-only peak extraction and dense controller continuity were previously audited from full arrays; they are not reconstructed here.',
    ODE_integrations=0)
(HERE/'diagnostics_crosscheck.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))
raise SystemExit(bool(errors))
