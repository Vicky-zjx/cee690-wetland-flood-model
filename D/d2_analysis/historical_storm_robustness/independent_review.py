"""Recheck completed historical-storm saved outputs; writes nothing and prints audit JSON.
Run from any directory with the existing Python environment, using -B.
"""

from pathlib import Path
import json, hashlib
import os
from datetime import datetime, timezone
ROOT = Path(__file__).resolve().parents[3]
os.chdir(ROOT)
import numpy as np
import pandas as pd
here=Path('D/d2_analysis/historical_storm_robustness')
plan=json.loads((here/'study_plan.json').read_text())
points={p['sample_id']:p for p in plan['points']}
expected_hash=hashlib.sha256(Path('D/feedback_control/policy.json').read_bytes()).hexdigest()
limits={'H_d1':101.2,'H_d2':101.,'H_w1':98.9,'H_w2':98.8}
maxbalance=maxmass=maxnet=maxrain=maxsolar=maxgates=0.
minmargin=minratio=np.inf
counts={};issues=[];allclass={};initial_counts={}
for event in plan['events']:
 if event['benchmark']:continue
 eid=event['event_id'];directory=here/'results/ensemble'/eid
 csv=directory/'paired_results.csv'
 if not csv.exists():counts[eid]=0;continue
 frame=pd.read_csv(csv,float_precision='round_trip')
 counts[eid]=len(frame)
 assert len(frame)==64 and set(frame.sample_id)==set(points), (eid,'missing or duplicate sample IDs')
 initial_counts[eid]={c:0 for c in ('wetland','no_wetland')}
 weather=pd.read_csv(event['forcing_file'])
 for row in frame.itertuples():
  sid=row.sample_id
  assert sid in points
  point=points[sid];dcs={};ncs={}
  for case in ('wetland','no_wetland'):
   name=f'{sid}__{case}__15s'
   d=json.loads((directory/(name+'.json')).read_text());dcs[case]=d
   cls=d['classification'];allclass[cls]=allclass.get(cls,0)+1
   record=d['parameter_evaluation'];fc=d['fixed_configuration']
   assert record==fc['parameter_evaluation'],(eid,sid,case,'nested record')
   assert record['effective_physical_values']==point['physical_values'],(eid,sid,case,'physical')
   assert d['historical_storm_study']['sample_physical_values']==point['physical_values']
   assert record['policy_sha256']==expected_hash
   assert record['prepared_forcing_sha256']==event['forcing_sha256']==d['forcing_file_sha256']
   assert pd.Timestamp(d['event_start_utc'])==pd.Timestamp(event['start_utc'])
   assert (pd.Timestamp(event['end_utc'])-pd.Timestamp(event['start_utc'])).total_seconds()==518400.
   assert d['control_period_s']==120.
   assert record['boundaries']=={'Q_ext_m3_s':0.,'h_d_m':97.15,'h_g_m':97.15,'B':0.}
   assert record['initial_actual_opening_m']==[.0015,.0015]
   for key,value in fc['stage_d_feedback_source_sha256'].items():
    assert hashlib.sha256(Path(key).read_bytes()).hexdigest()==value
   maxbalance=max(maxbalance,d.get('balance_relative',0.))
   maxmass=max(maxmass,d['all_rhs_observations']['derivative_physical_mass_identity_max_abs_m3_s'])
   cons=d['forcing_conservation']
   assert cons['pass_conservation']
   assert abs(cons['source_precipitation_m']-weather.precipitation_mm_per_h.sum()/1000)<1e-12
   assert abs(cons['source_solar_Wh_m2']-weather.solar_W_per_m2.sum())<1e-7
   maxrain=max(maxrain,abs(cons['precipitation_error_m']))
   maxsolar=max(maxsolar,abs(cons['solar_error_Wh_m2']))
   safety=d['safety']+d['all_rhs_observations']['safety_successful_rhs']
   for item in safety:
    assert item['threshold_m']==limits[item['variable']]
    minmargin=min(minmargin,item['margin_m'])
   ratios=[a['minimum_open_gate_H_over_a'] for a in d['actual_gate_successful_rhs_observations'] if a['minimum_open_gate_H_over_a'] is not None]
   minratio=min(minratio,*ratios)
   if not d['full_valid_trajectory_available']:
    issues.append((eid,sid,case,cls));continue
   assert d['full_integration_horizon_completed']
   assert d['balance_relative']<=d['balance_relative_tolerance']==1e-4
   with np.load(directory/(name+'.npz')) as z:
    z={k:z[k].copy() for k in z.files}
   ncs[case]=z
   assert z['time'][0]==0 and z['time'][-1]==518400. and np.all(np.diff(z['time'])==15.)
   assert np.array_equal(z['x'][:,0],d['initial_storage_m3'])
   assert float(z['downstream_stage'].max())==d['physical_river_peak_m']
   assert float(z['time'][np.argmax(z['downstream_stage'])])==d['physical_river_peak_time_s']
   assert float(z['downstream_discharge'].max())==d['physical_river_peak_discharge_m3_s']
   assert float(z['downstream_stage'].max())==getattr(row,case+'_peak_stage_m')
   q=np.minimum(3.,np.maximum(z['downstream_stage']-97.15,0)/point['physical_values']['R_r'])
   assert np.allclose(q,z['downstream_discharge'],rtol=0,atol=1e-12)
   assert np.allclose(z['actual_opening_m'][:,0],[.0015,.0015],rtol=0,atol=1e-16)
   assert np.max(abs(np.diff(z['actual_opening_m'],axis=1)))/15 <= .0005+1e-15
   assert np.min(z['actual_opening_m'])>=0 and np.max(z['actual_opening_m'])<=.3+1e-15
   assert abs(np.max(abs(z['balance']))-d['balance_max_m3'])<1e-8
   net=z['x'].sum(axis=0)-z['x'][:,0].sum()-z['balance']
   assert np.array_equal(net,z['external_net_integrated_m3'])
   maxnet=max(maxnet,abs(net[-1]-d['completed_blocks_mass_diagnostics']['external_net_volume_m3']))
   assert abs(net[-1]-getattr(row,case+'_external_net_integrated_m3'))<1e-8
   for where in ('node','rhs'):
    for item in safety:
     assert item['threshold_m']==limits[item['variable']]
   initial_counts[eid][case]+=int(d['physical_river_peak_time_s']==0.)
  if set(ncs)=={'wetland','no_wetland'}:
   assert dcs['wetland']['forcing_array_sha256']==dcs['no_wetland']['forcing_array_sha256']
   assert np.array_equal(ncs['wetland']['x'][[0,1,4,5,6,7],0],ncs['no_wetland']['x'][:,0])
   gates=float(np.max(abs(ncs['wetland']['actual_opening_m']-ncs['no_wetland']['actual_opening_m'])))
   maxgates=max(maxgates,gates)
   assert abs(gates-row.maximum_paired_gate_opening_difference_m)<1e-14

# Reused benchmark outputs retain every original value, not merely rounded summaries.
old=Path('D/d2_analysis/parameter_uncertainty/results/monte_carlo')
bench=here/'results/ensemble/benchmark_september_2022'
for name in ('paired_results.csv','model_results.csv'):
 a=pd.read_csv(old/name,float_precision='round_trip')
 b=pd.read_csv(bench/name,float_precision='round_trip')
 pd.testing.assert_frame_equal(a,b[a.columns],check_exact=True)
for p in old.glob('mc_*__*__15s.json'):
 a=json.loads(p.read_text());b=json.loads((bench/p.name).read_text())
 b.pop('historical_event');b.pop('historical_benchmark_reuse')
 assert a==b,p.name
source_plan=Path('D/d2_analysis/parameter_uncertainty/results/study_plan.json')
assert plan['points']==json.loads(source_plan.read_text())['mc']
assert plan['parameter_sampling_calls']==0
import collections
branch_margins={};new_branch_margins={};ratios={};regimes=collections.Counter()
classes=collections.Counter();events_summary=[];input_sources=[];allmaxbalance=0.;allmaxgate=0.;allmaxmass=0.;latch_mismatches=0
for event in plan['events']:
 eid=event['event_id'];directory=here/'results/ensemble'/eid;path=directory/'paired_results.csv'
 f=pd.read_csv(path,float_precision='round_trip')
 assert len(f)==64 and set(f.sample_id)==set(points)
 input_sources.append(dict(event_id=eid,paired_results_path=str(path),paired_results_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),forcing_sha256=event['forcing_sha256']))
 latch_mismatches+=int(f.paired_latch_mismatch_count.sum())
 allmaxgate=max(allmaxgate,float(f.maximum_paired_gate_opening_difference_m.max()))
 vals=f.delta_peak_stage_m
 events_summary.append(dict(event_id=eid,planned_pairs=64,valid_pairs=int(f.paired_valid.sum()),valid_unsafe_wetland=int(f.wetland_unsafe.sum()),valid_unsafe_no_wetland=int(f.no_wetland_unsafe.sum()),initial_peak_wetland=int((f.wetland_peak_stage_time_s==0).sum()),initial_peak_no_wetland=int((f.no_wetland_peak_stage_time_s==0).sum()),delta_peak_stage_positive=int((vals>0).sum()),delta_peak_stage_zero=int((vals==0).sum()),delta_peak_stage_negative=int((vals<0).sum()),median_delta_peak_stage_m=float(vals.median())))
 for sid in points:
  for case in ('wetland','no_wetland'):
   d=json.loads((directory/f'{sid}__{case}__15s.json').read_text());tag=[eid,sid,case]
   assert d['parameter_evaluation']['effective_physical_values']==points[sid]['physical_values']
   classes[d['classification']]+=1
   regimes.update(d['hydraulic_diagnostics']['regime_rhs_counts'])
   allmaxbalance=max(allmaxbalance,d['balance_relative'])
   allmaxmass=max(allmaxmass,d['all_rhs_observations']['derivative_physical_mass_identity_max_abs_m3_s'])
   for item in d['safety']+d['all_rhs_observations']['safety_successful_rhs']:
    key=case+'/'+item['variable']
    assert item['threshold_m']==limits[item['variable']]
    value=dict(margin_m=item['margin_m'],peak_m=item['peak_m'],threshold_m=item['threshold_m'],source_case=tag)
    if key not in branch_margins or item['margin_m']<branch_margins[key]['margin_m']:branch_margins[key]=value
    if not event['benchmark'] and (key not in new_branch_margins or item['margin_m']<new_branch_margins[key]['margin_m']):new_branch_margins[key]=value
   for item in d['actual_gate_successful_rhs_observations']:
    key=case+'/gate_'+str(item['branch']);ratio=item['minimum_open_gate_H_over_a']
    if key not in ratios or ratio<ratios[key]['H_over_a']:ratios[key]=dict(H_over_a=ratio,source_case=tag)
preservation_path=here/'results/ensemble/preservation_check.json'
preservation=json.loads(preservation_path.read_text())
assert preservation['passed'] and not preservation['changed'] and not preservation['missing']
regression=json.loads((here/'results/regression/regression_summary.json').read_text())
assert regression['passed'] and all(all(r['exact_arrays'].values()) for r in regression['checks'])
pilot=json.loads((here/'results/pilot/pilot_summary.json').read_text())
assert pilot['passed'] and len(pilot['checks'])==24
def digest(path):
 return hashlib.sha256(Path(path).read_bytes()).hexdigest()
report=dict(
 status='pass',reviewed_at_utc=datetime.now(timezone.utc).isoformat(),
 scope='Read-only saved-output audit; no sampling, fitting, ODE execution or modification of frozen inputs/results',
 model_runs_performed=0,parameter_sampling_calls=0,distinct_parameter_vectors=64,
 planned_event_parameter_pairs=320,retained_model_outcomes=640,
 new_full_payload_pairs_audited=counts,classifications=dict(classes),event_accounting=events_summary,
 original_benchmark_csv_values_exact=True,original_benchmark_json_fields_exact=True,
 original_2022_regression=regression,numerical_pilot_comparisons=24,numerical_pilot_passed=True,
 minimum_safety_margin_by_case_and_branch=branch_margins,
 new_event_minimum_safety_margin_by_case_and_branch=new_branch_margins,
 minimum_open_gate_ratio_by_case_and_branch=ratios,regime_rhs_counts=dict(regimes),
 hydraulic_coverage='Only closed and energy_gate regimes occurred; this experiment does not validate the weir/transition closure',
 maximum_relative_mass_residual=allmaxbalance,unchanged_relative_mass_tolerance=1e-4,
 maximum_rhs_mass_identity_error_m3_s=allmaxmass,
 maximum_new_block_net_recovery_error_m3=maxnet,
 external_net_interpretation='Recovered from stored block-integrated balance, compared to independently saved block diagnostics; not a separate re-integration of concatenated discontinuous forcing',
 maximum_new_rain_integral_error_m=maxrain,maximum_new_solar_integral_error_Wh_m2=maxsolar,
 maximum_paired_actual_gate_opening_difference_m=allmaxgate,paired_latch_mismatches=latch_mismatches,
 findings_review=dict(path=str(here/'FINDINGS.md'),sha256=digest(here/'FINDINGS.md'),
   result='No material scientific or quantitative discrepancy; event counts, means, signs, duration/volume caveats, matched effects, correlations and conditional interpretation checked against saved results',
   source_precision='round_trip parsing of source result CSVs and direct saved JSON/NPZ values'),
 inputs=dict(study_plan_sha256=digest(here/'study_plan.json'),original_parameter_plan_sha256=digest(source_plan),
   event_manifest_sha256=digest(here/'events/manifest.json'),paired_results=input_sources),
 preservation_marker_sha256=digest(preservation_path),preservation=preservation,
 limitations=['Audit verifies saved trajectories and declared tests, not real-system calibration or universal safety.',
 'No downstream flood threshold is available.',
 'The same 64 parameter vectors are reused across fixed events; the 320 pairs are not independent weather draws.',
 'Shared initial states are unchanged; the wetland model retains its accepted additional 5500 m3 initial water.'])
print(json.dumps(report,indent=2,allow_nan=False))

