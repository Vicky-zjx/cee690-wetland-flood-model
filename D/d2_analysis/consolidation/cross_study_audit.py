#!/usr/bin/env python3
"""Portable saved-JSON/CSV cross-study audit. No ODE, fitting or sampling.
Writes only this consolidation directory. Full-array checks remain in the
recorded independent audits; this script does not require ignored NPZ files.
"""
from collections import Counter
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import beta

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
BASE=ROOT/'D/d2_analysis'
inputs={};counts=Counter();failures=[]
def digest(p): return sha256(p.read_bytes()).hexdigest()
def record(name):
 p=Path(name);p=p if p.is_absolute() else BASE/p
 inputs[p]=digest(p);return p
def js(name): return json.loads(record(name).read_text())
def csv(name): return pd.read_csv(record(name),float_precision='round_trip')
def check(group,ok,detail=''):
 counts[group]+=1
 if not bool(ok): failures.append({'group':group,'detail':detail})
def near(group,a,b,tol=1e-10,detail=''):
 check(group,np.allclose(np.asarray(a,dtype=float),np.asarray(b,dtype=float),atol=tol,rtol=0,equal_nan=True),detail)
def clean(x):
 if isinstance(x,dict): return {str(k):clean(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [clean(v) for v in x]
 if isinstance(x,np.generic): return clean(x.item())
 if isinstance(x,float) and not np.isfinite(x):return None
 return x

def summary(frame,published,group):
 for _,r in published.iterrows():
  metric={'delta_maximum_rise_above_initial_m':'delta_maximum_rise_m'}.get(r['metric'],r['metric']);check(group,metric in frame,metric)
  if metric not in frame:continue
  y=frame.loc[frame.paired_valid,metric].to_numpy(float)
  vals={'n':len(y),'valid_pairs':len(y),'minimum':min(y),'maximum':max(y),'mean':np.mean(y),'median':np.median(y),'sample_sd':np.std(y,ddof=1),'q05':np.quantile(y,.05),'q95':np.quantile(y,.95),'positive_count':sum(y>0),'zero_count':sum(y==0),'negative_count':sum(y<0)}
  for k,v in vals.items():
   if k in r and pd.notna(r[k]):near(group,v,r[k],5e-10,metric+':'+k)

def pair_checks(f,label):
 check('pair_accounting',f.paired_valid.all(),label)
 for c in ['wetland','no_wetland']:
  check('pair_accounting',f[c+'_valid'].all() and not f[c+'_unsafe'].any(),label+c)
  check('no_flood_threshold',f[c+'_flood_threshold_m'].isna().all(),label+c)
  check('rating',all(f[c+'_peak_discharge_m3_s']<3),label+c)
  near('rating',f[c+'_peak_stage_m'],97.15+f.R_r*f[c+'_peak_discharge_m3_s'],4e-13,label+c)
 for m in ['peak_stage_m','peak_discharge_m3_s','total_outlet_volume_m3']:
  near('paired_sign',f['no_wetland_'+m]-f['wetland_'+m],f['delta_'+m],2e-10,label+m)
 near('rating',f.delta_peak_stage_m,f.R_r*f.delta_peak_discharge_m3_s,4e-13,label)
 check('pair_accounting',f.paired_weather_inputs_exact.all() and all(f.paired_latch_mismatch_count==0),label)
 check('pair_accounting',all(f.maximum_paired_gate_opening_difference_m<1e-12),label)

def effect(f):
 return {'pairs':len(f),'positive':int(sum(f.delta_peak_stage_m>0)),'zero':int(sum(f.delta_peak_stage_m==0)),'negative':int(sum(f.delta_peak_stage_m<0)),'mean_peak_reduction_m':float(f.delta_peak_stage_m.mean()),'median_peak_reduction_m':float(f.delta_peak_stage_m.median()),'minimum_peak_reduction_m':float(f.delta_peak_stage_m.min()),'maximum_peak_reduction_m':float(f.delta_peak_stage_m.max()),'mean_outlet_volume_difference_m3':float(f.delta_total_outlet_volume_m3.mean()),'longer_wetland_initial_reference_duration':int(sum(f.delta_initial_level_exceedance_h<0))}

# Retain all weather IDs and both resolutions; cold exclusions are not completed pairs.
w=csv('results/paired_metrics.csv');w=w[w.method.isin(['copula_var','block_bootstrap'])];ws=js('results/summary.json');w15=w[w.dt_s==15];wv=w15[w15.paired_valid]
for dt in [15,7.5]:
 f=w[w.dt_s==dt];check('weather',len(f)==20 and sum(f.paired_valid)==16,'20 planned / 16 valid')
 check('weather',(f.loc[~f.paired_valid,'wetland_classification']=='meteorological_liquid_water_applicability').all(),'four cold exclusions')
check('weather',sum((wv.wetland_river_peak_time_h==0)&(wv.no_wetland_river_peak_time_h==0))==15,'15 initial peaks')
near('weather',wv.no_wetland_river_peak_m-wv.wetland_river_peak_m,wv.difference_river_peak_m,1e-13)
for _,r in csv('results/summary_statistics.csv').query("dt_s == 15 and population == 'all_warm'").iterrows():
 col=r['series']+'_'+r['metric']
 if col not in wv:continue
 y=wv[col].to_numpy(float)
 for k,v in {'n':len(y),'mean':np.mean(y),'median':np.median(y),'minimum':min(y),'maximum':max(y)}.items():near('weather_summary',v,r[k],5e-10,col+k)

# Exact parameter lineage and published summaries.
pp=js('parameter_uncertainty/results/study_plan.json');hp=js('historical_storm_robustness/study_plan.json');ap=js('antecedent_conditions/study_plan.json')
points={p['sample_id']:p for p in pp['mc']};gpoints={p['sample_id']:p for p in pp['grid']}
check('parameter_lineage',hp['points']==pp['mc'],'same 64 historical vectors')
check('parameter_lineage',hp['parameter_sampling_calls']==0 and ap['parameter_sampling_calls']==0 and ap['weather_generation_calls']==0,'no redraw')
check('priors',pp['config']['grid']['probability_weights'] is None,'factorial no probabilities')
check('priors',pp['config']['monte_carlo']['lower_multiplier']==.5 and pp['config']['monte_carlo']['upper_multiplier']==2,'approved exploratory interval')
for p in pp['mc']+pp['grid']:
 v=p['physical_values'];m=p['m_routing'];n=pp['nominal_physical_values']
 check('parameter_lineage',v['k_ud']==[k*m for k in n['k_ud']] and v['k_ur']==n['k_ur']*m,p['sample_id']+' routing dependence')
keys={json.dumps(p['physical_values'],sort_keys=True) for p in pp['mc']+pp['grid']};check('parameter_lineage',len(keys)==91,'27 grid / 64 MC disjoint')
frames={};parameter={};pub=csv('parameter_uncertainty/tables/outcome_summaries.csv')
for phase,plan in [('sensitivity',gpoints),('monte_carlo',points)]:
 f=csv(f'parameter_uncertainty/results/{phase}/paired_results.csv');frames[phase]=f
 check('parameter_lineage',set(f.sample_id)==set(plan),phase)
 pair_checks(f,phase);summary(f,pub[pub.study==phase],'parameter_summaries');parameter[phase]=effect(f)
 parameter[phase]['reused_ids']=f.loc[f.reused_results,'sample_id'].tolist()
 for _,r in f.iterrows():
  p=plan[r.sample_id]
  near('parameter_lineage',[r.K_u,r.k_ud_1,r.k_ud_2,r.k_ur,r.R_r],[p['physical_values']['K_u'],*p['physical_values']['k_ud'],p['physical_values']['k_ur'],p['physical_values']['R_r']],0,r.sample_id)
pilot=js('parameter_uncertainty/results/pilot/pilot_summary.json');check('numerical_evidence',pilot['passed'] and pilot['planned_points']==4 and pilot['paired_evaluations']==8,'parameter pilot')
check('reuse',set(p['source_sample_id'] for p in pp['pilot'])=={'mc_00000','mc_00001','grid_000','grid_026'},'pilot subset')

# Five matched fixed storms, including exact benchmark copies.
hframes={};historical={};hpub=csv('historical_storm_robustness/tables/outcome_summaries.csv');za=csv('robustness_checks/tables/statistical_zero_atoms.csv')
for e in hp['events']:
 eid=e['event_id'];f=csv(f'historical_storm_robustness/results/ensemble/{eid}/paired_results.csv');hframes[eid]=f
 check('historical',len(f)==64 and set(f.sample_id)==set(points),eid)
 pair_checks(f,eid);summary(f,hpub[hpub.event_id==eid],'historical_summaries');historical[eid]=effect(f)
 for _,r in f.iterrows():
  p=points[r.sample_id]['physical_values'];near('parameter_lineage',[r.K_u,r.k_ud_1,r.k_ud_2,r.k_ur,r.R_r],[p['K_u'],*p['k_ud'],p['k_ur'],p['R_r']],0,eid+r.sample_id)
 both=(f.wetland_peak_stage_time_s==0)&(f.no_wetland_peak_stage_time_s==0)
 check('zero_atoms',np.array_equal(f.delta_peak_stage_m==0,both),eid)
 near('zero_atoms',sum(f.delta_peak_stage_m==0),za.loc[za.event_id==eid,'exact_zero_count'].iloc[0],0,eid)
 if e.get('benchmark'):
  original=frames['monte_carlo'];cols=[c for c in original if c in f and c not in ['elapsed_s','reused_results']]
  try:pd.testing.assert_frame_equal(original[cols],f[cols],check_exact=True)
  except AssertionError as err:check('reuse',False,str(err)[:500])
  else:check('reuse',True,'benchmark exact physical and numerical columns')
check('zero_atoms',sum(v['zero'] for v in historical.values())==87,'total 87')
hr=js('historical_storm_robustness/independent_review.json');check('recorded_array_audits',hr['status']=='pass' and hr['numerical_pilot_passed'],'historical')

# ECDF counts preserve ties. Clopper-Pearson is pointwise at fixed thresholds.
cp=csv('historical_storm_robustness/tables/fixed_grid_pointwise_cdf_intervals.csv.gz');max_cp=0
for (eid,metric),g in cp.groupby(['event_id','metric']):
 y=np.sort(hframes[eid][metric]);n=len(y);k=np.searchsorted(y,g.threshold,side='right')
 near('ecdf',k,g.count_leq,0,eid+metric);near('ecdf',k/n,g.empirical_cdf,0,eid+metric)
 lo=np.where(k==0,0,beta.ppf(.025,k,n-k+1));hi=np.where(k==n,1,beta.ppf(.975,k+1,n-k))
 err=max(float(max(abs(lo-g.pointwise_lower))),float(max(abs(hi-g.pointwise_upper))));max_cp=max(max_cp,err);check('ecdf',err<1e-10,eid+metric)

# All 320 endpoint pairs are present. 258 have retained full waveforms.
sp=csv('robustness_checks/tables/storage_paired_decomposition.csv');sm=csv('robustness_checks/tables/storage_model_balances.csv');se=csv('robustness_checks/tables/storage_event_summary.csv')
near('storage',sp.initial_stock_term_m3,-5500,0)
recon=sp.initial_stock_term_m3+sp.final_stock_term_m3+sp.nonoutlet_external_term_m3+sp.residual_term_m3
near('storage',sp.delta_outlet_no_minus_wet_m3,recon,3e-10,'paired identity')
near('storage',sm.final_total_storage_m3-sm.initial_total_storage_m3-sm.external_net_integrated_m3,sm.final_balance_residual_m3,3e-10,'endpoints')
near('storage',sp.final_stock_term_m3,sp.final_groundwater_stock_difference_m3+sp.final_shared_non_groundwater_stock_difference_m3+sp.wetland_compartments_final_m3,3e-10,'compartment partition')
for eid,g in sp.groupby('event_id'):
 source=hframes[eid].set_index('sample_id').loc[g.sample_id]
 near('storage',g.delta_outlet_no_minus_wet_m3,source.delta_total_outlet_volume_m3,0,eid)
 r=se[se.event_id==eid].iloc[0]
 for col,sc in [('delta_outlet_no_minus_wet_m3','mean_delta_outlet_no_minus_wet_m3'),('final_stock_term_m3','mean_final_stock_term_m3'),('nonoutlet_external_term_m3','mean_nonoutlet_external_term_m3'),('wetland_compartments_final_m3','mean_wetland_compartments_final_m3')]:near('storage_summary',g[col].mean(),r[sc],3e-10,eid+col)
sa=js('robustness_checks/tables/storage_audit.json');ta=js('robustness_checks/tables/trajectory_validation.json')
check('recorded_array_audits',sa['status']=='pass' and sa['endpoint_pairs']==320 and sa['full_trajectory_pairs']==258,'storage coverage')

# Verify independent-uniform density ratios; all other laws remain reweightings.
weights=csv('robustness_checks/tables/statistical_law_weights.csv');wo=csv('robustness_checks/tables/statistical_law_outcomes.csv')
u=weights[weights.law.str.contains('uniform') & ~weights.law.str.contains('loguniform')]
check('weights',len(u)==64,'one uniform alternative')
raw=np.array([np.prod([points[s][c]*np.log(4)/1.5 for c in ['m_K_u','m_routing','m_R_r']]) for s in u.sample_id]);near('weights',raw,u.raw_density_ratio,1e-14);near('weights',raw/raw.sum(),u.normalized_weight,1e-15)
for law,g in weights.groupby('law'):
 near('weights',sum(g.normalized_weight),1,1e-14,law)
 for _,r in wo[wo.law==law].iterrows():
  metric={'delta_maximum_rise_above_initial_m':'delta_maximum_rise_m'}.get(r.metric,r.metric)
  y=hframes[r.event_id].set_index('sample_id').loc[g.sample_id,metric].to_numpy();wgt=g.normalized_weight.to_numpy()
  near('weights',sum(wgt*y),r.self_normalized_weighted_mean,5e-10,law+r.metric)
  near('weights',sum(wgt[y==0]),r.weighted_zero_mass,1e-14,law+r.metric)
wr=js('antecedent_conditions/statistics/weight_summary.json')

# Full horizon and event windows are separate, with continuous state carry-over.
am=csv('antecedent_conditions/tables/model_window_metrics.csv');af=csv('antecedent_conditions/tables/paired_window_metrics.csv');ag=csv('antecedent_conditions/tables/event_conditioning_summary.csv')
aa=js('antecedent_conditions/independent_audit.json');apt=js('antecedent_conditions/results/pilot_summary.json')
check('antecedent',aa['status']=='pass' and aa['primary_physical_pairs']==12 and aa['retained_model_records']==26 and aa['valid_model_records']==26,'12 primary + 1 refinement')
check('antecedent',aa['unsafe_valid_model_records']==0 and aa['invalid_model_records']==0,'valid vs safety')
check('antecedent',apt['passed'] and aa['maximum_state_seam_error_m3']==0,'numerical and seam')
check('antecedent',am.flood_threshold_m.isna().all(),'no flood threshold')
for _,r in af.iterrows():
 g=am[(am.event_id==r.event_id)&(am.sample_id==r.sample_id)&(am.antecedent_hours==r.antecedent_hours)&(am.scope==r.scope)].set_index('case_name')
 for m in ['peak_stage_m','peak_discharge_m3_s','outlet_volume_m3','initial_total_storage_m3','final_total_storage_m3']:
  near('antecedent_paired',g.loc['no_wetland',m]-g.loc['wetland',m],r['delta_'+m],1e-10,m)
 near('antecedent_mass',r.delta_outlet_volume_m3,r.delta_initial_total_storage_m3-r.delta_final_total_storage_m3+r.delta_nonoutlet_external_net_volume_m3+r.delta_end_mass_residual_m3,1e-9,r.scope)
 near('antecedent_rating',r.delta_peak_stage_m,points[r.sample_id]['physical_values']['R_r']*r.delta_peak_discharge_m3_s,4e-13)
for (eid,sid,ah,case),g in am.groupby(['event_id','sample_id','antecedent_hours','case_name']):
 g=g.set_index('scope');full=g.loc['full_horizon'];ev=g.loc['event_window']
 near('window_duration',ev.window_duration_h,144,0)
 near('window_duration',full.window_duration_h,144+ah,0)
 if ah:
  pre=g.loc['antecedent'];near('window_accounting',full.outlet_volume_m3,pre.outlet_volume_m3+ev.outlet_volume_m3,1e-8)
  near('window_accounting',full.peak_stage_m,max(pre.peak_stage_m,ev.peak_stage_m),0)
  near('window_accounting',pre.final_total_storage_m3,ev.initial_total_storage_m3,0)
for _,r in ag.iterrows():
 g=af[(af.event_id==r.event_id)&(af.antecedent_hours==r.antecedent_hours)&(af.scope==r.scope)]
 near('antecedent_summary',g.delta_peak_stage_m.mean()*1000,r.mean_peak_reduction_mm,1e-10)
 near('antecedent_summary',g.delta_outlet_volume_m3.mean(),r.mean_outlet_difference_m3,1e-9)

ledger=[
 ['weather_pilot','primary',20,16,4,16,0,20,'20 generated IDs; four cold preflight exclusions; 16 complete warm pairs'],
 ['factorial','primary',27,27,0,27,3,0,'No probability weights; two pilot corners and prior deterministic center reused'],
 ['fixed_storm_MC','primary',64,64,0,64,2,64,'64 iid exploratory parameter vectors; first two pilot results reused'],
 ['historical_storms','primary',320,320,0,256,76,0,'64 benchmark pairs duplicate MC; 12 new-event primary outcomes copied from pilot; no redraw'],
 ['robustness','reanalysis',320,320,0,0,320,0,'87 zero peak contrasts; 320 endpoint pairs and 258 retained waveform pairs'],
 ['antecedent','primary',12,12,0,12,0,0,'Two events x two prefixes x three saved vectors; six no-prefix baseline pairs reused separately'],
 ['importance_verification','reanalysis',320,320,0,0,320,0,'Existing 64 draws reweighted across five events; no new 32- or 64-draw ensemble'],
 ['weather_refinement','numerical_repeat',16,16,0,0,0,0,'Same 16 warm configurations at 7.5 s; deterministic controls separate'],
 ['parameter_pilot','numerical_repeat',8,8,0,0,0,0,'Four primary configurations x two resolutions; not eight new samples'],
 ['historical_pilot','numerical_repeat',24,24,0,0,0,0,'Twelve primary configurations x two resolutions; 15 s outcomes reused'],
 ['historical_regression','exact_repeat',1,1,0,0,0,0,'Fresh benchmark mc_00000 pair; no new physical configuration'],
 ['antecedent_refinement','numerical_repeat',1,1,0,0,0,0,'One of twelve configurations repeated at 7.5 s; two additional model runs']]
ld=pd.DataFrame(ledger,columns=['study_id','role','planned_pairs_or_evaluations','valid_pairs_or_evaluations','cold_excluded_pairs','additional_primary_configurations_in_reviewed_union','primary_result_copies','new_random_draws_or_weather_ids','interpretation'])
check('ledger',sum(ld.additional_primary_configurations_in_reviewed_union)==375,'scoped unique union')
changed=[str(p.relative_to(ROOT)) for p,h in inputs.items() if digest(p)!=h];check('preservation',not changed,'inputs unchanged during audit')
manifest='\n'.join(f'{p.relative_to(ROOT)}\t{h}' for p,h in sorted(inputs.items(),key=lambda x:str(x[0])))
result={'status':'pass' if not failures else 'fail','audit_utc':datetime.now(timezone.utc).isoformat(),'scope':'Saved JSON/CSV recomputation only. No ODE, fitting, sampling, Git or previous result writes. Dense-array validation is separately recorded in existing independent audits.','checks_total':sum(counts.values()),'checks_by_group':dict(counts),'failures':failures,'new_ODE_runs':0,'new_parameter_draws':0,'unique_reviewed_primary_configurations':375,'distinct_grid_or_MC_parameter_vectors':91,'iid_parameter_draws':64,
 'weather_pilot':{'planned':20,'valid_pairs':16,'cold_excluded':4,'both_initial_peaks':15,'maximum_peak_reduction_m':float(wv.difference_river_peak_m.max())},'parameter_study':parameter,'historical_storms':historical,
 'ecdf':{'rows_recomputed':len(cp),'max_CP95_endpoint_error':max_cp,'scope':'Pointwise fixed-grid intervals under the declared iid parameter priors'},
 'storage':{'endpoint_pairs':len(sp),'full_waveform_pairs':int(sp.full_trajectory_pair_available.sum()),'initial_extra_wetland_water_m3':5500,'maximum_recomputed_identity_error_m3':float(max(abs(recon-sp.delta_outlet_no_minus_wet_m3))),'event_means':se[['event_id','mean_final_stock_term_m3','mean_nonoutlet_external_term_m3','mean_delta_outlet_no_minus_wet_m3','mean_final_groundwater_stock_difference_m3','mean_wetland_compartments_final_m3']].to_dict('records')},
 'importance_reweighting':{'new_draws':0,'uniform_weight_ESS':float(1/sum((raw/raw.sum())**2)),'maximum_uniform_weight':float(max(raw/raw.sum()))},
 'antecedent':{'primary_pairs':12,'baseline_pairs_reused':6,'refinement_pairs':1,'model_executions':26,'event_window_summaries':ag[(ag.scope=='event_window')&(ag.antecedent_hours>0)].to_dict('records'),'minimum_safety_margin_m':aa['minimum_safety_margin_m'],'max_phase_relative_balance':aa['maximum_phase_relative_mass_residual'],'state_seam_error_m3':aa['maximum_state_seam_error_m3']},
 'input_files':len(inputs),'input_collection_sha256':sha256(manifest.encode()).hexdigest(),'changed_inputs':changed,'script_sha256':digest(Path(__file__)),'key_plan_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items() if p.name=='study_plan.json'},
 'limitations':['375 is not a pooled iid sample count or all historical CPU executions.','No accepted downstream flood threshold; initial-reference duration/excess is not flood risk.','64 parameter draws and 27 grid points use explicit hypothetical assumptions, not calibrated empirical distributions.','Five fixed storms are matched repeated evaluations, not 320 independent storms or parameter draws.','Initial wetland model water exceeds no-wetland by 5500 m3; accounting is not causal water tracing.','62 benchmark full waveforms are unavailable; all 64 benchmark endpoint records exist.','Numerical evidence applies to tested points/resolutions only; no universal convergence or real-system certification.','Accepted trajectories exercise energy-gate/closed regimes, not unused transition or frozen-weather physics.','Benchmark and added-event forcing preparation differ as documented.','Importance ESS is not an iid sample count or a binomial interval denominator.','Full-prefix and event-window metrics address different questions; event-start states/controllers are carried, not reset.']}
ld.to_csv(HERE/'study_ledger.csv',index=False)
(HERE/'cross_study_audit.json').write_text(json.dumps(clean(result),indent=2,allow_nan=False)+'\n')
print(json.dumps({'status':result['status'],'checks':result['checks_total'],'failures':failures,'input_files':len(inputs)},indent=2))
raise SystemExit(bool(failures))
