"""Saved-trajectory D.2 metric/safety/controller analysis; never runs an ODE."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from D.d2_analysis import preservation

OUT = HERE/'results'
DET = 'deterministic_september_2022'
BASE = ROOT/'D/feedback_control/results/fine_nominal'
FINE = OUT/'fine_7p5'
POLICY = json.loads((ROOT/'D/feedback_control/policy.json').read_text())
LIMITS = {'H_d1': 101.2, 'H_d2': 101.0, 'H_w1': 98.9, 'H_w2': 98.8}
METRICS = ['river_peak_m', 'outlet_peak_m3_s', 'maximum_rise_above_initial_m',
           'river_drawup_m', 'outlet_drawup_m3_s', 'initial_level_exceedance_h',
           'initial_discharge_excess_m3', 'total_outlet_volume_m3', 'river_peak_time_h']


def clean(x):
    if isinstance(x, dict): return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)): return [clean(v) for v in x]
    if isinstance(x, np.ndarray): return clean(x.tolist())
    if isinstance(x, np.generic): return clean(x.item())
    if isinstance(x, float) and not np.isfinite(x): return None
    return x


def write_json(path, data):
    path.write_text(json.dumps(clean(data), indent=2, allow_nan=False)+'\n')


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def table(name, rows):
    df = pd.DataFrame(rows)
    df.to_csv(OUT/(name+'.csv'), index=False)
    return df


def exceedance(t, value, threshold):
    """Strict exceedance duration and positive area of a piecewise linear curve."""
    t = np.asarray(t, float); z = np.asarray(value, float)-float(threshold)
    if len(t) != len(z) or len(t) < 2 or not np.all(np.isfinite(z)) or np.any(np.diff(t) <= 0):
        raise ValueError('Finite values on a strictly increasing time grid required')
    left, right, width = z[:-1], z[1:], np.diff(t)
    duration = np.zeros_like(width); area = np.zeros_like(width)
    both = (left > 0) & (right > 0)
    duration[both] = width[both]
    area[both] = width[both]*(left[both]+right[both])/2
    crossing = (left > 0) ^ (right > 0)
    positive = np.maximum(left[crossing], right[crossing])
    active = width[crossing]*positive/np.abs(right[crossing]-left[crossing])
    duration[crossing] = active
    area[crossing] = active*positive/2
    return float(duration.sum()), float(area.sum())


def metric_checks():
    checks = []
    for label, values, expected in [('below', [-2,-1],(0,0)), ('above',[1,3],(2,4)),
                                    ('upcross',[-1,3],(1.5,2.25)),('downcross',[3,-1],(1.5,2.25)),
                                    ('equal',[0,0],(0,0)),('zero_to_above',[0,2],(2,2)),
                                    ('above_to_zero',[2,0],(2,2))]:
        actual = exceedance([0,2],values,0)
        checks.append(dict(check=label, passed=np.allclose(actual,expected,rtol=0,atol=1e-14),actual=actual,expected=expected))
    coarse = exceedance([0,2],[-1,3],0)
    fine = exceedance([0,.5,1,1.5,2],[-1,0,1,2,3],0)
    checks.append(dict(check='Crossing integration invariant to linear subdivision',passed=coarse==fine))
    vals=np.array([5,3,4,2,4.5]); rise=float(vals.max()-vals[0]); draw=float(np.max(vals-np.minimum.accumulate(vals)))
    checks.append(dict(check='Initial maximum retained and distinguished from drawup',passed=rise==0 and draw==2.5))
    write_json(OUT/'metric_checks.json',dict(passed=all(c['passed'] for c in checks),checks=checks))
    if not all(c['passed'] for c in checks):raise AssertionError('Metric integration checks failed')


def load_runs(folder):
    if not (folder/'counts.json').exists():raise RuntimeError(f'Incomplete suite: {folder}')
    rows=[]
    for p in sorted(folder.glob('*.json')):
        d=json.loads(p.read_text())
        if not isinstance(d,dict) or 'scenario_id' not in d or 'case_name' not in d:continue
        if p.name.endswith('.error.json'):raise RuntimeError(f'Unexpected execution error must be investigated: {p}')
        valid=bool(d.get('full_valid_trajectory_available',False))
        z=None
        if valid:
            with np.load(p.with_suffix('.npz'),allow_pickle=False) as npz:z={k:npz[k] for k in npz.files}
        rows.append(dict(path=p,d=d,z=z,valid=valid,sid=d['scenario_id'],case=d['case_name'],
                         design=d['design_source_label'],dt=float(d['dt_s'])))
    if len(rows)!=43:raise RuntimeError(f'Expected 43 original case records in {folder}, found {len(rows)}')
    return rows


def identity(r):
    return dict(scenario_id=r['sid'],method='deterministic' if r['sid']==DET else r['sid'].rsplit('_',1)[0],
                case_name=r['case'],design=r['design'],dt_s=r['dt'])


def ycol(z,name): return z['y'][list(z['output_names']).index(name)]
def xcol(z,name): return z['x'][list(z['state_names']).index(name)]


def weather_summary(runs):
    manifest=json.loads((BASE/'execution_plan.json').read_text())['weather']; rows=[]
    for w in manifest:
        f=pd.read_csv(w['forcing_file']); rain=f.precipitation_mm_per_h
        if sha(w['forcing_file'])!=w['forcing_sha256']:raise ValueError('Original weather hash changed')
        positive=np.flatnonzero(rain.to_numpy()>0)
        rows.append(dict(scenario_id=w['scenario_id'],method=w['method'],nonfreezing=w['original_liquid_domain_admissible'],
            total_rain_mm=float(rain.sum()),maximum_hourly_rain_mm_h=float(rain.max()),
            maximum_6h_rain_mm=float(rain.rolling(6).sum().max()),maximum_24h_rain_mm=float(rain.rolling(24).sum().max()),
            wet_hours=int((rain>0).sum()),first_wet_interval_start_h=float(positive[0]) if len(positive) else None,
            last_wet_interval_end_h=float(positive[-1]+1) if len(positive) else None,
            minimum_temperature_C=float(f.temperature_C.min()),forcing_sha256=w['forcing_sha256'],
            original_seed=json.dumps(w['seed_sequence_entropy']),rainfall_origin='Saved stochastic weather; hourly interval rates',
            aggregation='Exact whole-hour source accumulation; no interpolated hourly labels'))
    r=next(r for r in runs if r['sid']==DET and r['case']=='wetland' and r['design']=='best' and r['dt']==15)
    z=r['z']; t=z['time']; rate=z['u'][2]*1000
    # Integrate original continuous piecewise-linear deterministic forcing into clock hours.
    acc=np.concatenate(([0.],np.cumsum(np.diff(t)*(rate[:-1]+rate[1:])/2)))
    hour=np.diff(np.interp(np.arange(145)*3600,t,acc)); rr=pd.Series(hour)
    rows.append(dict(scenario_id=DET,method='deterministic',nonfreezing=True,total_rain_mm=float(hour.sum()),
        maximum_hourly_rain_mm_h=float(hour.max()),maximum_6h_rain_mm=float(rr.rolling(6).sum().max()),
        maximum_24h_rain_mm=float(rr.rolling(24).sum().max()),wet_hours=int((hour>0).sum()),
        first_wet_interval_start_h=float(np.flatnonzero(hour>0)[0]),last_wet_interval_end_h=float(np.flatnonzero(hour>0)[-1]+1),
        minimum_temperature_C=float(z['u'][3].min()),rainfall_origin='Original deterministic regional event forcing',
        aggregation='Clock-hour integrals of preserved piecewise-linear event input; not direct hourly observations'))
    return table('weather_summary',rows)


def evaluate_metrics(r):
    d=r['d']; row=dict(**identity(r),valid=r['valid'],classification=d['classification'],
        ode_started=d.get('ode_started',False),completed=d.get('full_integration_horizon_completed',False),
        safety_feasible=(not d['safety_violation']) if d.get('safety_violation') is not None else None,
        balance_relative=d.get('balance_relative'),balance_tolerance=d.get('balance_relative_tolerance'),
        flood_threshold_m=None,flood_duration_h=None,flood_excess_volume_m3=None,
        flood_metric_status='Undefined: no accepted downstream flood threshold; initial-reference metrics are separate',
        source=str(r['path'].relative_to(ROOT)),source_sha256=sha(r['path']))
    gate_records=d.get('actual_gate_successful_rhs_observations',[])
    ratios=[g['minimum_open_gate_H_over_a'] for g in gate_records if g.get('minimum_open_gate_H_over_a') is not None]
    row.update(minimum_open_gate_H_over_a=min(ratios) if ratios else None,
        hydraulic_regimes=json.dumps(d.get('hydraulic_diagnostics',{}).get('regime_rhs_counts',{}),sort_keys=True),
        derivative_mass_identity_m3_s=d.get('all_rhs_observations',{}).get('derivative_physical_mass_identity_max_abs_m3_s'),
        solver_termination=d.get('solver_termination'),failure_reason=d.get('failure_reason'))
    if not r['valid']:return row
    z=r['z'];t=z['time'];h=z['downstream_stage'];q=z['downstream_discharge']
    duration,area=exceedance(t,h,h[0]); qduration,volume=exceedance(t,q,q[0])
    draw=h-np.minimum.accumulate(h);qdraw=q-np.minimum.accumulate(q)
    row.update(river_initial_stage_m=float(h[0]),river_initial_storage_m3=float(xcol(z,'S_r')[0]),
        outlet_initial_discharge_m3_s=float(q[0]),river_initial_derivative_m3_s=float(z['dx'][list(z['state_names']).index('S_r'),0]),
        river_initial_inflow_m3_s=float(z['dx'][list(z['state_names']).index('S_r'),0]+q[0]+ycol(z,'E_r')[0]),
        river_peak_m=float(h.max()),outlet_peak_m3_s=float(q.max()),
        river_peak_time_h=float(t[np.argmax(h)]/3600),outlet_peak_time_h=float(t[np.argmax(q)]/3600),
        river_peak_at_initial=bool(np.argmax(h)==0),outlet_peak_at_initial=bool(np.argmax(q)==0),
        maximum_rise_above_initial_m=float(h.max()-h[0]),river_drawup_m=float(draw.max()),
        river_drawup_peak_time_h=float(t[np.argmax(draw)]/3600),outlet_drawup_m3_s=float(qdraw.max()),
        initial_level_exceedance_h=duration/3600,initial_level_excess_area_m_h=area/3600,
        initial_discharge_exceedance_h=qduration/3600,initial_discharge_excess_m3=volume,
        total_outlet_volume_m3=float(np.trapezoid(q,t)),final_river_stage_m=float(h[-1]),
        total_initial_water_m3=float(z['x'][:,0].sum()),total_final_water_m3=float(z['x'][:,-1].sum()),
        stage_outlet_linear_relation_max_error_m3_s=float(np.max(np.abs(q-np.maximum(h-97.15,0)/.45))),
        minimum_physical_storage_m3=float(z['x'].min()))
    return row


def safety_records(r):
    d=r['d']; obs={s['variable']:s for s in d.get('all_rhs_observations',{}).get('safety_successful_rhs',[])}
    rows=[]
    for variable,limit in LIMITS.items():
        applicable=not(r['case']=='no_wetland' and variable.startswith('H_w'))
        row=dict(**identity(r),variable=variable,applicable=applicable,threshold_m=limit,
            completed=d.get('full_integration_horizon_completed',False),valid=r['valid'],
            status='not_applicable_no_wetland' if not applicable else 'not_simulated' if not d.get('ode_started') else 'incomplete')
        if applicable and variable in obs:
            s=obs[variable]
            if s['threshold_m']!=limit:raise AssertionError('Safety threshold changed')
            peak=float(s['peak_m']); row.update(successful_RK_peak_m=peak,successful_RK_peak_time_s=s['peak_time_s'])
            if r['valid']:
                y=ycol(r['z'],variable); accepted=float(y.max());peak=max(peak,accepted)
                row.update(accepted_output_peak_m=accepted,accepted_exceedance_h=exceedance(r['z']['time'],y,limit)[0]/3600)
            row.update(peak_m=peak,margin_m=limit-peak,exceeded=peak>limit,
                status='full_horizon_checked' if row['completed'] else 'partial_horizon_checked')
        rows.append(row)
    return rows


def controller_audit(r):
    history=pd.read_csv(r['path'].with_name(r['path'].stem+'_controller.csv.gz'))
    rows=[]
    for branch in (1,2):
        f=history[history.branch==branch].sort_values('time_s') if len(history) else history
        row=dict(**identity(r),branch=branch,records=len(f),evaluated=bool(len(f)))
        if not len(f):rows.append(row);continue
        times=f.time_s.to_numpy();head=f.head_at_update_m.to_numpy();enabled=f.enabled.to_numpy(bool)
        before=np.concatenate(([POLICY['initial_enabled'][branch-1]],enabled[:-1]))
        expected=np.where(head<=.2,False,np.where(head>=.3,True,before))
        targets=np.where(~enabled,0.,np.where(head<=.3,.1,np.where(head>=.8,.3,.1+.2*(head-.3)/.5)))
        start=f.start_opening_m.to_numpy();end=f.planned_end_opening_m.to_numpy();rate=f.rate_m_per_s.to_numpy()
        endpoint=np.where(np.abs(targets-start)<=.0005*120,targets,start+np.sign(targets-start)*.0005*120)
        continuity=float(np.max(abs(start[1:]-end[:-1]))) if len(f)>1 else 0.
        reconstructed_error=None;head_error=None
        if r['valid']:
            z=r['z'];at=np.searchsorted(z['time'],times)
            actual=z['u'][branch-1]*.3
            reconstruction=np.interp(z['time'],np.r_[times,f.block_end_s.iloc[-1]],np.r_[start,end[-1]])
            reconstructed_error=float(np.max(abs(actual-reconstruction)))
            head_error=float(np.max(abs(ycol(z,f'H_d{branch}')[at]-(100. if branch==1 else 99.8)-head)))
        checks=dict(clock=np.array_equal(times,np.arange(len(f))*120),hysteresis=np.array_equal(enabled,expected),
            mode_changes=np.array_equal(f.latch_changed.to_numpy(bool),enabled!=before),
            target=np.max(abs(f.target_opening_m.to_numpy()-targets))<1e-12,
            endpoint=np.max(abs(end-endpoint))<1e-12,rate=np.max(abs(rate-(end-start)/120))<1e-14,
            limits=bool(np.max(abs(rate))<=.0005+1e-14 and np.min(start)>=0 and np.max(start)<=.3 and np.min(end)>=0 and np.max(end)<=.3),
            continuous=continuity<1e-12,actual_history=reconstructed_error is None or reconstructed_error<1e-12,
            accepted_heads=head_error is None or head_error<1e-12)
        switches=times[enabled!=before]
        row.update(passed=all(checks.values()),failed_checks=','.join(k for k,v in checks.items() if not v),
            maximum_actual_history_error_m=reconstructed_error,maximum_head_record_error_m=head_error,
            maximum_actual_rate_m_s=float(np.max(abs(rate))),continuity_error_m=continuity,
            latch_changes=len(switches),latch_change_times_s=json.dumps(switches.tolist()))
        rows.append(row)
    return rows


def pair_records(runs,metrics):
    index={(r['sid'],r['case'],r['design'],r['dt']):r for r in runs};rows=[]
    mindex={(m['scenario_id'],m['case_name'],m['design'],m['dt_s']):m for m in metrics}
    for w in runs:
        if w['case']!='wetland' or w['design']!='best':continue
        n=index[(w['sid'],'no_wetland','best',w['dt'])]
        row=dict(scenario_id=w['sid'],method=identity(w)['method'],dt_s=w['dt'],paired_valid=False,
            wetland_classification=w['d']['classification'],no_wetland_classification=n['d']['classification'],
            paired_hydraulic_completion=all(r['d'].get('full_integration_horizon_completed',False) for r in (w,n)))
        if w['valid'] and n['valid']:
            zw,zn=w['z'],n['z'];shared=list(zn['state_names']);iw=[list(zw['state_names']).index(s) for s in shared]
            fw,fn=w['d']['fixed_configuration'],n['d']['fixed_configuration']
            checks=dict(scenario_and_seed=w['d']['original_seed']==n['d']['original_seed'],
                forcing_bytes=w['d'].get('forcing_file_sha256')==n['d'].get('forcing_file_sha256'),
                time_grid=np.array_equal(zw['time'],zn['time']),weather=np.array_equal(zw['u'][2:],zn['u'][2:]),
                shared_initial=np.array_equal(zw['x'][iw,0],zn['x'][:,0]),
                policy=fw['controller_policy']==fn['controller_policy'],
                hydraulics=w['d']['hydraulic_variant_settings']==n['d']['hydraulic_variant_settings'],
                source_hashes=fw['stage_d_feedback_source_sha256']==fn['stage_d_feedback_source_sha256'])
            row.update(paired_valid=all(checks.values()),failed_pair_checks=','.join(k for k,v in checks.items() if not v),
                actual_gate_max_abs_difference_m=float(np.max(abs(zw['u'][:2]-zn['u'][:2]))*.3),
                actual_gate_arrays_identical=np.array_equal(zw['u'][:2],zn['u'][:2]),
                extra_initial_wetland_water_m3=float(zw['x'][:,0].sum()-zn['x'][:,0].sum()),
                both_safety_feasible=not(w['d']['safety_violation'] or n['d']['safety_violation']))
            if row['paired_valid']:
                mw=mindex[(w['sid'],w['case'],w['design'],w['dt'])];mn=mindex[(n['sid'],n['case'],n['design'],n['dt'])]
                for metric in METRICS:
                    row['wetland_'+metric]=mw[metric];row['no_wetland_'+metric]=mn[metric]
                    row['difference_'+metric]=mn[metric]-mw[metric]
                hdiff=zn['downstream_stage']-zw['downstream_stage'];qdiff=zn['downstream_discharge']-zw['downstream_discharge'];t=zw['time']
                row.update(paired_stage_difference_integral_m_h=float(np.trapezoid(hdiff,t)/3600),
                    paired_stage_difference_max_m=float(hdiff.max()),paired_stage_difference_min_m=float(hdiff.min()),
                    paired_stage_difference_max_time_h=float(t[np.argmax(hdiff)]/3600),
                    paired_stage_difference_min_time_h=float(t[np.argmin(hdiff)]/3600),
                    paired_discharge_difference_integral_m3=float(np.trapezoid(qdiff,t)),
                    difference_total_final_water_m3=float(zn['x'][:,-1].sum()-zw['x'][:,-1].sum()),
                    wetland_flood_duration_h=None,no_wetland_flood_duration_h=None,difference_flood_duration_h=None,
                    wetland_flood_excess_volume_m3=None,no_wetland_flood_excess_volume_m3=None,difference_flood_excess_volume_m3=None)
        rows.append(row)
    return rows


def convergence(runs,metrics):
    idx={(r['sid'],r['case'],r['design'],r['dt']):r for r in runs};rows=[]
    midx={(r['scenario_id'],r['case_name'],r['design'],r['dt_s']):r for r in metrics}
    for a in runs:
        if a['dt']!=15:continue
        b=idx[(a['sid'],a['case'],a['design'],7.5)]
        row=dict(**identity(a),fine_dt_s=7.5,coarse_classification=a['d']['classification'],
            fine_classification=b['d']['classification'],comparable=a['valid'] and b['valid'])
        if row['comparable']:
            za,zb=a['z'],b['z'];positions=np.searchsorted(zb['time'],za['time'])
            if not np.array_equal(za['time'],zb['time'][positions]):raise ValueError('Nested grids expected')
            row.update(common_weather_identical=np.array_equal(za['u'][2:],zb['u'][2:,positions]),
                initial_states_identical=np.array_equal(za['x'][:,0],zb['x'][:,0]),
                controller_policy_identical=a['d']['fixed_configuration']['controller_policy']==b['d']['fixed_configuration']['controller_policy'],
                hydraulic_settings_identical=a['d']['hydraulic_variant_settings']==b['d']['hydraulic_variant_settings'],
                original_seed_identical=a['d']['original_seed']==b['d']['original_seed'],
                forcing_bytes_identical=a['d'].get('forcing_file_sha256')==b['d'].get('forcing_file_sha256'),
                river_common_grid_max_abs_m=float(np.max(abs(za['downstream_stage']-zb['downstream_stage'][positions]))),
                discharge_common_grid_max_abs_m3_s=float(np.max(abs(za['downstream_discharge']-zb['downstream_discharge'][positions]))),
                state_common_grid_max_abs_m3=float(np.max(abs(za['x']-zb['x'][:,positions]))),
                gate_common_grid_max_abs_m=float(np.max(abs(za['u'][:2]-zb['u'][:2,positions]))*.3),
                coarse_balance_relative=a['d']['balance_relative'],fine_balance_relative=b['d']['balance_relative'],
                coarse_balance_max_m3=a['d']['balance_max_m3'],fine_balance_max_m3=b['d']['balance_max_m3'])
            for name in LIMITS:
                if name in za['output_names']:row[name+'_common_grid_max_abs_m']=float(np.max(abs(ycol(za,name)-ycol(zb,name)[positions])))
            ma=midx[(a['sid'],a['case'],a['design'],15.)];mb=midx[(b['sid'],b['case'],b['design'],7.5)]
            for metric in METRICS:row['fine_minus_coarse_'+metric]=mb[metric]-ma[metric]
            ca=pd.read_csv(a['path'].with_name(a['path'].stem+'_controller.csv.gz'));cb=pd.read_csv(b['path'].with_name(b['path'].stem+'_controller.csv.gz'))
            row['controller_clocks_identical']=np.array_equal(ca[['time_s','branch']],cb[['time_s','branch']])
            mismatch=ca.enabled.to_numpy()!=cb.enabled.to_numpy();row['latch_record_mismatch_count']=int(mismatch.sum())
            row['first_latch_mismatch_s']=float(ca.time_s.iloc[np.flatnonzero(mismatch)[0]]) if mismatch.any() else None
            for branch in (1,2):
                ta=ca.loc[(ca.branch==branch)&ca.latch_changed,'time_s'].to_numpy();tb=cb.loc[(cb.branch==branch)&cb.latch_changed,'time_s'].to_numpy()
                row[f'branch{branch}_switches_coarse']=json.dumps(ta.tolist());row[f'branch{branch}_switches_fine']=json.dumps(tb.tolist())
                row[f'branch{branch}_max_switch_time_change_s']=float(np.max(abs(ta-tb))) if len(ta)==len(tb) and len(ta) else 0. if len(ta)==len(tb) else None
                row[f'branch{branch}_switch_counts_equal']=len(ta)==len(tb)
        rows.append(row)
    return rows


def distribution_tables(metrics,pairs):
    rows=[];ecdf=[]
    for dt in (15.,7.5):
        p=pairs[(pairs.dt_s==dt)&(pairs.scenario_id!=DET)&pairs.paired_valid]
        for population in ['all_warm','copula_var','block_bootstrap']:
            g=p if population=='all_warm' else p[p.method==population]
            for model,prefix in [('wetland','wetland_'),('no_wetland','no_wetland_'),('paired_no_minus_wet','difference_')]:
                for metric in METRICS:
                    values=g[prefix+metric].dropna().to_numpy(float)
                    if not len(values):continue
                    n=len(values);unique,count=np.unique(values,return_counts=True)
                    rows.append(dict(dt_s=dt,population=population,series=model,metric=metric,n=n,
                        mean=float(values.mean()),sample_std=float(values.std(ddof=1)) if n>1 else None,
                        minimum=float(values.min()),q25=float(np.quantile(values,.25)),median=float(np.median(values)),
                        q75=float(np.quantile(values,.75)),q95=float(np.quantile(values,.95)),maximum=float(values.max()),
                        negative_count=int((values<0).sum()),zero_count=int((values==0).sum()),positive_count=int((values>0).sum()),
                        quantile_method='NumPy linear sample interpolation; descriptive, not a return level'))
                    for v,c,cum in zip(unique,count,np.cumsum(count)):
                        ecdf.append(dict(dt_s=dt,population=population,series=model,metric=metric,value=float(v),
                            tied_count=int(c),sample_n=n,ecdf=float(cum/n)))
    table('summary_statistics',rows);table('ecdf_values',ecdf)


def main():
    OUT.mkdir(exist_ok=True);metric_checks()
    config=json.loads((ROOT/'B/no_spillway_config.json').read_text())
    if config['reporting']['flood_threshold_m'] is not None:raise ValueError('Reassess documented threshold before changing metrics')
    runs=load_runs(BASE)+load_runs(FINE)
    weather=weather_summary(runs);metricrows=[evaluate_metrics(r) for r in runs]
    metrics=table('metrics',metricrows)
    safety=table('safety_audit',[row for r in runs for row in safety_records(r)])
    controls=table('controller_audit',[row for r in runs for row in controller_audit(r)])
    pairs=table('paired_metrics',pair_records(runs,metricrows))
    conv=table('convergence',convergence(runs,metricrows));distribution_tables(metrics,pairs)
    accounting=[]
    for dt in (15.,7.5):
        f=metrics[(metrics.dt_s==dt)&(metrics.scenario_id!=DET)];p=pairs[(pairs.dt_s==dt)&(pairs.scenario_id!=DET)]
        accounting.append(dict(dt_s=dt,original_weather_ids=20,model_records=len(f),warm_weather_ids=16,
            cold_excluded_ids=int((f.classification=='meteorological_liquid_water_applicability').sum()/2),
            wetland_started=int(f.loc[f.case_name=='wetland','ode_started'].sum()),no_wetland_started=int(f.loc[f.case_name=='no_wetland','ode_started'].sum()),
            wetland_completed=int(f.loc[f.case_name=='wetland','completed'].sum()),no_wetland_completed=int(f.loc[f.case_name=='no_wetland','completed'].sum()),
            complete_valid_pairs=int(p.paired_valid.sum()),unsafe_valid_model_runs=int((f.valid&f.safety_feasible.eq(False)).sum()),
            incomplete_warm_model_runs=int((f.ode_started&~f.completed).sum()),
            completed_numerically_invalid_runs=int((f.completed&~f.valid).sum())))
    table('accounting',accounting)
    # Preserve previous controlled attribution; no failed fixed-weather peak is manufactured.
    attribution=pd.read_csv(ROOT/'D/feedback_control/results/deterministic_attribution.csv')
    attribution.to_csv(OUT/'policy_hydraulics_attribution.csv',index=False)
    valid=metrics[metrics.valid];warm=metrics[(metrics.dt_s==15)&(metrics.scenario_id!=DET)&metrics.valid]
    c=conv[conv.comparable];sw=safety[(safety.dt_s==7.5)&(safety.scenario_id!=DET)&safety.applicable]
    summary=dict(scope='Preliminary fixed-policy paired analysis of saved original weather; no new weather or full Monte Carlo',
        primary_dt_s=15,fine_dt_s=7.5,control_period_s=120,accounting=accounting,
        flood_threshold_m=None,flood_probability_estimated=False,initial_conditions_changed=False,hydrologic_spinup_performed=False,
        primary_initial_peak_model_runs=int(warm.river_peak_at_initial.sum()),
        primary_warm_method_counts=weather[(weather.scenario_id!=DET)&weather.nonfreezing].groupby('method').size().to_dict(),
        controller_audits_passed=bool(controls.loc[controls.evaluated,'passed'].all()),
        safety_violations_successful_evaluations=int(safety.exceeded.eq(True).sum()),
        minimum_fine_weather_safety_margin_m=float(sw.margin_m.min()),
        max_15_to_7p5_river_difference_m=float(c.river_common_grid_max_abs_m.max()),
        max_15_to_7p5_discharge_difference_m3_s=float(c.discharge_common_grid_max_abs_m3_s.max()),
        max_15_to_7p5_gate_difference_m=float(c.gate_common_grid_max_abs_m.max()),
        total_latch_record_mismatches=int(c.latch_record_mismatch_count.sum()),
        maximum_fine_balance_relative=float(valid.loc[valid.dt_s==7.5,'balance_relative'].max()),
        all_common_weather_identical=bool(c.common_weather_identical.all()),all_initial_states_identical=bool(c.initial_states_identical.all()),
        all_controller_policies_identical=bool(c.controller_policy_identical.all()),all_hydraulic_settings_identical=bool(c.hydraulic_settings_identical.all()),
        minimum_fine_open_gate_H_over_a=float(valid.loc[valid.dt_s==7.5,'minimum_open_gate_H_over_a'].min()),
        maximum_pair_gate_difference_m=float(pairs.actual_gate_max_abs_difference_m.max()),
        plotting_denominator_note='16 warm pairs:10 copula/VAR and6 block-bootstrap;4 cold original IDs remain in accounting; no representative flood probability asserted')
    write_json(OUT/'summary.json',summary)
    input_files=[ROOT/'D/feedback_control/policy.json',ROOT/'B/no_spillway_config.json',BASE/'execution_plan.json',FINE/'execution_plan.json']
    write_json(OUT/'analysis_manifest.json',dict(analysis_sha256=sha(Path(__file__)),
        sources={str(p.relative_to(ROOT)):sha(p) for p in input_files},
        primary_metrics='Full144h maxima includingt0; maximumriseaboveinitial separately; nofloodthreshold',
        reference_excess='Use exactsavedinitialH/Q, notroundeddisplayvalues; benchmarknotflood',
        pair_difference_sign='no_wetland minus wetland; not automatically a causal wetland effect',
        same_policy=True,common_physical_initial_states=True,extra_initial_wetland_water_m3=5500,
        draws_regenerated=False,ode_executed_by_analysis=False))
    write_json(OUT/'analysis_status.json',dict(status='Completed saved-output analysis at15s and7.5s',final=True))
    check=preservation.verify();write_json(HERE/'preservation_check.json',check)
    if not check['passed']:raise RuntimeError('Protected files changed')
    print(json.dumps(clean(summary),indent=2))


if __name__=='__main__':main()
