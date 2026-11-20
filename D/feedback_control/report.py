"""Analyze saved feedback runs only. Never import or execute an ODE runner."""
from __future__ import annotations
import hashlib
from functools import lru_cache
import json
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/cee690-feedback-report-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / 'results'
PLOTS = HERE / 'plots'
DET = 'deterministic_september_2022'
CASE_LABEL = {('wetland', 'reference'): 'Wetland A/A reference',
              ('wetland', 'best'): 'Wetland B/B design',
              ('no_wetland', 'best'): 'No wetland'}


def clean(value):
    if isinstance(value, dict): return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [clean(v) for v in value]
    if isinstance(value, np.ndarray): return clean(value.tolist())
    if isinstance(value, np.generic): return clean(value.item())
    if isinstance(value, float) and not np.isfinite(value): return None
    return value


def dump(path, data):
    path.write_text(json.dumps(clean(data), indent=2, allow_nan=False) + '\n')


@lru_cache(maxsize=None)
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def archived_source(expected):
    for folder in (HERE/'source_snapshots',HERE/'source_archive'):
        if folder.exists():
            for p in folder.rglob('*'):
                if p.is_file() and sha(p)==expected:return p
    return None


def array_sha(*arrays):
    h = hashlib.sha256()
    for a in arrays:
        a = np.ascontiguousarray(a)
        h.update(str(a.shape).encode()); h.update(str(a.dtype).encode()); h.update(a.tobytes())
    return h.hexdigest()


def write_table(name, rows):
    frame = pd.DataFrame(rows)
    frame.to_csv(RESULTS / name, index=False)
    return frame


def identity(r):
    return dict(tag=r['tag'], scenario_id=r['sid'], case_name=r['case'], design=r['design'],
                hydraulic_arm=r['arm'], Cw_US=r['cw'] if r['arm']=='D' else None, dt_s=r['dt'])


@lru_cache(maxsize=8)
def load_history(path):
    base=path.with_suffix('')
    choices=[Path(str(base)+'_controller.csv.gz'),Path(str(base)+'_controller.csv')]
    found=next((p for p in choices if p.exists()),None)
    return pd.read_csv(found) if found else pd.DataFrame()

def history(r):return load_history(r['path'])

@lru_cache(maxsize=8)
def load_arrays(path):
    if not path.exists():raise FileNotFoundError(f'Valid result has no saved trajectory: {path}')
    with np.load(path,allow_pickle=False) as z:return {k:z[k] for k in z.files}

def arrays(r):return load_arrays(r['path'].with_suffix('.npz'))


def out(z, name):
    names = list(z['output_names'])
    return z['y'][names.index(name)]


def load_runs():
    runs=[]; audits=[]
    plans=sorted(RESULTS.glob('*/execution_plan.json'))
    if not plans: raise RuntimeError('No saved execution plans; report does not run simulations')
    for planfile in plans:
        directory=planfile.parent; plan=json.loads(planfile.read_text())
        if not (directory/'counts.json').exists():
            raise RuntimeError(f'Suite still running or incomplete: {directory.name}; wait for counts.json')
        for p in plan['plans']:
            sid=p['scenario_id'] or DET
            stem=f'{sid}__{p["case"]}__{p["design"]}__{p["dt"]}s'
            path=directory/(stem+'.json')
            if not path.exists(): path=directory/(stem+'.error.json')
            if not path.exists(): raise RuntimeError(f'Missing planned result {directory.name}/{stem}')
            d=json.loads(path.read_text()); valid=bool(d.get('full_valid_trajectory_available',False))
            r=dict(tag=directory.name,path=path,d=d,sid=sid,case=p['case'],design=p['design'],
                   dt=int(p['dt']),arm=plan['hydraulics'],cw=float(plan['Cw_US']),valid=valid)
            r['classification']=d.get('classification','unclassified_execution_error')
            if valid and not path.with_suffix('.npz').exists():raise FileNotFoundError(f'Valid result has no trajectory: {path}')
            runs.append(r)
        guardpath=directory/'preservation_check.json'
        guard=json.loads(guardpath.read_text()) if guardpath.exists() else {}
        audits.append(dict(tag=directory.name,planned_cases=len(plan['plans']),
                          completed_result_records=len(plan['plans']),policy_sha256=plan.get('policy_sha256'),
                          preservation_passed=guard.get('passed'),execution_plan_sha256=sha(planfile)))
    # A duplicate completed execution must not masquerade as an independent draw.
    keys=[(r['sid'],r['case'],r['design'],r['dt'],r['arm'],r['cw']) for r in runs]
    if len(keys)!=len(set(keys)): raise RuntimeError('Duplicate analysis keys across saved suites; select an unambiguous run set')
    return runs,audits


def safety_items(r):
    d=r['d']; combined={}
    for scope,items in [('output_nodes',d.get('safety',[])),
                        ('successful_RHS',d.get('all_rhs_observations',{}).get('safety_successful_rhs',[]))]:
        for v in items:
            key=v['variable']
            if key not in combined or v['peak_m']>combined[key]['peak_m']:
                combined[key]={**v,'peak_source':scope}
    return list(combined.values())


def min_margin(r, prefix=None):
    a=[v['margin_m'] for v in safety_items(r) if prefix is None or v['variable'].startswith(prefix)]
    return min(a) if a else None


def basic_tables(runs):
    cases=[]; safety=[]; controllers=[]; outcomes=[]; sources=[]
    for r in runs:
        d=r['d']; obs=d.get('all_rhs_observations',{}); violation=d.get('first_violation') or {}
        full=bool(d.get('full_integration_horizon_completed',False))
        fc=d.get('fixed_configuration',{}); sw=safety_items(r)
        row=dict(**identity(r),classification=r['classification'],ode_started=bool(d.get('ode_started',False)),
                 horizon_completed=full,valid_trajectory=r['valid'],
                 hydraulic_horizon_completed=full, numerical_validation_passed=r['valid'],
                 safety_feasible=(not bool(d.get('safety_violation'))) if full else None,
                 any_observed_safety_exceedance=any(v['exceeded'] for v in sw) if sw else None,
                 safety_scope='full event at output and successful RK nodes' if full else 'partial successfully evaluated states only',
                 completed_until_s=d.get('completed_until_s'),failure_stage=d.get('failure_stage'),
                 failure_reason=d.get('failure_reason'),first_failure_time_s=d.get('first_failure_time_s'),
                 first_failure_timestamp_utc=d.get('first_failure_timestamp_utc'),
                 first_failure_branch=violation.get('branch'),failure_H_m=violation.get('water_depth_m'),
                 failure_a_m=violation.get('opening_height_m'),failure_H_over_a=violation.get('H_over_a',violation.get('head_opening_ratio')),
                 failure_message=d.get('failure_message'),exception_type=d.get('root_exception_type',d.get('exception_type')),
                 failure_location=json.dumps(d.get('failure_location')),solver_termination=d.get('solver_termination'),
                 controller_updates=d.get('controller_updates',0),controller_completed_blocks=d.get('controller_completed_blocks',0),
                 balance_max_m3=d.get('balance_max_m3'),balance_relative=d.get('balance_relative'),
                 balance_relative_tolerance=d.get('balance_relative_tolerance'),
                 derivative_mass_identity_m3_s=obs.get('derivative_physical_mass_identity_max_abs_m3_s'),
                 rainfall_mapping_error_m=d.get('forcing_conservation',{}).get('precipitation_error_m'),
                 solar_mapping_error_Wh_m2=d.get('forcing_conservation',{}).get('solar_error_Wh_m2'),
                 forcing_sha256=d.get('forcing_file_sha256'),original_seed=json.dumps(d.get('original_seed')),
                 source_result=str(r['path'].relative_to(ROOT)))
        row['mass_balance_pass']=(row['balance_relative']<row['balance_relative_tolerance']) if row['balance_relative'] is not None and row['balance_relative_tolerance'] is not None else None
        cases.append(row)
        for v in sw: safety.append(dict(**identity(r),**v,full_event=full,scope=row['safety_scope']))
        for kind,mapping in [('historical_source',fc.get('source_sha256',{})),
                             ('feedback_source',fc.get('stage_d_feedback_source_sha256',{}))]:
            for path,expected in mapping.items():
                f=ROOT/path
                current_match=f.exists() and sha(f)==expected
                archived=archived_source(expected) if not current_match else None
                sources.append(dict(**identity(r),source_kind=kind,source_path=path,recorded_sha256=expected,
                                    current_sha256=sha(f) if f.exists() else None,
                                    bytes_match=current_match,verified_archive=str(archived.relative_to(ROOT)) if archived else None,
                                    recorded_bytes_available=bool(current_match or archived)))
        h=history(r)
        if not h.empty:
            for branch,g in h.groupby('branch'):
                g=g.sort_values('time_s'); completed=g[g.completed_block.astype(bool)]
                actual=next((a for a in d.get('actual_gate_successful_rhs_observations',[]) if a['branch']==int(branch)),{})
                ratio=actual.get('minimum_open_gate_ratio_record',{})
                changes=g[g.latch_changed.astype(bool)]
                controllers.append(dict(**identity(r),branch=int(branch),updates=len(g),completed_blocks=len(completed),
                    latch_changes=len(changes),latch_change_times_s=json.dumps(changes.time_s.tolist()),
                    maximum_attained_opening_m=float(max(completed.start_opening_m.max(),completed.attained_end_opening_m.max())) if len(completed) else None,
                    minimum_attained_opening_m=float(min(completed.start_opening_m.min(),completed.attained_end_opening_m.min())) if len(completed) else None,
                    total_attained_travel_m=float(np.abs(completed.attained_end_opening_m-completed.start_opening_m).sum()),
                    maximum_attained_opening_rate_m_s=float(max(0,completed.rate_m_per_s.max())) if len(completed) else None,
                    maximum_attained_closing_rate_m_s=float(max(0,-completed.rate_m_per_s.min())) if len(completed) else None,
                    fully_closed_completed_intervals=int(((completed.start_opening_m==0)&(completed.attained_end_opening_m==0)).sum()),
                    full_open_target_requests=int((g.target_opening_m==.3).sum()),
                    minimum_head_at_updates_m=float(g.head_at_update_m.min()),
                    maximum_head_at_updates_m=float(g.head_at_update_m.max()),full_event=full,
                    minimum_successful_rhs_head_m=actual.get('minimum_head_m'),
                    minimum_successful_rhs_open_gate_H_over_a=actual.get('minimum_open_gate_H_over_a'),
                    minimum_ratio_time_s=ratio.get('time_s'),
                    minimum_ratio_head_m=ratio.get('water_depth_m'),minimum_ratio_opening_m=ratio.get('opening_height_m'),
                    hydraulic_regime_rhs_counts=json.dumps(d.get('hydraulic_diagnostics',{}).get('regime_rhs_counts',{}))))
        result=dict(**identity(r),valid_trajectory=r['valid'],safety_feasible=row['safety_feasible'],
                    physical_river_peak_m=d.get('physical_river_peak_m') if r['valid'] else None,
                    river_peak_at_initial_time=(d.get('physical_river_peak_time_s')==0.) if r['valid'] else None,
                    river_peak_time_s=d.get('physical_river_peak_time_s') if r['valid'] else None,
                    peak_outlet_m3_s=d.get('physical_river_peak_discharge_m3_s') if r['valid'] else None,
                    outlet_peak_time_s=d.get('physical_river_peak_discharge_time_s') if r['valid'] else None,
                    minimum_safety_margin_m=min_margin(r) if full else None)
        if r['valid']:
            z=arrays(r); t=z['time']
            for variable in ['Q_G1','Q_G2','Q_wr1','Q_wr2','Q_wsp1','Q_wsp2','Q_out']:
                if variable in z['output_names']:
                    q=out(z,variable); result[variable+'_volume_m3']=float(np.trapezoid(q,t))
                    result[variable+'_peak_time_s']=float(t[np.argmax(q)])
            for v in sw: result[v['variable']+'_peak_m']=v['peak_m']; result[v['variable']+'_margin_m']=v['margin_m']
        outcomes.append(result)
    return (write_table('case_results.csv',cases),write_table('safety_results.csv',safety),
            write_table('controller_results.csv',controllers),write_table('physical_outcomes.csv',outcomes),
            write_table('source_checks.csv',sources))


def compare_histories(a,b):
    ha,hb=history(a),history(b); result={}
    if ha.empty or hb.empty: return dict(controller_comparison_available=False)
    keys=['time_s','branch']; joined=ha.merge(hb,on=keys,suffixes=('_a','_b'),how='inner')
    result.update(controller_comparison_available=True,common_controller_records=len(joined),
                  controller_time_records_identical=len(joined)==len(ha)==len(hb),
                  latch_state_mismatched_records=int((joined.enabled_a!=joined.enabled_b).sum()),
                  controller_start_opening_max_abs_difference_m=float(np.max(np.abs(joined.start_opening_m_a-joined.start_opening_m_b))) if len(joined) else None)
    for branch in (1,2):
        ta=ha[(ha.branch==branch)&ha.latch_changed.astype(bool)].time_s.to_numpy()
        tb=hb[(hb.branch==branch)&hb.latch_changed.astype(bool)].time_s.to_numpy()
        result[f'gate{branch}_latch_change_count_a']=len(ta);result[f'gate{branch}_latch_change_count_b']=len(tb)
        result[f'gate{branch}_latch_event_times_a_s']=json.dumps(ta.tolist());result[f'gate{branch}_latch_event_times_b_s']=json.dumps(tb.tolist())
        result[f'gate{branch}_max_order_matched_switch_shift_s']=float(np.max(np.abs(ta-tb))) if len(ta)==len(tb) and len(ta) else 0. if len(ta)==len(tb) else None
        result[f'gate{branch}_switch_matching_note']='Events matched by ordinal order, not an exact crossing estimate' if len(ta)==len(tb) else 'Unequal switch counts: ordinal timing comparison unavailable'
    return result


def paired_table(runs):
    rows=[]
    for w in runs:
        if w['case']!='wetland' or w['design']!='best':continue
        matches=[r for r in runs if r['sid']==w['sid'] and r['case']=='no_wetland' and r['dt']==w['dt'] and r['arm']==w['arm'] and r['cw']==w['cw']]
        n=matches[0] if matches else None
        row=dict(**identity(w),wetland_classification=w['classification'],no_wetland_classification=n['classification'] if n else 'missing',
                 both_complete_valid=bool(n and w['valid'] and n['valid']),paired_complete_valid=False,
                 both_hydraulic_horizons_completed=bool(n and w['d'].get('full_integration_horizon_completed') and n['d'].get('full_integration_horizon_completed')),
                 wetland_full_horizon_safety_feasible=not bool(w['d'].get('safety_violation')) if w['d'].get('full_integration_horizon_completed') else None,
                 no_wetland_full_horizon_safety_feasible=not bool(n['d'].get('safety_violation')) if n and n['d'].get('full_integration_horizon_completed') else None)
        if n:
            row.update(same_source_weather_hash=w['d'].get('forcing_file_sha256')==n['d'].get('forcing_file_sha256'),
                       same_original_seed=w['d'].get('original_seed')==n['d'].get('original_seed'),
                       same_model_source_hashes=w['d'].get('fixed_configuration',{}).get('source_sha256')==n['d'].get('fixed_configuration',{}).get('source_sha256'),
                       same_feedback_source_hashes=w['d'].get('fixed_configuration',{}).get('stage_d_feedback_source_sha256')==n['d'].get('fixed_configuration',{}).get('stage_d_feedback_source_sha256'),
                       same_controller_policy=w['d'].get('fixed_configuration',{}).get('controller_policy')==n['d'].get('fixed_configuration',{}).get('controller_policy'))
        if row['both_complete_valid']:
            zw,zn=arrays(w),arrays(n); row['identical_time_grid']=np.array_equal(zw['time'],zn['time'])
            row['identical_mapped_weather']=row['identical_time_grid'] and np.array_equal(zw['u'][2:],zn['u'][2:])
            shared=list(zn['state_names']); wi=[list(zw['state_names']).index(s) for s in shared]
            row['identical_shared_initial_states']=np.array_equal(zw['x'][wi,0],zn['x'][:,0])
            row['extra_initial_wetland_water_m3']=float(zw['x'][:,0].sum()-zn['x'][:,0].sum())
            pair_checks=('identical_time_grid','identical_mapped_weather','identical_shared_initial_states','same_source_weather_hash','same_original_seed','same_controller_policy','same_model_source_hashes','same_feedback_source_hashes')
            row['paired_complete_valid']=all(row[k] for k in pair_checks)
            row['pairing_failure_checks']=','.join(k for k in pair_checks if not row[k])
            if row['identical_time_grid']:
                diff=np.abs((zw['u'][:2]-zn['u'][:2])*.3)
                row.update(actual_gates_identical=np.array_equal(zw['u'][:2],zn['u'][:2]),actual_gate_max_difference_m=float(diff.max()),
                           wetland_actual_gate_sha256=array_sha(zw['time'],zw['u'][:2]),no_wetland_actual_gate_sha256=array_sha(zn['time'],zn['u'][:2]))
            row.update(wetland_peak_m=w['d']['physical_river_peak_m'],no_wetland_peak_m=n['d']['physical_river_peak_m'],
                       wetland_peak_time_s=w['d']['physical_river_peak_time_s'],no_wetland_peak_time_s=n['d']['physical_river_peak_time_s'],
                       both_peak_at_initial_time=w['d']['physical_river_peak_time_s']==n['d']['physical_river_peak_time_s']==0.,
                       no_wetland_minus_wetland_peak_m=(n['d']['physical_river_peak_m']-w['d']['physical_river_peak_m']) if row['paired_complete_valid'] else None,
                       both_safety_feasible=not(w['d']['safety_violation'] or n['d']['safety_violation']))
            row.update(compare_histories(w,n))
        rows.append(row)
    return write_table('paired_results.csv',rows)


def resolution_table(runs):
    rows=[]
    for coarse in runs:
        if coarse['arm']!='D' or coarse['cw']!=2.8 or coarse['dt'] not in (120,60,30):continue
        fine=next((r for r in runs if r['sid']==coarse['sid'] and r['case']==coarse['case'] and r['design']==coarse['design'] and r['dt']==coarse['dt']//2 and r['arm']=='D' and r['cw']==2.8),None)
        if not fine:continue
        row=dict(**identity(coarse),fine_dt_s=fine['dt'],coarse_classification=coarse['classification'],fine_classification=fine['classification'],
                 comparable_full_trajectories=coarse['valid'] and fine['valid'],
                 failure_time_change_s=(fine['d']['first_failure_time_s']-coarse['d']['first_failure_time_s']) if fine['d'].get('first_failure_time_s') is not None and coarse['d'].get('first_failure_time_s') is not None else None)
        row.update(compare_histories(coarse,fine))
        if row['comparable_full_trajectories']:
            a,b=arrays(coarse),arrays(fine); idx=np.searchsorted(b['time'],a['time'])
            if np.any(idx>=len(b['time'])) or not np.array_equal(a['time'],b['time'][idx]):raise ValueError('Resolution grids are not exactly nested')
            row.update(common_grid_river_max_abs_m=float(np.max(np.abs(out(a,'H_r')-out(b,'H_r')[idx]))),
                       common_grid_outlet_max_abs_m3_s=float(np.max(np.abs(out(a,'Q_out')-out(b,'Q_out')[idx]))),
                       common_grid_state_max_abs_m3=float(np.max(np.abs(a['x']-b['x'][:,idx]))),
                       common_grid_actual_gate_max_abs_m=float(np.max(np.abs(a['u'][:2]-b['u'][:2,idx]))*.3),
                       physical_peak_change_m=fine['d']['physical_river_peak_m']-coarse['d']['physical_river_peak_m'],
                       physical_peak_time_change_s=fine['d']['physical_river_peak_time_s']-coarse['d']['physical_river_peak_time_s'],
                       min_safety_margin_change_m=min_margin(fine)-min_margin(coarse))
            for v in ['H_d1','H_d2','H_w1','H_w2']:
                if v in a['output_names']:row[v+'_common_grid_max_abs_m']=float(np.max(np.abs(out(a,v)-out(b,v)[idx])))
        rows.append(row)
    return write_table('resolution_comparison.csv',rows)


def sensitivity_table(runs):
    rows=[]
    for r in runs:
        if r['arm']!='D':continue
        nominal=next((n for n in runs if n['sid']==r['sid'] and n['case']==r['case'] and n['design']==r['design'] and n['dt']==r['dt'] and n['arm']=='D' and n['cw']==2.8),None)
        row=dict(**identity(r),classification=r['classification'],valid_trajectory=r['valid'],
                 minimum_safety_margin_m=min_margin(r) if r['valid'] else None,
                 physical_river_peak_m=r['d'].get('physical_river_peak_m') if r['valid'] else None,
                 physical_peak_outlet_m3_s=r['d'].get('physical_river_peak_discharge_m3_s') if r['valid'] else None)
        if nominal and r['valid'] and nominal['valid']:
            a,b=arrays(r),arrays(nominal)
            if not np.array_equal(a['time'],b['time']):raise ValueError('Sensitivity grids differ')
            row.update(river_peak_change_from_nominal_m=r['d']['physical_river_peak_m']-nominal['d']['physical_river_peak_m'],
                       common_grid_actual_gate_change_from_nominal_m=float(np.max(np.abs(a['u'][:2]-b['u'][:2]))*.3),
                       weather_inputs_unchanged=np.array_equal(a['u'][2:],b['u'][2:]),
                       shared_policy_unchanged=r['d']['fixed_configuration']['controller_policy']==nominal['d']['fixed_configuration']['controller_policy'])
        rows.append(row)
    return write_table('coefficient_sensitivity.csv',rows)


def attribution(runs):
    rows=[]; plotdata={}
    for case,design in [('wetland','reference'),('wetland','best'),('no_wetland','best')]:
        feedback=next((r for r in runs if r['sid']==DET and r['case']==case and r['design']==design and r['dt']==120 and r['arm']=='D' and r['cw']==2.8),None)
        bfeedback=next((r for r in runs if r['sid']==DET and r['case']==case and r['design']==design and r['dt']==120 and r['arm']=='B'),None)
        fixedpath=ROOT/f'D/gate_hydraulics_revision/results/nominal/{DET}__{case}__{design}__120s.json'
        fixed=json.loads(fixedpath.read_text()); compare=fixed['saved_B_baseline_comparison']; bpath=ROOT/compare['saved_trajectory']
        if sha(bpath)!=compare['saved_sha256']:raise ValueError('Historical B trajectory hash changed')
        with np.load(fixedpath.with_suffix('.npz'),allow_pickle=False) as z: zd={k:z[k] for k in z.files}
        with np.load(bpath,allow_pickle=False) as z: zb={k:z[k] for k in z.files if k!='metadata_json'}
        zb['output_names']=zd['output_names'];zb['state_names']=zd['state_names']
        data={'B_fixed':zb,'D_fixed':zd}
        if feedback and feedback['valid']:data['D_feedback']=arrays(feedback)
        if bfeedback and bfeedback['valid']:data['B_feedback']=arrays(bfeedback)
        plotdata[(case,design)]=data
        row=dict(case_name=case,design=design,dt_s=120,B_fixed_source=compare['saved_trajectory'],
                 B_fixed_sha256=sha(bpath),D_fixed_source=str(fixedpath.relative_to(ROOT)),
                 Bfixed_Dfixed_initial_identical=np.array_equal(zb['x'][:,0],zd['x'][:,0]),
                 Bfixed_Dfixed_time_identical=np.array_equal(zb['time'],zd['time']),
                 Bfixed_Dfixed_all_inputs_identical=np.array_equal(zb['u'],zd['u']),
                 D_fixed_safety_feasible=not bool(fixed['safety_violation']),
                 D_feedback_classification=feedback['classification'] if feedback else 'missing',
                 B_feedback_classification=bfeedback['classification'] if bfeedback else 'missing')
        if feedback and feedback['valid']:
            row['D_feedback_safety_feasible']=not feedback['d']['safety_violation']
            row['Dfixed_Dfeedback_initial_identical']=np.array_equal(zd['x'][:,0],data['D_feedback']['x'][:,0])
            row['Dfixed_Dfeedback_weather_identical']=np.array_equal(zd['u'][2:],data['D_feedback']['u'][2:])
        for label,z in data.items():
            for metric,var in [('river_peak_m','H_r'),('outlet_peak_m3_s','Q_out')]:row[label+'_'+metric]=float(out(z,var).max())
        if 'D_feedback' in data:
            for metric in ['river_peak_m','outlet_peak_m3_s']:
                bh=row['B_fixed_'+metric];dh=row['D_fixed_'+metric];df=row['D_feedback_'+metric]
                row['hydraulic_change_under_fixed_'+metric]=dh-bh
                row['control_change_under_D_'+metric]=df-dh
                row['total_path_change_'+metric]=df-bh
                row['additive_identity_error_'+metric]=(dh-bh)+(df-dh)-(df-bh)
        if 'D_feedback' in data and 'B_feedback' in data:
            row['Bfeedback_Dfeedback_gate_max_difference_m']=float(np.max(np.abs(data['D_feedback']['u'][:2]-data['B_feedback']['u'][:2]))*.3)
            row['Bfeedback_Dfeedback_actual_gates_identical']=np.array_equal(data['D_feedback']['u'][:2],data['B_feedback']['u'][:2])
            row['DminusB_under_common_feedback_peak_m']=row['D_feedback_river_peak_m']-row['B_feedback_river_peak_m']
        row['interpretation']='Signed path differences (new minus prior); arithmetic decomposition, not optimization or unconditional causal attribution. Common feedback policy can produce different actual gates.'
        rows.append(row)
    return write_table('deterministic_attribution.csv',rows),plotdata


def savefig(fig,name):
    fig.savefig(PLOTS/name,dpi=180,bbox_inches='tight');plt.close(fig)


def figures(runs,plotdata,resolution,sensitivity,canonical_dt):
    PLOTS.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.2})
    colors={'B_fixed':'#777777','D_fixed':'#d47b25','D_feedback':'#177e73','B_feedback':'#725b9b'}
    labels={'B_fixed':'B hydraulics + fixed gates','D_fixed':'D hydraulics + fixed gates','D_feedback':'D hydraulics + feedback','B_feedback':'B hydraulics + feedback'}
    data=plotdata[('wetland','reference')]
    fig,ax=plt.subplots(2,2,figsize=(12,7.5),sharex=True)
    for i in (0,1):
        for arm,z in data.items():
            if arm!='D_fixed':ax[0,i].plot(z['time']/3600,z['u'][i]*.3,label=labels[arm],color=colors[arm],lw=1.4,ls='--' if arm=='B_fixed' else '-')
            ax[1,i].plot(z['time']/3600,out(z,f'H_d{i+1}')-[100,99.8][i],label=labels[arm],color=colors[arm],lw=1.5,ls='--' if arm=='B_fixed' else '-')
        ax[0,i].set_title(f'Gate {i+1}: actual opening');ax[0,i].set_ylabel('Opening a (m)');ax[0,i].set_ylim(-.01,.32)
        ax[1,i].axhline(1.2,color='black',ls=':',label='Unchanged safety depth 1.20 m')
        ax[1,i].set_title(f'Dam pool {i+1}');ax[1,i].set_ylabel('Depth above gate sill (m)');ax[1,i].set_xlabel('Hours since 2022-09-09 00:00 UTC')
    ax[0,0].legend(fontsize=8,loc='upper left');ax[1,0].legend(fontsize=8,loc='upper left')
    fig.suptitle('Deterministic A/A reference (120 s): operating policy and upstream response',fontsize=14)
    fig.tight_layout();savefig(fig,'01_gate_operations_and_dam_levels.png')
    fig,ax=plt.subplots(2,3,figsize=(15,7.5),sharex=True)
    for j,(key,data) in enumerate(plotdata.items()):
        for arm in ('B_fixed','D_fixed','D_feedback'):
            if arm not in data:continue
            z=data[arm]
            for row,var in [(0,'H_r'),(1,'Q_out')]:ax[row,j].plot(z['time']/3600,out(z,var),label=labels[arm],color=colors[arm],lw=1.4,ls='--' if arm=='B_fixed' else '-')
        ax[0,j].set_title(CASE_LABEL[key]);ax[0,j].set_ylabel('River elevation (m, model datum)');ax[1,j].set_ylabel('Outlet discharge (m³/s)');ax[1,j].set_xlabel('Event time (h)')
    ax[0,0].legend(fontsize=8);fig.suptitle('Deterministic response (120 s): change hydraulics, then change control',fontsize=14)
    fig.tight_layout();savefig(fig,'02_deterministic_river_comparison.png')
    selected=[r for r in runs if r['arm']=='D' and r['cw']==2.8 and r['dt']==canonical_dt and r['sid']!=DET]
    ids=sorted({r['sid'] for r in selected if r['d'].get('ode_started')})
    fig,ax=plt.subplots(2,1,figsize=(13,7.5),sharex=True)
    for case,color in [('wetland','#177e73'),('no_wetland','#a64b36')]:
        records=[next((r for r in selected if r['sid']==sid and r['case']==case),None) for sid in ids]
        peaks=[r['d']['physical_river_peak_m'] if r and r['valid'] else np.nan for r in records]
        margins=[min_margin(r) if r and r['valid'] else np.nan for r in records]
        ax[0].plot(range(len(ids)),peaks,'o-',color=color,label=case.replace('_',' '),lw=1,ms=5)
        ax[1].plot(range(len(ids)),margins,'o-',color=color,label=case.replace('_',' '),lw=1,ms=5)
        for i,r in enumerate(records):
            if r and not r['valid']:ax[0].annotate('invalid',(i,.05),xycoords=('data','axes fraction'),rotation=90,fontsize=8,color=color)
    ax[0].set_ylabel('Physical river peak (m)');ax[0].legend();ax[1].set_ylabel('Minimum dam/wetland safety margin (m)');ax[1].axhline(0,color='black',ls=':')
    ax[1].set_xticks(range(len(ids)),[s.replace('copula_var_','CV ').replace('block_bootstrap_','BB ') for s in ids],rotation=50,ha='right')
    ax[1].set_xlabel('Original scenario ID; four cold IDs remain in the accounting tables')
    fig.suptitle(f'Nominal feedback pilot ({canonical_dt} s): physical peaks and separate safety assessment',fontsize=14)
    fig.tight_layout();savefig(fig,'03_weather_peaks_and_safety.png')
    fig,ax=plt.subplots(2,1,figsize=(13,7.5),sharex=True)
    for cw,color in [(2.6,'#516bb5'),(3.,'#bf7441')]:
        d=sensitivity[(sensitivity.case_name=='wetland')&(sensitivity.Cw_US==cw)&(sensitivity.scenario_id!=DET)].sort_values('dt_s')
        d=d.dropna(subset=['river_peak_change_from_nominal_m']).drop_duplicates('scenario_id').set_index('scenario_id')
        vals=[d.loc[s,'river_peak_change_from_nominal_m']*1000 if s in d.index and 'river_peak_change_from_nominal_m' in d else np.nan for s in ids]
        ax[0].plot(range(len(ids)),vals,'o-',lw=1,color=color,label=f'Cw customary {cw:g} minus nominal 2.8')
    ax[0].axhline(0,color='black',ls=':');ax[0].set_ylabel('Wetland river-peak change (mm)');ax[0].legend(fontsize=9)
    for dt,color in [(120,'#516bb5'),(60,'#177e73'),(30,'#ab6c28')]:
        d=resolution[(resolution.case_name=='wetland')&(resolution.dt_s==dt)&(resolution.scenario_id!=DET)].set_index('scenario_id')
        vals=[d.loc[s,'common_grid_river_max_abs_m']*1000 if s in d.index and 'common_grid_river_max_abs_m' in d else np.nan for s in ids]
        ax[1].plot(range(len(ids)),vals,'o-',lw=1,color=color,label=f'{dt} s vs {dt//2} s')
    ax[1].set_ylabel('Max common-grid river difference (mm)');ax[1].legend(fontsize=9)
    ax[1].set_xticks(range(len(ids)),[s.replace('copula_var_','CV ').replace('block_bootstrap_','BB ') for s in ids],rotation=50,ha='right')
    fig.suptitle('Coefficient and resolution checks under the same feedback policy',fontsize=14)
    fig.tight_layout();savefig(fig,'04_sensitivity_and_resolution.png')


def frame_markdown(frame,cols,limit=None):
    # Avoid an optional tabulate dependency.
    selected=frame.loc[:,[c for c in cols if c in frame.columns]]
    if limit:selected=selected.head(limit)
    def fmt(v):
        if pd.isna(v):return '—'
        if isinstance(v,(float,np.floating)):return f'{v:.8g}'
        return str(v).replace('|','/').replace('\n',' ')
    return '\n'.join(['| '+' | '.join(selected.columns)+' |','| '+' | '.join(['---']*len(selected.columns))+' |']+['| '+' | '.join(fmt(v) for v in row)+' |' for row in selected.itertuples(index=False,name=None)])


def main():
    runs,audits=load_runs()
    cases,safety,controller,outcomes,sources=basic_tables(runs)
    pairs=paired_table(runs);resolution=resolution_table(runs);sensitivity=sensitivity_table(runs)
    attr,plotdata=attribution(runs)
    counts=cases.groupby(['hydraulic_arm','Cw_US','dt_s','case_name','classification'],dropna=False).size().reset_index(name='case_count')
    write_table('outcome_counts.csv',counts.to_dict('records'))
    write_table('suite_checks.csv',audits)
    # Event/weather accounting remains separate: deterministic reference is not a weather draw.
    weather_pairs=pairs[pairs.scenario_id!=DET]
    nominal=weather_pairs[(weather_pairs.hydraulic_arm=='D')&(weather_pairs.Cw_US==2.8)]
    nominal_counts=[]
    for dt,g in nominal.groupby('dt_s'):
        nominal_counts.append(dict(dt_s=int(dt),original_ids=len(g),hydraulically_completed_pairs=int(g.both_hydraulic_horizons_completed.sum()),complete_valid_pairs=int(g.paired_complete_valid.sum()),
            safe_complete_pairs=int((g.paired_complete_valid & g.get('both_safety_feasible',pd.Series(False,index=g.index)).eq(True)).sum()),
            cold_excluded_ids=int((g.wetland_classification=='meteorological_liquid_water_applicability').sum()),
            invalid_or_incomplete_non_cold_ids=int((~g.paired_complete_valid & (g.wetland_classification!='meteorological_liquid_water_applicability')).sum()),
            identical_actual_gate_pairs=int((g.paired_complete_valid & g.get('actual_gates_identical',pd.Series(False,index=g.index)).eq(True)).sum()),
            negative_wetland_peak_reduction_pairs=int((g.get('no_wetland_minus_wetland_peak_m',pd.Series(dtype=float))<0).sum()),
            zero_wetland_peak_reduction_pairs=int((g.get('no_wetland_minus_wetland_peak_m',pd.Series(dtype=float))==0).sum()),
            positive_wetland_peak_reduction_pairs=int((g.get('no_wetland_minus_wetland_peak_m',pd.Series(dtype=float))>0).sum()),
            both_peak_at_initial_time_pairs=int((g.paired_complete_valid & g.get('both_peak_at_initial_time',pd.Series(False,index=g.index)).eq(True)).sum())))
    finite_resolution=resolution[resolution.comparable_full_trajectories]
    common_dt=[int(dt) for dt,g in nominal.groupby('dt_s') if len(g)==20]
    canonical_dt=min(common_dt) if common_dt else None
    if canonical_dt is None: raise RuntimeError('No complete uniform nominal weather suite')
    canonical_pairs=nominal[nominal.dt_s==canonical_dt]
    initial_peak_pairs=int(canonical_pairs.get('both_peak_at_initial_time',pd.Series(False,index=canonical_pairs.index)).eq(True).sum())
    positive_pairs=canonical_pairs[canonical_pairs.no_wetland_minus_wetland_peak_m>0]
    summary=dict(scope='Saved-output analysis only; no model fitting, ODE execution, optimization or weather generation',
        run_records=len(runs),independent_saved_weather_ids=len({r['sid'] for r in runs if r['sid']!=DET}),
        suite_checks=audits,nominal_weather_pair_counts=nominal_counts,
        all_recorded_source_hashes_match=bool(sources.bytes_match.all()),
        all_recorded_source_bytes_available=bool(sources.recorded_bytes_available.all()),
        canonical_uniform_weather_dt_s=canonical_dt,
        canonical_initial_peak_pairs=initial_peak_pairs,
        canonical_positive_peak_contrasts=positive_pairs[['scenario_id','wetland_peak_m','no_wetland_peak_m','wetland_peak_time_s','no_wetland_peak_time_s','no_wetland_minus_wetland_peak_m']].to_dict('records'),
        maximum_completed_balance_relative=cases.balance_relative.max(),
        maximum_completed_balance_m3=cases.balance_max_m3.max(),
        maximum_accepted_balance_relative=cases.loc[cases.valid_trajectory,'balance_relative'].max(),
        canonical_nominal_maximum_balance_relative=cases.loc[(cases.dt_s==canonical_dt)&(cases.hydraulic_arm=='D')&(cases.Cw_US==2.8),'balance_relative'].max(),
        completed_numerical_validation_failure_count=int((cases.horizon_completed&~cases.valid_trajectory).sum()),
        minimum_successful_rhs_open_gate_H_over_a=controller.minimum_successful_rhs_open_gate_H_over_a.min(),
        minimum_successful_rhs_head_m=controller.minimum_successful_rhs_head_m.min(),
        canonical_nominal_minimum_successful_rhs_open_gate_H_over_a=controller.loc[(controller.dt_s==canonical_dt)&(controller.hydraulic_arm=='D')&(controller.Cw_US==2.8),'minimum_successful_rhs_open_gate_H_over_a'].min(),
        maximum_derivative_mass_identity_m3_s=cases.derivative_mass_identity_m3_s.max(),
        maximum_60_to_30_common_grid_river_difference_m=finite_resolution[finite_resolution.dt_s==60].get('common_grid_river_max_abs_m',pd.Series(dtype=float)).max(),
        maximum_30_to_15_common_grid_river_difference_m=finite_resolution[finite_resolution.dt_s==30].get('common_grid_river_max_abs_m',pd.Series(dtype=float)).max(),
        maximum_60_to_30_actual_gate_difference_m=finite_resolution[finite_resolution.dt_s==60].get('common_grid_actual_gate_max_abs_m',pd.Series(dtype=float)).max(),
        counts=counts.to_dict('records'),deterministic_attribution=attr.to_dict('records'),
        interpretation=['Safety exceedances are retained separately from trajectory validity.',
                        'Negative modeled wetland benefit, if present, is retained.',
                        'Shared feedback policy is not automatically identical actual gate operation.',
                        'Coefficient changes can alter actual feedback gate histories as well as hydraulics.',
                        'Grid comparisons are finite numerical consistency checks, not empirical validation or a convergence proof.',
                        'Four original cold scenarios remain in each 20-ID accounting denominator.',
                        'Wetland initial states retain 5500 m3 more water; no redistribution or spin-up.'])
    figures(runs,plotdata,resolution,sensitivity,canonical_dt)
    dump(RESULTS/'summary.json',summary)
    countframe=pd.DataFrame(nominal_counts)
    det=outcomes[(outcomes.scenario_id==DET)&(outcomes.dt_s==120)&(outcomes.hydraulic_arm=='D')&(outcomes.Cw_US==2.8)]
    text='''# Feedback-control saved-results report

This report analyzes saved trajectories and controller records only. The hydraulic candidate is unchanged; the operating policy and its engineering assumptions are documented in [policy_basis.md](policy_basis.md). Results describe the tested hypothetical system, not a calibrated installation or certified controller.

## Completion and safety accounting

'''+frame_markdown(countframe,['dt_s','original_ids','hydraulically_completed_pairs','complete_valid_pairs','safe_complete_pairs','cold_excluded_ids','invalid_or_incomplete_non_cold_ids'])+'''

All 20 original weather IDs are retained. The four cold scenarios are preflight exclusions under the unchanged liquid-water assumption; they have no fabricated trajectories or river peaks. The 16 warm scenarios constitute a conditional non-freezing pilot. Repeated timesteps, hydraulic arms and coefficients are repeated evaluations, not additional independent weather draws. A complete unsafe trajectory remains a valid simulation with an unfavorable safety outcome. Full ODE/hydraulic completion is separate from passing the post-integration quadrature balance check; a complete trajectory can fail that numerical validation without a hydraulic failure. The original failed checks remain recorded, and refinement does not change the tolerance.

## Deterministic nominal feedback results (120 s ODE grid)

'''+frame_markdown(det,['case_name','design','physical_river_peak_m','river_peak_time_s','peak_outlet_m3_s','minimum_safety_margin_m','safety_feasible'])+'''

![Actual gate operation and dam levels](plots/01_gate_operations_and_dam_levels.png)

*September 2022 A/A reference, 120 s ODE grid. Gate panels show actual positions, including rate-limited travel. Pool depths use the unchanged sill datum; dotted lines are unchanged safety limits. Fixed openings are shared by the historical B and prior D cases. The additional B-feedback arm uses the same controller but may realize different openings.*

![River response comparison](plots/02_deterministic_river_comparison.png)

*The three deterministic physical cases use a common 120 s ODE grid. They use the same event and physical initial states within each case. Gray-to-orange changes the hydraulic law under fixed gates; orange-to-green changes the control under D hydraulics. This is a signed path comparison, not an optimized benefit. Unsafe fixed cases are retained. Feedback increases the downstream peak in these deterministic cases while keeping upstream levels below the safety limits.*

## Attribution and realized operation

The exact arithmetic path is `D feedback − B fixed = (D fixed − B fixed) + (D feedback − D fixed)`, applied to physical river and outlet peaks. [deterministic_attribution.csv](results/deterministic_attribution.csv) records the terms, additive residual, source hash and input checks. Negative signed changes mean a lower peak; they must not be confused with a penalized objective.

[paired_results.csv](results/paired_results.csv) checks scenario IDs, original seeds, weather bytes, mapped weather, initial shared states, time grids, policy and model source hashes before reporting a paired difference. It separately compares actual gates and controller latches. If actual gates differ, the comparison is under a common feedback policy, not identical operation. The wetland case retains its accepted additional 5,500 m³ initial wetland water, so this structural contrast is conditional on those accepted initial states, not a comparison with matched total initial water volume. Negative, zero and positive peak reductions are counted separately below; a negative modeled wetland effect is a normal risk result, not an invalid trajectory.

'''+frame_markdown(countframe,['dt_s','identical_actual_gate_pairs','negative_wetland_peak_reduction_pairs','zero_wetland_peak_reduction_pairs','positive_wetland_peak_reduction_pairs'])+'''

'''+f'In the canonical {canonical_dt} s suite, {initial_peak_pairs} pairs attain their physical river maximum at the shared initial time. Their zero peak contrasts do not establish identical downstream trajectories or absence of a wetland routing effect; the whole-window maximum is dominated by the accepted initial river level. The remaining positive contrasts are listed below.\n\n'+frame_markdown(positive_pairs,['scenario_id','wetland_peak_m','no_wetland_peak_m','wetland_peak_time_s','no_wetland_peak_time_s','no_wetland_minus_wetland_peak_m'])+'''

The B-feedback arm is an additional closed-loop hydraulic-model comparison. A difference from D feedback is not a pure discharge-law effect when feedback changes the realized gate history. No counterfactual wetland improvement is interpreted as optimization.

![Original weather outcomes](plots/03_weather_peaks_and_safety.png)

*Nominal D feedback at the finest uniformly completed reported weather resolution (see summary.json). Every started original weather ID appears; unavailable full-event values remain missing. Safety margin is the smallest dam/wetland limit minus the observed maximum, including successful RK evaluations. A positive margin indicates feasibility at the checked evaluations, not continuous-time certification. The four cold exclusions remain in the accounting table. Lines connect scenario summaries for readability; they are not a time series.*

## Resolution, coefficient and controller checks

![Sensitivity and resolution](plots/04_sensitivity_and_resolution.png)

*Top: change in the wetland physical peak for the declared lower/upper coefficient cases, with the policy held fixed. The finest available common numerically valid comparison is shown per scenario; its timestep is explicit in the sensitivity table, including targeted refinements. Failed coarse checks remain in the tables. Actual gate histories can respond differently. Bottom: maximum river-level difference over common output nodes under ODE refinement; the controller clock remains 120 s. These are finite comparisons, not a formal convergence proof.*

[resolution_comparison.csv](results/resolution_comparison.csv) reports 120→60, 60→30 and available 30→15 differences in full trajectories, physical peaks, actual gates and latch states. Latch-event times are matched by ordinal event number only when counts agree. Unequal counts are reported explicitly. [coefficient_sensitivity.csv](results/coefficient_sensitivity.csv) retains all outcomes for all three declared coefficients; failed cases are not assigned substitute peaks.

## Numerical balance and limits

'''+f"Maximum recorded relative integrated balance residual: {summary['maximum_completed_balance_relative']:.8g}; maximum recorded absolute residual: {summary['maximum_completed_balance_m3']:.8g} m³. Maximum instantaneous physical mass-identity mismatch: {summary['maximum_derivative_mass_identity_m3_s']:.8g} m³/s. Each run retains its own tolerance and conservation status in [case_results.csv](results/case_results.csv). These maxima can belong to different runs. Maximum relative residual among accepted trajectories: {summary['maximum_accepted_balance_relative']:.8g}; maximum in the canonical {canonical_dt} s nominal suite: {summary['canonical_nominal_maximum_balance_relative']:.8g}. The unchanged acceptance tolerance is 1e-4.\n\n"+'''
Safety, completion, mass balance and controller behavior are distinct tables. A stopped hydraulic run has partial-state scope; a completed run that fails the post-integration balance check retains full-horizon safety diagnostics but is excluded from accepted physical-peak comparisons until a documented finer-grid run passes. Gate, wetland-drainage and overflow volumes in [physical_outcomes.csv](results/physical_outcomes.csv) are individual internal transfers; they must not be added together as independent external water inputs. No water level, opening command, discharge or weather value was clipped by this report.

The target-opening rule aims for H/a ≥ 2 while enabled; actual rate-limited travel can produce a slightly smaller ratio. [controller_results.csv](results/controller_results.csv) records the minimum ratio over successful RHS evaluations, including its time, head and opening, and the executed hydraulic regimes. The actual ratio is assessed separately from the target ratio; no opening is clipped to force it to two. Remaining in a supported gate regime avoids invoking the unvalidated low-head branches; it does not validate those branches. No downstream flood threshold has been invented: river elevation, outlet discharge and timing are physical outputs.

Hydraulic geometry, the assumed transition, coefficient transfer to small openings, ideal sensors/actuators and omitted freezing remain limitations. Computational completion does not validate those assumptions or establish safety for untested weather. The main README supplies the project readiness decision.

## Reproduce the report

```bash
cd /Users/zhangjiaxuan/Documents/Codex/Modeling
/opt/miniconda3/bin/python3 -B D/feedback_control/report.py
```

This command refuses unfinished suites and reads existing saved outputs; it does not run simulations. It rewrites only this report, its CSV/JSON summaries and four figures. [suite_checks.csv](results/suite_checks.csv) records plan counts and preservation results; [source_checks.csv](results/source_checks.csv) compares recorded source hashes with current bytes and, for a preserved earlier runner revision, verified source snapshots. A current-file mismatch is not hidden by the snapshot column; both are reported separately.
'''
    (HERE/'README_results.md').write_text(text)
    print(json.dumps(clean({k:summary[k] for k in ('run_records','independent_saved_weather_ids','nominal_weather_pair_counts','all_recorded_source_hashes_match','all_recorded_source_bytes_available','canonical_uniform_weather_dt_s','maximum_60_to_30_common_grid_river_difference_m','maximum_30_to_15_common_grid_river_difference_m')}),indent=2))


if __name__=='__main__':main()
