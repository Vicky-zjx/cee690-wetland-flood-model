"""Matched fixed-storm parameter grid, numerical pilot and approved iid MCS."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import get_context
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
import numpy as np
import pandas as pd
from D.d2_analysis.parameter_uncertainty import adapter,preservation,sampling
from D.feedback_control import run_policy as frozen
from D.d2_analysis.analysis import exceedance
CASES=('wetland','no_wetland')
ARRAY_KEYS=('time','x','dx','y','u','balance','downstream_stage','downstream_discharge','state_names','output_names')
METRIC_KEYS=('peak_stage_m','peak_discharge_m3_s','peak_stage_time_s','initial_stage_m','initial_discharge_m3_s',
    'maximum_rise_above_initial_m','initial_level_exceedance_h','initial_level_excess_area_m_h',
    'initial_discharge_exceedance_h','initial_discharge_excess_m3','total_outlet_volume_m3',
    'minimum_safety_margin_m','minimum_geometry_margin_m','balance_relative')


def label(dt):return f'{dt:g}'.replace('.','p')+'s'
def stem(sample,case,dt):return f'{sample}__{case}__{label(dt)}'
def write(path,value):frozen.write_json(path,value)
def history(result):
    return result['controller_dataframe'] if 'controller_dataframe' in result else frozen.controller_table(result.get('controller_history',[]))


def _failure(case,exc,dt):
    return dict(case_name=case,valid=False,controller_history=[],diagnostics=dict(case_name=case,dt_s=float(dt),
        classification='unexpected_execution_error',failure_reason='unexpected_execution_error',
        failure_message=str(exc),exception_type=type(exc).__name__,exception_traceback=traceback.format_exc(),
        full_valid_trajectory_available=False,full_integration_horizon_completed=False,
        ode_started=None,safety_violation=None,
        attempt_note='This case raised before its result object returned; full trajectory unavailable. The other model is still evaluated.'))


def evaluate(point,dt):
    """Separate model catches preserve the second case after any first-case error."""
    results={}
    for case in CASES:
        try:
            runner,records=adapter._private_runner(point['physical_values'],dt)
            r=runner(case=case,scenario_id=None,design='best',dt=float(dt),hydraulics='D',Cw_US=2.8)
            if len(records)!=1:raise RuntimeError('Expected exactly one private parameter object')
            d=r['diagnostics'];d['parameter_evaluation']=records[0]
            d['fixed_configuration']['conditional_on']=adapter.CONDITIONING
            expected=records[0]['initial_wetland_m3' if case=='wetland' else 'initial_no_wetland_m3']
            if d['initial_storage_m3']!=expected:raise RuntimeError('Physical initial state changed')
            results[case]=r
        except Exception as exc:results[case]=_failure(case,exc,dt)
    return results


def _read_saved(directory,name,case):
    d=json.loads((directory/(name+'.json')).read_text())
    with np.load(directory/(name+'.npz'),allow_pickle=False) as z:arrays={k:z[k].copy() for k in ARRAY_KEYS}
    return dict(case_name=case,valid=True,diagnostics=d,**arrays,
        controller_dataframe=pd.read_csv(directory/(name+'_controller.csv.gz'),float_precision='round_trip'),
        reused_from=dict(directory=str(directory.relative_to(ROOT)),name=name,
            npz_sha256=frozen.sha(directory/(name+'.npz')),diagnostics_sha256=frozen.sha(directory/(name+'.json'))))


def reuse(point,dt,plan,phase):
    if phase=='pilot':return None
    if all(point[k]==1 for k in ('m_K_u','m_routing','m_R_r')):
        directory=ROOT/('D/feedback_control/results/fine_nominal' if dt==15 else 'D/d2_analysis/results/fine_7p5')
        results={case:_read_saved(directory,f'deterministic_september_2022__{case}__best__{label(dt)}',case) for case in CASES}
        regression=json.loads((HERE/'results/adapter_regression.json').read_text())
        if not regression['passed']:raise RuntimeError('Nominal adapter regression has not passed')
        for case,r in results.items():r['diagnostics']['parameter_evaluation']=copy.deepcopy(regression['cases'][case]['parameter_evaluation'])
        return results
    for old in plan['pilot']:
        if old['physical_values']==point['physical_values']:
            directory=HERE/'results/pilot'
            names={case:stem(old['sample_id'],case,dt) for case in CASES}
            if all((directory/(name+'.npz')).exists() for name in names.values()):
                return {case:_read_saved(directory,name,case) for case,name in names.items()}
    return None


def save_case(result,point,dt,directory,full,plan):
    case=result['case_name'];name=stem(point['sample_id'],case,dt)
    d=result['diagnostics']
    d['parameter_study']=dict(sample_id=point['sample_id'],study_kind=point['study_kind'],
        multiplier_values={k:point[k] for k in ('m_K_u','m_routing','m_R_r')},physical_values=point['physical_values'],
        study_plan_sha256=frozen.sha(HERE/'results/study_plan.json'),config_sha256=plan['config_sha256'],
        driver_sha256=frozen.sha(__file__),source_sample_id=point.get('source_sample_id'),
        reused_from=result.get('reused_from'),probabilistic_interpretation=point['study_kind']=='explicit_exploratory_mc')
    if 'fixed_configuration' in d:d['fixed_configuration']['conditional_on']=adapter.CONDITIONING
    write(directory/(name+'.json'),d)
    history(result).to_csv(directory/(name+'_controller.csv.gz'),index=False,compression='gzip')
    if full and result['valid']:np.savez_compressed(directory/(name+'.npz'),**{k:result[k] for k in ARRAY_KEYS})


def metrics(result):
    d=result['diagnostics'];obs=d.get('all_rhs_observations',{})
    row=dict(classification=d['classification'],valid=bool(result['valid']),unsafe=d.get('safety_violation'),
        failure_reason=d.get('failure_reason') or d.get('failure_message'),**{k:None for k in METRIC_KEYS})
    row['balance_relative']=d.get('balance_relative')
    safety=obs.get('safety_successful_rhs',[])
    if safety:row['minimum_safety_margin_m']=min(v['margin_m'] for v in safety)
    for item in safety:
        row[item['variable']+'_rhs_peak_m']=item['peak_m'];row[item['variable']+'_rhs_margin_m']=item['margin_m']
    for item in d.get('safety',[]):
        row[item['variable']+'_node_peak_m']=item['peak_m'];row[item['variable']+'_node_margin_m']=item['margin_m']
    row.update(flood_threshold_m=None,flood_duration_h=None,flood_excess_volume_m3=None,
        flood_metrics_status='Unavailable: no accepted downstream flood threshold; initial-reference metrics are separate')
    geometry=obs.get('geometry_attempted_rhs',obs.get('geometry',[]))
    if isinstance(geometry,dict):geometry=list(geometry.values())
    if geometry:row['minimum_geometry_margin_m']=min(v['margin_m'] for v in geometry)
    if not result['valid']:return row
    t=np.asarray(result['time']);h=np.asarray(result['downstream_stage']);q=np.asarray(result['downstream_discharge'])
    duration,area=exceedance(t,h,h[0]);qduration,volume=exceedance(t,q,q[0])
    row.update(peak_stage_m=float(h.max()),peak_discharge_m3_s=float(q.max()),peak_stage_time_s=float(t[np.argmax(h)]),
        initial_stage_m=float(h[0]),initial_discharge_m3_s=float(q[0]),maximum_rise_above_initial_m=float(h.max()-h[0]),
        initial_level_exceedance_h=duration/3600,initial_level_excess_area_m_h=area/3600,
        initial_discharge_exceedance_h=qduration/3600,initial_discharge_excess_m3=volume,
        total_outlet_volume_m3=float(np.trapezoid(q,t)))
    return row


def _history_compare(a,b):
    if len(a)!=len(b) or not np.array_equal(a[['time_s','branch']].to_numpy(),b[['time_s','branch']].to_numpy()):
        return dict(clock_equal=False,latch_mismatch_count=None)
    return dict(clock_equal=True,latch_mismatch_count=int(np.count_nonzero(a.enabled.to_numpy()!=b.enabled.to_numpy())))


def paired_row(point,dt,results):
    p=point['physical_values'];row={k:v for k,v in point.items() if k!='physical_values'}
    row.update(dt_s=float(dt),K_u=p['K_u'],k_ud_1=p['k_ud'][0],k_ud_2=p['k_ud'][1],k_ur=p['k_ur'],R_r=p['R_r'])
    for case,r in results.items():row.update({case+'_'+key:value for key,value in metrics(r).items()})
    valid=all(results[c]['valid'] for c in CASES)
    row.update(paired_valid=False,delta_peak_stage_m=None,delta_peak_discharge_m3_s=None,
        paired_weather_inputs_exact=None,paired_gate_history_agreement=None,maximum_paired_gate_opening_difference_m=None,
        paired_latch_mismatch_count=None,delta_maximum_rise_m=None,delta_initial_level_exceedance_h=None,
        delta_initial_discharge_excess_m3=None,delta_total_outlet_volume_m3=None,failure_reason=None,reused_results=all('reused_from' in results[c] for c in CASES))
    if valid:
        a,b=(results[c] for c in CASES)
        clocks=np.array_equal(a['time'],b['time']);weather=clocks and np.array_equal(a['u'][2:],b['u'][2:])
        common_fields=('effective_physical_values','policy_sha256','prepared_forcing_sha256','boundaries','gate_parameters','design_categories','design_covers')
        am=a['diagnostics']['parameter_evaluation'];bm=b['diagnostics']['parameter_evaluation']
        common=all(am[k]==bm[k] for k in common_fields)
        hc=_history_compare(history(a),history(b))
        row.update(paired_valid=bool(clocks and weather and common),paired_weather_inputs_exact=bool(weather),
            paired_latch_mismatch_count=hc['latch_mismatch_count'])
        if clocks:
            gmax=np.array([g['gate_max_opening'] for g in am['gate_parameters']])[:,None]
            difference=float(np.max(np.abs((a['u'][:2]-b['u'][:2])*gmax)))
            row.update(maximum_paired_gate_opening_difference_m=difference,
                paired_gate_history_agreement=bool(hc['clock_equal'] and hc['latch_mismatch_count']==0 and difference<=1e-12))
        if row['paired_valid']:
            row.update(delta_peak_stage_m=row['no_wetland_peak_stage_m']-row['wetland_peak_stage_m'],
                delta_peak_discharge_m3_s=row['no_wetland_peak_discharge_m3_s']-row['wetland_peak_discharge_m3_s'])
            for out,metric in [('delta_maximum_rise_m','maximum_rise_above_initial_m'),('delta_initial_level_exceedance_h','initial_level_exceedance_h'),('delta_initial_discharge_excess_m3','initial_discharge_excess_m3'),('delta_total_outlet_volume_m3','total_outlet_volume_m3')]:
                row[out]=row['no_wetland_'+metric]-row['wetland_'+metric]

        else:row['failure_reason']='Completed outputs have incompatible pairing provenance, time or weather inputs'
    else:row['failure_reason']='; '.join(c+': '+str(row[c+'_failure_reason'] or row[c+'_classification']) for c in CASES if not results[c]['valid'])
    return row


def convergence(point,coarse,fine,cfg):
    rows=[];module=frozen.historical.adapter();fixed=frozen.historical.load_fixed(dt=15,design='best')
    p,_,_=module._parameters(fixed)
    denominators8=np.r_[p.A_u,p.A_u,p.A_w,p.A_d,p.S_y*p.A_g,p.A_r]
    limits=cfg['pilot']
    for case in CASES:
        a,b=coarse[case],fine[case]
        record=dict(sample_id=point['sample_id'],case_name=case,coarse_dt_s=15,refined_dt_s=7.5,
            both_valid=bool(a['valid'] and b['valid']),coarse_classification=a['diagnostics']['classification'],
            refined_classification=b['diagnostics']['classification'],passed=False)
        if not record['both_valid']:
            record['reason']='At least one trajectory incomplete/invalid; no complete numerical convergence claim'
            rows.append(record);continue
        indices=np.searchsorted(b['time'],a['time']);common=np.array_equal(b['time'][indices],a['time'])
        record['common_grid_exact']=bool(common)
        if not common:record['reason']='No identical common time grid';rows.append(record);continue
        names=list(a['output_names']);hrows=[i for i,name in enumerate(names) if str(name).startswith('H_')]
        den=denominators8 if case=='wetland' else denominators8[module.dry.SHARED]
        state=float(np.max(np.abs(a['x']-b['x'][:,indices])/den[:,None]))
        stages=float(np.max(np.abs(a['y'][hrows]-b['y'][hrows][:,indices])))
        gmax=np.array([g['gate_max_opening'] for g in a['diagnostics']['parameter_evaluation']['gate_parameters']])[:,None]
        gate=float(np.max(np.abs((a['u'][:2]-b['u'][:2,indices])*gmax)))
        peak=abs(float(np.max(a['downstream_stage'])-np.max(b['downstream_stage'])))
        hcheck=_history_compare(history(a),history(b))
        times=history(a).time_s.to_numpy()[::2]
        fixed_clock=bool(hcheck['clock_equal'] and len(times)==4320 and np.array_equal(times,np.arange(4320)*120.))
        mass=bool(all(r['diagnostics']['balance_relative']<cfg['numerics']['mass_balance_relative_tolerance'] for r in (a,b)))
        weather_equal=bool(np.array_equal(a['u'][2:],b['u'][2:,indices]))
        all_peaks=float(np.max(np.abs(np.max(a['y'][hrows],axis=1)-np.max(b['y'][hrows],axis=1))))
        quantitative=bool(weather_equal and all_peaks<=limits['max_physical_peak_stage_difference_m'] and state<=limits['max_state_equivalent_depth_difference_m'] and stages<=limits['max_stage_difference_m'] and
            peak<=limits['max_physical_peak_stage_difference_m'] and gate<=limits['max_actual_gate_opening_difference_m'] and fixed_clock and mass)
        record.update(max_state_equivalent_depth_difference_m=state,max_stage_difference_m=stages,
            full_physical_peak_stage_difference_m=peak,max_all_stage_native_peak_difference_m=all_peaks,
            full_physical_peak_discharge_difference_m3_s=abs(float(np.max(a['downstream_discharge'])-np.max(b['downstream_discharge']))),
            max_actual_gate_opening_difference_m=gate,latch_mismatch_count=hcheck['latch_mismatch_count'],
            controller_clock_120s_exact=fixed_clock,mass_pass=mass,quantitative_pass=quantitative,
            passed=bool(quantitative and hcheck['latch_mismatch_count']==0),
            coarse_unsafe=a['diagnostics'].get('safety_violation'),refined_unsafe=b['diagnostics'].get('safety_violation'),
            weather_inputs_at_common_nodes_exact=bool(np.array_equal(a['u'][2:],b['u'][2:,indices])))
        rows.append(record)
    return rows


def run_point(point,dts,phase,output,plan):
    rows=[];by_dt={};model_rows=[]
    for dt in dts:
        print('START',phase,point['sample_id'],dt,flush=True);started=time.monotonic()
        results=reuse(point,dt,plan,phase) or evaluate(point,dt)
        for result in results.values():
            save_case(result,point,dt,output,phase=='pilot',plan)
            model_rows.append(dict(sample_id=point['sample_id'],study_kind=point['study_kind'],dt_s=dt,
                case_name=result['case_name'],**metrics(result)))
        row=paired_row(point,dt,results);row['elapsed_s']=time.monotonic()-started;rows.append(row)
        by_dt[dt]=results
        print('DONE',point['sample_id'],dt,row['paired_valid'],row['wetland_classification'],row['no_wetland_classification'],round(row['elapsed_s'],2),flush=True)
    checks=convergence(point,by_dt[15.],by_dt[7.5],plan['config']) if phase=='pilot' else []
    return rows,checks,model_rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phase',choices=['prepare','pilot','sensitivity','monte_carlo'],required=True)
    parser.add_argument('--workers',type=int,choices=[1,2],default=1)
    args=parser.parse_args();guard=preservation.verify()
    if not guard['passed']:raise RuntimeError('Protected files changed: '+json.dumps(guard))
    plan=sampling.prepare();cfg=plan['config']
    if args.phase=='prepare':print('Saved immutable 27-grid/64-MC/4-pilot plan; no ODE calls');return
    if args.phase=='monte_carlo':
        summary=json.loads((HERE/'results/pilot/pilot_summary.json').read_text())
        if not summary['passed']:raise RuntimeError('Pilot has unresolved numerical/hysteresis checks; full MCS is blocked pending audit')
    output=HERE/'results'/args.phase
    if output.exists():raise FileExistsError(output)
    output.mkdir(parents=True)
    selected=plan[{'sensitivity':'grid','monte_carlo':'mc'}.get(args.phase,args.phase)]
    dts=[15.,7.5] if args.phase=='pilot' else [15.]
    write(output/'execution_plan.json',dict(study_plan_sha256=frozen.sha(HERE/'results/study_plan.json'),
        driver_sha256=frozen.sha(__file__),points=selected,dt_s=dts,workers=args.workers,
        full_trajectories_saved=args.phase=='pilot',valid_unsafe_retained=True,no_samples_resampled_or_dropped=True))
    rows=[];comparisons=[];model_rows=[];order={p['sample_id']:i for i,p in enumerate(selected)}
    def collect(value):
        new_rows,new_checks,new_models=value
        rows.extend(new_rows);comparisons.extend(new_checks);model_rows.extend(new_models)
        rows.sort(key=lambda r:(order[r['sample_id']],-r['dt_s']))
        model_rows.sort(key=lambda r:(order[r['sample_id']],-r['dt_s'],r['case_name']))
        pd.DataFrame(rows).to_csv(output/'paired_results.csv',index=False)
        pd.DataFrame(model_rows).to_csv(output/'model_results.csv',index=False)
        if comparisons:pd.DataFrame(comparisons).to_csv(output/'numerical_comparison.csv',index=False)
    if args.workers==1:
        for point in selected:collect(run_point(point,dts,args.phase,output,plan))
    else:
        # Independent processes own their parameter/configuration objects and
        # write distinct sample files. Only this parent writes combined tables.
        with ProcessPoolExecutor(max_workers=2,mp_context=get_context('spawn')) as pool:
            pending=[pool.submit(run_point,p,dts,args.phase,output,plan) for p in selected]
            for future in as_completed(pending):collect(future.result())
    if args.phase=='pilot':
        write(output/'pilot_summary.json',dict(passed=all(r['passed'] for r in comparisons),
            planned_points=4,paired_evaluations=8,case_evaluations=16,valid_pairs=sum(r['paired_valid'] for r in rows),
            checks=comparisons,tolerances=cfg['pilot'],
            interpretation='Numerical pilot only; unsafe valid results are retained. Latch mismatches require explicit audit; no full64 started by this command.'))
    write(output/'counts.json',dict(planned_points=len(selected),paired_rows=len(rows),valid_pairs=sum(r['paired_valid'] for r in rows),
        classifications={c:pd.Series([r[c+'_classification'] for r in rows]).value_counts().to_dict() for c in CASES},
        valid_unsafe={c:sum(bool(r[c+'_valid']) and r[c+'_unsafe'] is True for r in rows) for c in CASES}))
    guard=preservation.verify();write(output/'preservation_check.json',guard)
    if not guard['passed']:raise RuntimeError('Protected files changed during execution')


if __name__=='__main__':main()
