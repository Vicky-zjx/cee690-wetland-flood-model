"""Sampled feedback with exact linear actuator motion and official ode4u.

Only the accepted 8/6 water states enter the ODE. Two gate positions and two
latches are controller memory, updated once per accepted 120 s boundary. No
RHS/RK trial state may update the controller. No hydraulic law is changed here.
"""
from __future__ import annotations
import argparse
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import time
import traceback

sys.dont_write_bytecode=True
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/cee690-feedback-mpl')
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
import numpy as np
import pandas as pd
from D.g5_diagnostics import runner as historical
from D.gate_hydraulics_revision import run_study as revision
from D.gate_hydraulics_revision import model as revised_model
from D.feedback_control.policy import Controller, CONTROL_PERIOD_S, POLICY_PATH
from D.feedback_control import preservation


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write_json(path,value):Path(path).write_text(json.dumps(historical.json_clean(value),indent=2,allow_nan=False)+'\n')


@lru_cache(maxsize=1)
def original_weather():
    return revision.weather_manifest()


def _join(chunks):
    return {key:np.concatenate([part[key][...,:-1] for part in chunks[:-1]]+[chunks[-1][key]],axis=-1)
            for key in ('time','x','dx','y','u','balance')}


def run_case(case,scenario_id=None,design='best',dt=120,hydraulics='D',Cw_US=2.8):
    """Execute one frozen event under the same declared sampled policy.

    B is the unchanged original no-dam-spillway law; D is the prior opt-in
    hydraulic variant. Coefficients are not adjusted according to run outcomes.
    Controller memory planned for a failed block is never reported as attained.
    """
    if case not in ('wetland','no_wetland') or hydraulics not in ('B','D'):
        raise ValueError('Choose wetland/no_wetland and B/D hydraulics')
    if dt not in (120,60,30,15) or CONTROL_PERIOD_S!=120 or CONTROL_PERIOD_S%dt:
        raise ValueError('Keep a 120 s controller clock and 120/60/30/15 s ODE resolution')
    if hydraulics=='B' and Cw_US!=2.8:
        raise ValueError('Cw sensitivity is defined only for D hydraulics')
    module=historical.adapter()
    fixed=historical.load_fixed(dt=dt,gate_scale=1,design=design,constant_opening=None)
    params,initial8,gates=module._parameters(fixed)
    initial=initial8.copy() if case=='wetland' else initial8[module.dry.SHARED].copy()
    states=list(module.a2.STATE_NAMES if case=='wetland' else module.dry.STATE_NAMES)
    names=list(module.a2.OUTPUT_NAMES if case=='wetland' else module.dry.OUTPUT_NAMES)
    Y={name:i for i,name in enumerate(names)}
    X={name:i for i,name in enumerate(states)}
    max_open=np.array([g.gate_max_opening for g in gates])
    variant=json.loads((ROOT/'D/gate_hydraulics_revision/variant.json').read_text())
    variant['Cw_US']=float(Cw_US)
    policy=Controller.from_initial(actual_opening_m=[.0015,.0015])
    actual=np.asarray(policy.actual,dtype=float).copy()
    fixed['manifest']['historical_fixed_gate_values']=fixed['manifest'].pop('fixed_gate_values')
    fixed['manifest'].update(schedule_kind='sampled_feedback_exact_linear_actuator',schedule_is_original=False,
        schedule_scale=None,constant_opening=None,controller_policy=policy.config,
        control_period_s=CONTROL_PERIOD_S,initial_actual_opening_m=actual.tolist(),
        controller_state_location='Outside hydrologic ODE: two gate positions and two latches',
        hydraulic_arm=hydraulics,hydraulic_variant_settings=variant if hydraulics=='D' else None,
        conditional_on='Fixed physics, design, initial water, controller policy, actuator assumptions and boundary paths; weather realization varies',
        stage_d_feedback_source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
            (HERE/'run_policy.py',HERE/'policy.py',POLICY_PATH,ROOT/'D/gate_hydraulics_revision/model.py',ROOT/'D/gate_hydraulics_revision/hydraulics.py')})
    sid=scenario_id or 'deterministic_september_2022'
    seed=None;weather_meta=None
    start=fixed['model_config']['geography']['event_start_utc']
    diagnostics=dict(scenario_id=sid,case_name=case,design_source_label=design,dt_s=dt,
        hydraulic_arm=hydraulics,hydraulic_variant_settings=variant if hydraulics=='D' else None,
        fixed_configuration=fixed['manifest'],state_names=states,initial_storage_m3=initial.tolist(),
        original_seed=None,control_period_s=CONTROL_PERIOD_S,hydrologic_state_count=len(states),
        controller_memory_is_not_water_storage=True,ode_started=False,completed_until_s=0.,
        rainfall_area_m2=float(params.A_rain),former_wetland_rainfall_area_m2=float(sum(params.A_w)),
        former_wetland_rainfall_routing='wetland storages' if case=='wetland' else 'direct to river',
        safety_is_a_reported_outcome=True,physical_river_peak_m=None,physical_river_peak_discharge_m3_s=None,
        full_valid_trajectory_available=False,full_integration_horizon_completed=False)
    hourly=None;forcing=None
    if scenario_id is not None:
        weather_meta=next((r for r in original_weather() if r['scenario_id']==scenario_id),None)
        if weather_meta is None:raise ValueError('Unknown original scenario ID')
        if sha(weather_meta['forcing_file'])!=weather_meta['forcing_sha256']:raise ValueError('Original weather bytes changed')
        frame,source_diagnostics=module.read_forcing(weather_meta['forcing_file'])
        diagnostics.update(source_diagnostics)
        start=frame['timestamp_utc'].iloc[0].isoformat();seed=weather_meta['seed_sequence_entropy']
        diagnostics['original_seed']=seed
        diagnostics['event_start_utc']=start
        diagnostics['forcing_conservation']=module.forcing_conservation(frame,fixed)
        if not diagnostics['forcing_conservation']['pass_conservation']:raise ValueError('Weather mapping does not conserve rain/solar input')
        cold=frame.temperature_C<=0
        if cold.any():
            index=int(np.flatnonzero(cold.to_numpy())[0]);failure_time=index*3600.
            diagnostics.update(classification='meteorological_liquid_water_applicability',
                failure_reason='nonfreezing_liquid_water_domain',failure_stage='preflight_model_applicability',
                failure_message='Unchanged liquid-water domain requires T > 0 C; no ODE or controller update occurred',
                first_violation=dict(time_s=failure_time,timestamp_utc=frame.timestamp_utc.iloc[index].isoformat(),
                    temperature_C=float(frame.temperature_C.iloc[index]),time_s_meaning='First excluded hourly interval; no ODE started'),
                first_failure_time_s=failure_time,first_failure_timestamp_utc=frame.timestamp_utc.iloc[index].isoformat(),
                nonpositive_temperature_hours=int(cold.sum()),controller_updates=0,safety_violation=None)
            result=module._empty_result(case,'applicability_excluded',diagnostics,states)
            result['controller_history']=[]
            return result
        hourly=list(module.hourly_input_segments(frame,fixed))
    else:
        forcing=module.wet.fixed_forcing(fixed['model_config'],params)
        diagnostics.update(forcing_mapping='Original deterministic weather knots; no resampling of source data',
            forcing_conservation=dict(source_precipitation_m=float(np.trapezoid(forcing[1][0],forcing[0])),
                source_solar_Wh_m2=float(np.trapezoid(forcing[1][2],forcing[0])/3600)))
    diagnostics['event_start_utc']=start
    hydraulic_diagnostics={}
    raw=(revised_model.make_rhs(case,gates,variant,hydraulic_diagnostics) if hydraulics=='D' else
         module.wet.make_b2_rhs(gates) if case=='wetland' else module.dry.make_rhs(gates))
    observer=historical.RHSObserver(case,fixed)
    solver=observer.wrap_solver(fixed['official']['ode4u'])
    last={}
    gate_observations=[dict(branch=i+1,successful_rhs_evaluations=0,
        minimum_head_m=None,maximum_head_m=None,minimum_actual_opening_m=None,
        maximum_actual_opening_m=None,minimum_open_gate_H_over_a=None,
        minimum_open_gate_ratio_record=None) for i in range(2)]
    def rhs(t,x,u,p):
        last.update(time_s=float(t),state=np.asarray(x).copy(),inputs=np.asarray(u).copy())
        answer=raw(t,x,u,p)
        # Observe successful actual RK calls, without evaluating a second law or
        # changing controller memory. Ratios include every positive opening,
        # including dry/reverse states if such a successful state ever occurs.
        for i,gate in enumerate(gates):
            head=float(x[X[f'S_d{i+1}']]/p.A_d[i]+p.z_d[i]-gate.gate_sill_z)
            opening=float(u[i]*gate.gate_max_opening)
            record=gate_observations[i];record['successful_rhs_evaluations']+=1
            for prefix,value in [('head',head),('actual_opening',opening)]:
                low='minimum_'+prefix+'_m';high='maximum_'+prefix+'_m'
                record[low]=value if record[low] is None else min(record[low],value)
                record[high]=value if record[high] is None else max(record[high],value)
            if opening>0:
                ratio=head/opening
                if record['minimum_open_gate_H_over_a'] is None or ratio<record['minimum_open_gate_H_over_a']:
                    tail=float(x[X[f'S_w{i+1}']]/p.A_w[i]+p.z_w[i] if case=='wetland' else x[X['S_r']]/p.A_r+p.z_r)
                    record['minimum_open_gate_H_over_a']=ratio
                    record['minimum_open_gate_ratio_record']=dict(time_s=float(t),H_over_a=ratio,
                        water_depth_m=head,opening_height_m=opening,tailwater_stage_m=tail,
                        upstream_stage_m=head+gate.gate_sill_z)
        return answer
    chunks=[];history=[];net_total=0.;inflow_total=0.;max_balance=0.;rain_total=0.;solar_total=0.
    current=initial.copy();duration=float(fixed['manifest']['duration_s']);period=CONTROL_PERIOD_S
    current_record=None;failure=None;failure_stage='controller_or_input_preparation';completed=0.
    try:
        for block in range(int(round(duration/period))):
            t0=block*period;t1=t0+period
            heads=[float(current[X[f'S_d{i+1}']]/params.A_d[i]+params.z_d[i]-gates[i].gate_sill_z) for i in range(2)]
            failure_stage='controller_update'
            command=policy.update(heads_m=heads,time_s=t0)
            a0=np.asarray(command['start_opening_m'],dtype=float);a1=np.asarray(command['end_opening_m'],dtype=float)
            target=np.asarray(command['target_opening_m'],dtype=float);rate=np.asarray(command['rate_m_per_s'],dtype=float)
            if any(v.shape!=(2,) or not np.all(np.isfinite(v)) for v in (a0,a1,target,rate)):
                raise ValueError('Controller must return two finite branch values per field')
            if not np.array_equal(a0,actual):raise ValueError('Actuator starting position differs from last attained endpoint')
            if command.get('period_s',period)!=period or np.any(a0<0) or np.any(a1<0) or np.any(a0>max_open) or np.any(a1>max_open):
                raise ValueError('Controller plan violates its fixed period or physical opening bounds')
            if not np.allclose(rate,(a1-a0)/period,atol=1e-15,rtol=0):raise ValueError('Actuator rate does not match planned linear travel')
            if np.any(rate>policy.config['opening_rate_limit_m_per_s']+1e-15) or np.any(-rate>policy.config['closing_rate_limit_m_per_s']+1e-15):
                raise ValueError('Actual actuator plan exceeds declared speed')
            current_record=dict(command,hydraulic_arm=hydraulics,completed_block=False,
                attained_end_opening_m=None,block_end_s=t1,gate_max_opening_m=max_open.tolist())
            history.append(current_record)
            grid=t0+np.arange(int(round(period/dt))+1)*dt
            if hourly is not None:
                hour=int(t0//3600);offset=int(round((t0-hour*3600)/dt));n=len(grid)
                source_grid,source_u=hourly[hour]
                if not np.array_equal(grid,source_grid[offset:offset+n]):raise ValueError('Control and weather clocks do not align exactly')
                weather=source_u[2:,offset:offset+n].copy()
            else:
                weather=np.vstack([np.interp(grid,forcing[0],row) for row in forcing[1]])
            actuator=np.vstack([np.interp(grid,[t0,t1],[a0[i],a1[i]]) for i in range(2)])
            inputs=np.vstack([actuator/max_open[:,None],weather])
            if block==0:
                # Reuse unchanged physical parameter checks with closed gates;
                # the actual selected hydraulic law validates actual g on RHS.
                failure_stage='preflight_model_validation';original_duration=params.duration
                params.duration=period;structural=inputs.copy();structural[:2]=0.
                try:
                    if case=='wetland':module.wet.validate(params,structural,gates)
                    else:module.dry.validate_six(params,structural,current,gates)
                finally:params.duration=original_duration
            failure_stage='ode_rhs_or_output_validation';diagnostics['ode_started']=True
            t,x,dx,y=solver(rhs,grid,current,inputs,params)
            if x.shape!=(len(states),len(t)) or not all(np.all(np.isfinite(v)) for v in (x,dx,y)):
                raise FloatingPointError('Official ode4u returned incomplete/nonfinite arrays')
            if np.any(x<0) or np.any(x[X['S_m']]>params.C_m):raise FloatingPointError('Computed physical storage is inadmissible; no clipping')
            net=module.a2.cumulative_trapezoid(y[Y['external_net']],t)
            balance=x.sum(axis=0)-initial.sum()-(net_total+net)
            max_balance=max(max_balance,float(np.max(np.abs(balance))));net_total+=float(net[-1])
            inflow_total+=float(np.trapezoid(y[Y['P_total']]+inputs[module.a2.U['Q_ext']]+y[Y['Q_gx_in']],t))
            rain_total+=float(np.trapezoid(inputs[module.a2.U['P']],t));solar_total+=float(np.trapezoid(inputs[module.a2.U['R_s']],t)/3600)
            chunks.append(dict(time=t,x=x,dx=dx,y=y,u=inputs,balance=balance))
            current=x[:,-1].copy();actual=a1.copy();completed=float(t[-1])
            current_record.update(completed_block=True,attained_end_opening_m=actual.tolist(),
                endpoint_is_plan_until_ODE_block_completes=False)
    except Exception as exc:
        failure=exc
        if current_record is not None and not current_record['completed_block']:
            current_record.update(failure_time_s=last.get('time_s'),failure_message=str(exc))
    diagnostics.update(completed_until_s=completed,controller_updates=len(history),
        controller_completed_blocks=sum(h['completed_block'] for h in history),
        controller_last_attained_opening_m=actual.tolist(),all_rhs_observations=observer.summary(),
        hydraulic_diagnostics=hydraulic_diagnostics,last_completed_state_m3=current.tolist(),
        actual_gate_successful_rhs_observations=gate_observations,
        completed_blocks_mass_diagnostics=dict(scope='Accepted completed control blocks only; not a full-horizon balance for an incomplete trajectory',
            through_time_s=completed,external_net_volume_m3=net_total,external_inflow_volume_m3=inflow_total,
            maximum_balance_residual_m3=max_balance if chunks else None,
            final_balance_residual_m3=float(chunks[-1]['balance'][-1]) if chunks else None,
            integrated_precipitation_m=rain_total,integrated_solar_Wh_m2=solar_total),
        expected_output_samples=int(duration/dt)+1,
        original_G5_margin_role='Active unchanged hydraulic applicability for B arm; descriptive historical diagnostic for D arm',
        maximum_planned_opening_rate_m_per_s=max([max(h['rate_m_per_s']) for h in history]+[0.]),
        maximum_planned_closing_rate_m_per_s=max([-min(h['rate_m_per_s']) for h in history]+[0.]))
    if failure is not None:
        cause=failure
        while cause.__cause__ is not None:cause=cause.__cause__
        stack=traceback.extract_tb(cause.__traceback__);loc=stack[-1] if stack else None
        failure_time=last.get('time_s') if failure_stage.startswith('ode_') else (current_record or {}).get('time_s',completed)
        diagnostics.update(failure_stage=failure_stage,exception_type=type(failure).__name__,
            root_exception_type=type(cause).__name__,root_exception_message=str(cause),
            failure_message=str(failure),exception_traceback=''.join(traceback.format_exception(failure)),
            failure_location=None if loc is None else dict(file=loc.filename,line=loc.lineno,function=loc.name),
            first_failure_time_s=failure_time,
            first_failure_timestamp_utc=None if failure_time is None else (pd.Timestamp(start)+pd.to_timedelta(failure_time,unit='s')).isoformat(),
            failure_rhs_state_m3=None if not last else last['state'].tolist(),
            failure_rhs_inputs=None if not last else last['inputs'].tolist(),
            failure_actual_opening_m=None if not last else (last['inputs'][:2]*max_open).tolist(),
            attempted_until_s=last.get('time_s'),safety_violation=None)
        if isinstance(failure,module.wet.HydraulicApplicabilityFailure):
            violation=dict(failure.diagnostic)
            H=violation.get('water_depth_m');opening=violation.get('opening_height_m')
            if H is not None and opening is not None:violation['H_over_a']=float(H/opening) if opening else None
            violation.update(scenario_id=sid,case_name=case)
            diagnostics.update(classification='unsupported_gate_hydraulics',failure_reason='unsupported_gate_regime',first_violation=violation)
            status='applicability_excluded'
        elif isinstance(failure,module.wet.TrajectoryFailure):
            diagnostics.update(classification='other_physical_geometry_applicability',failure_reason='hydrologic_geometry_applicability')
            status='applicability_excluded'
        elif isinstance(failure,FloatingPointError):
            diagnostics.update(classification='numerical_or_state_failure',failure_reason='nonfinite_or_inadmissible_storage')
            status='numerical_failure'
        else:
            diagnostics.update(classification='unexpected_execution_error',failure_reason='unexpected_execution_error')
            status='execution_error'
        result=module._empty_result(case,status,diagnostics,states)
        if hydraulics=='D':revision._revised_failure_reporting(result,start,hydraulic_diagnostics)
        if diagnostics.get('first_violation'):
            diagnostics['first_violation'].update(timestamp_utc=diagnostics['first_failure_timestamp_utc'],
                time_s_meaning='First failing official RHS/RK evaluation, not a continuous-time crossing')
    else:
        arrays=_join(chunks);scale=inflow_total if inflow_total>0 else float(initial.sum())
        relative=max_balance/scale if scale>0 else 0. if max_balance==0 else float('inf')
        safety=module._safety_summary(case,arrays['y'],names,fixed['model_config'],arrays['time'])
        rhs_safety=observer.summary()['safety_successful_rhs']
        unsafe=any(r['exceeded'] for r in safety) or any(r['exceeded'] for r in rhs_safety)
        fc=diagnostics['forcing_conservation']
        fc.update(integrated_precipitation_m=rain_total,integrated_solar_Wh_m2=solar_total,
            precipitation_error_m=rain_total-fc['source_precipitation_m'],solar_error_Wh_m2=solar_total-fc['source_solar_Wh_m2'],
            mapping='Exact 120 s controller blocks contained in original weather intervals; P/Rs jumps are not averaged across hour boundaries')
        fc['pass_conservation']=abs(fc['precipitation_error_m'])<=1e-12 and abs(fc['solar_error_Wh_m2'])<=1e-7
        diagnostics.update(balance_max_m3=max_balance,balance_relative=relative,balance_relative_tolerance=float(params.balance_rtol),
            balance_method='External net flux integrated per control block before concatenating; internal transfers cancel once',
            completed_until_s=duration,full_integration_horizon_completed=True,
            safety=safety,safety_violation=bool(unsafe),safety_violation_at_output_nodes=any(r['exceeded'] for r in safety),
            safety_violation_any_successful_rhs=any(r['exceeded'] for r in rhs_safety),
            min_storage_m3=float(arrays['x'].min()),mapped_input_array_sha256=module._array_hash(arrays['time'],arrays['u']),
            solver_termination='Completed fixed interval under fixed-clock feedback')
        if relative>=params.balance_rtol or not fc['pass_conservation']:
            diagnostics.update(classification='numerical_or_state_failure',failure_reason='excessive_mass_balance_residual' if relative>=params.balance_rtol else 'forcing_conservation_failure',
                failure_stage='post_integration_balance_validation',full_valid_trajectory_available=False)
            result=module._empty_result(case,'numerical_failure',diagnostics,states)
        else:
            diagnostics.update(classification='completed_valid_but_safety_exceeded' if unsafe else 'completed_valid_safety_feasible',
                full_valid_trajectory_available=True,
                physical_river_peak_m=float(arrays['y'][Y['H_r']].max()),
                physical_river_peak_time_s=float(arrays['time'][np.argmax(arrays['y'][Y['H_r']])]),
                physical_river_peak_discharge_m3_s=float(arrays['y'][Y['Q_out']].max()),
                physical_river_peak_discharge_time_s=float(arrays['time'][np.argmax(arrays['y'][Y['Q_out']])]))
            result=dict(case_name=case,**arrays,state_names=states,output_names=names,
                downstream_stage=arrays['y'][Y['H_r']],downstream_discharge=arrays['y'][Y['Q_out']],
                solver_status='completed',valid=True,diagnostics=diagnostics)
    result['controller_history']=history
    return result


def controller_table(history):
    rows=[]
    for h in history:
        for i in range(2):
            rows.append(dict(time_s=h['time_s'],block_end_s=h['block_end_s'],branch=i+1,
                head_at_update_m=h['heads_m'][i],start_opening_m=h['start_opening_m'][i],
                target_opening_m=h['target_opening_m'][i],planned_end_opening_m=h['end_opening_m'][i],
                rate_m_per_s=h['rate_m_per_s'][i],enabled=h['enabled'][i],mode=h['modes'][i],
                latch_changed=h['latch_changed'][i],completed_block=h['completed_block'],
                attained_end_opening_m=None if h['attained_end_opening_m'] is None else h['attained_end_opening_m'][i],
                failure_time_s=h.get('failure_time_s')))
    columns=['time_s','block_end_s','branch','head_at_update_m','start_opening_m','target_opening_m',
             'planned_end_opening_m','rate_m_per_s','enabled','mode','latch_changed','completed_block','attained_end_opening_m','failure_time_s']
    return pd.DataFrame(rows,columns=columns)


def save_run(result,directory,name):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    write_json(directory/(name+'.json'),result['diagnostics'])
    controller_table(result.get('controller_history',[])).to_csv(directory/(name+'_controller.csv.gz'),index=False,compression='gzip')
    if result['valid']:
        np.savez_compressed(directory/(name+'.npz'),**{k:result[k] for k in
            ('time','x','dx','y','u','balance','downstream_stage','downstream_discharge','state_names','output_names')})
    d=result['diagnostics'];o=d.get('all_rhs_observations',{})
    row=dict(scenario_id=d['scenario_id'],case_name=result['case_name'],design=d['design_source_label'],
        dt_s=d['dt_s'],hydraulic_arm=d['hydraulic_arm'],classification=d['classification'],
        valid=result['valid'],ode_started=d['ode_started'],controller_updates=d['controller_updates'],
        failure_reason=d.get('failure_reason'),first_failure_time_s=d.get('first_failure_time_s'),
        first_failure_timestamp_utc=d.get('first_failure_timestamp_utc'),
        physical_river_peak_m=d.get('physical_river_peak_m'),physical_river_peak_discharge_m3_s=d.get('physical_river_peak_discharge_m3_s'),
        safety_violation=d.get('safety_violation'),balance_relative=d.get('balance_relative'),
        derivative_mass_identity_m3_s=o.get('derivative_physical_mass_identity_max_abs_m3_s'),
        forcing_sha256=d.get('forcing_file_sha256'),original_seed=json.dumps(d.get('original_seed')))
    violation=d.get('first_violation') or {}
    row.update({k:violation.get(k) for k in ('branch','water_depth_m','opening_height_m','H_over_a','upstream_stage_m','tailwater_stage_m','legacy_G5_residual_m','applicability_residual','applicability_residual_units')})
    for item in o.get('safety_successful_rhs',[]):row[item['variable']+'_peak_m']=item['peak_m'];row[item['variable']+'_margin_m']=item['margin_m']
    return row


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite',choices=['deterministic','weather','all'],default='all')
    parser.add_argument('--dt',type=int,nargs='+',choices=[120,60,30,15],default=[120])
    parser.add_argument('--hydraulics',choices=['B','D'],default='D')
    parser.add_argument('--Cw-US',type=float,choices=[2.6,2.8,3.0],default=2.8)
    parser.add_argument('--tag',required=True)
    args=parser.parse_args()
    if not re.fullmatch('[A-Za-z0-9_-]+',args.tag):raise ValueError('Use a simple output tag')
    guard=preservation.verify()
    if not guard['passed']:raise RuntimeError('Protected files changed: '+json.dumps(guard))
    output=HERE/'results'/args.tag
    if output.exists():raise FileExistsError(output)
    weather=original_weather();plans=[]
    for dt in args.dt:
        if args.suite in ('deterministic','all'):plans.extend(dict(case=c,design=d,dt=dt,scenario_id=None) for c,d in [('wetland','reference'),('wetland','best'),('no_wetland','best')])
        if args.suite in ('weather','all'):plans.extend(dict(case=c,design='best',dt=dt,scenario_id=w['scenario_id']) for w in weather for c in ('wetland','no_wetland'))
    output.mkdir(parents=True)
    write_json(output/'execution_plan.json',dict(plans=plans,weather=weather,policy=json.loads(POLICY_PATH.read_text()),
        policy_sha256=sha(POLICY_PATH),hydraulics=args.hydraulics,Cw_US=args.Cw_US,control_period_s=CONTROL_PERIOD_S))
    rows=[]
    for plan in plans:
        sid=plan['scenario_id'] or 'deterministic_september_2022';name=f'{sid}__{plan["case"]}__{plan["design"]}__{plan["dt"]}s'
        print('START',name,flush=True);start=time.monotonic()
        try:
            result=run_case(**plan,hydraulics=args.hydraulics,Cw_US=args.Cw_US)
            row=save_run(result,output,name)
        except Exception as exc:
            error=dict(scenario_id=sid,case_name=plan['case'],dt_s=plan['dt'],classification='unexpected_execution_error',
                valid=False,exception_type=type(exc).__name__,failure_message=str(exc),exception_traceback=traceback.format_exc())
            write_json(output/(name+'.error.json'),error);row={k:v for k,v in error.items() if k!='exception_traceback'}
        row['elapsed_s']=time.monotonic()-start;rows.append(row)
        pd.DataFrame(rows).to_csv(output/'scenario_diagnostics.csv',index=False)
        print('DONE',name,row.get('classification'),round(row['elapsed_s'],2),flush=True)
    counts=pd.DataFrame(rows).groupby(['case_name','classification'],dropna=False).size().to_dict()
    write_json(output/'counts.json',{f'{k[0]} / {k[1]}':int(v) for k,v in counts.items()})
    guard=preservation.verify();write_json(output/'preservation_check.json',guard)
    if not guard['passed']:raise RuntimeError('Protected files changed during execution')


if __name__=='__main__':main()
