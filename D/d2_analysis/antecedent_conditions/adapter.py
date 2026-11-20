"""Opt-in continuous antecedent-plus-event integration with observational audits.

The frozen feedback loop, physical RHS, official solver, controller and original
event input mapping are reused. No state or controller is restarted at the event
boundary. Only private duration/forcing providers and diagnostic observers vary.
"""
from __future__ import annotations
import copy
import json
from pathlib import Path
import sys
from types import FunctionType, SimpleNamespace
sys.dont_write_bytecode = True
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
import numpy as np
import pandas as pd
from D.d2_analysis.parameter_uncertainty import adapter as parameter_adapter
from D.feedback_control import run_policy as frozen
from D.d2_analysis.antecedent_conditions.forcing import provider

EVENT_HOURS=144
CONDITIONING=('Accepted physical initial storages and actuator/latch initialization at the antecedent start, '
    'including 5500 m3 additional wetland water. One continuous unchanged physical/controller run '
    'through complete existing antecedent forcing and the unchanged 144-hour event. Event-start '
    'states and controller memory are evolved conditions, not observations or reset states.')


def _resolve(path):
    path=Path(path)
    return path if path.is_absolute() else ROOT/path


class BlockAudit:
    """Observe completed official-solver blocks; never call the RHS a second time."""
    def __init__(self,case,fixed,module,offset):
        self.case=case;self.fixed=fixed;self.module=module;self.offset=float(offset)
        self.total=float(fixed['manifest']['duration_s'])
        self.names=list(module.a2.OUTPUT_NAMES if case=='wetland' else module.dry.OUTPUT_NAMES)
        self.states=list(module.a2.STATE_NAMES if case=='wetland' else module.dry.STATE_NAMES)
        self.Y={name:i for i,name in enumerate(self.names)};self.X={name:i for i,name in enumerate(self.states)}
        self.gates=module.wet.gate_parameters(fixed['model_config'])
        self.max_open=np.asarray([g.gate_max_opening for g in self.gates])
        self.observers={key:frozen.historical.RHSObserver(case,fixed) for key in ('antecedent','event')}
        self.phases={key:dict(start_time_s=0. if key=='antecedent' else self.offset,
            target_end_time_s=self.offset if key=='antecedent' else self.total,
            attempted_blocks=0,completed_blocks=0,initial_storage_m3=None,final_storage_m3=None,
            completed_until_s=0. if key=='antecedent' else self.offset,
            external_net_volume_m3=0.,outlet_volume_m3=0.,external_inflow_volume_m3=0.,
            integrated_precipitation_m=0.,integrated_solar_Wh_m2=0.,
            maximum_balance_residual_m3=0.,final_balance_residual_m3=None,
            minimum_open_gate_H_over_a=[None,None]) for key in ('antecedent','event')}
        self.last_end=None;self.last_state=None;self.last_gate=None;self.rows=[]
        self.boundary=dict(reached=False,state_exact=None,actual_opening_exact=None,
            boundary_time_s=self.offset,controller_reinitialized=False)
        self.maximum_state_continuity_error_m3=0.;self.maximum_gate_continuity_error_m=0.

    def wrap(self,solver):
        def observed(rhs,time,initial,inputs,params):
            t0=float(time[0]);t1=float(time[-1])
            if t1-t0!=120. or t0%120.:raise ValueError('Antecedent observer requires the unchanged 120 s controller blocks')
            phase='antecedent' if t0<self.offset else 'event'
            if phase=='antecedent' and t1>self.offset:raise ValueError('Control block crosses the event boundary')
            p=self.phases[phase];p['attempted_blocks']+=1
            initial=np.asarray(initial);inputs=np.asarray(inputs)
            start_gate=inputs[:2,0]*self.max_open
            state_error=gate_error=0.
            if self.last_end is not None:
                if t0!=self.last_end:raise ValueError('Gap or duplicated block in continuous physical integration')
                state_error=float(np.max(np.abs(initial-self.last_state)))
                gate_error=float(np.max(np.abs(start_gate-self.last_gate)))
                self.maximum_state_continuity_error_m3=max(self.maximum_state_continuity_error_m3,state_error)
                self.maximum_gate_continuity_error_m=max(self.maximum_gate_continuity_error_m,gate_error)
                if not np.array_equal(initial,self.last_state):raise ValueError('Physical state was reset between control blocks')
                if not np.array_equal(start_gate,self.last_gate):raise ValueError('Actual gate opening was reset between control blocks')
            if p['initial_storage_m3'] is None:p['initial_storage_m3']=initial.tolist()
            if t0==self.offset:
                self.boundary.update(reached=True,state_before_m3=self.last_state.tolist(),state_after_m3=initial.tolist(),
                    actual_opening_before_m=self.last_gate.tolist(),actual_opening_after_m=start_gate.tolist(),
                    state_exact=bool(np.array_equal(initial,self.last_state)),
                    actual_opening_exact=bool(np.array_equal(start_gate,self.last_gate)))
            def ratio_rhs(t,x,u,c):
                answer=rhs(t,x,u,c)
                for i,g in enumerate(self.gates):
                    opening=float(u[i]*g.gate_max_opening)
                    if opening>0:
                        head=float(x[self.X[f'S_d{i+1}']]/c.A_d[i]+c.z_d[i]-g.gate_sill_z)
                        ratio=head/opening;old=p['minimum_open_gate_H_over_a'][i]
                        p['minimum_open_gate_H_over_a'][i]=ratio if old is None else min(old,ratio)
                return answer
            # The phase observer calls ratio_rhs once per ordinary RK evaluation;
            # the unchanged full-run observer remains outside this wrapper.
            t,x,dx,y=self.observers[phase].wrap_solver(solver)(ratio_rhs,time,initial,inputs,params)
            net=self.module.a2.cumulative_trapezoid(y[self.Y['external_net']],t)
            balance=x.sum(axis=0)-sum(p['initial_storage_m3'])-(p['external_net_volume_m3']+net)
            block_net=float(net[-1]);block_out=float(np.trapezoid(y[self.Y['Q_out']],t))
            block_in=float(np.trapezoid(y[self.Y['P_total']]+inputs[self.module.a2.U['Q_ext']]+y[self.Y['Q_gx_in']],t))
            block_rain=float(np.trapezoid(inputs[self.module.a2.U['P']],t))
            block_solar=float(np.trapezoid(inputs[self.module.a2.U['R_s']],t)/3600.)
            p['external_net_volume_m3']+=block_net;p['outlet_volume_m3']+=block_out;p['external_inflow_volume_m3']+=block_in
            p['integrated_precipitation_m']+=block_rain;p['integrated_solar_Wh_m2']+=block_solar
            p['maximum_balance_residual_m3']=max(p['maximum_balance_residual_m3'],float(np.max(np.abs(balance))))
            p['final_balance_residual_m3']=float(balance[-1]);p['final_storage_m3']=x[:,-1].tolist()
            p['completed_blocks']+=1;p['completed_until_s']=t1
            self.rows.append(dict(phase=phase,start_time_s=t0,end_time_s=t1,external_net_m3=block_net,
                outlet_m3=block_out,nonoutlet_external_m3=block_net+block_out,external_inflow_m3=block_in,
                precipitation_m=block_rain,solar_Wh_m2=block_solar,
                phase_end_balance_residual_m3=float(balance[-1]),phase_block_max_balance_residual_m3=float(np.max(np.abs(balance))),
                state_continuity_error_m3=state_error,actual_opening_continuity_error_m=gate_error))
            self.last_end=t1;self.last_state=x[:,-1].copy();self.last_gate=inputs[:2,-1]*self.max_open
            return t,x,dx,y
        return observed

    def finish(self,result,frame,antecedent_hours,tolerance):
        history=result.get('controller_history',[])
        before=next((h for h in history if h['block_end_s']==self.offset and h['completed_block']),None)
        after=next((h for h in history if h['time_s']==self.offset),None)
        if self.boundary['reached']:
            if before is None or after is None:raise RuntimeError('Reached boundary lacks the uninterrupted controller history')
            self.boundary.update(latches_before=before['enabled'],latches_before_next_update=after['previous_enabled'],
                latches_after_boundary_update=after['enabled'],latch_memory_exact=before['enabled']==after['previous_enabled'],
                raw_actual_opening_before_m=before['attained_end_opening_m'],raw_actual_opening_after_m=after['start_opening_m'],
                raw_actual_opening_exact=before['attained_end_opening_m']==after['start_opening_m'],
                previous_block_end_s=before['block_end_s'],next_update_time_s=after['time_s'],
                next_planned_end_opening_m=after['end_opening_m'],
                interpretation='Boundary update may legitimately change a latch from the carried state; compare previous_enabled, not post-update enabled')
            if not all(self.boundary[k] for k in ('state_exact','actual_opening_exact','latch_memory_exact','raw_actual_opening_exact')):
                raise RuntimeError('State, actuator or latch memory discontinuity at event boundary')
        for key,p in self.phases.items():
            segment=frame.iloc[:antecedent_hours] if key=='antecedent' else frame.iloc[antecedent_hours:]
            p['source_precipitation_m']=float(segment.precipitation_mm_per_h.sum()/1000.)
            p['source_solar_Wh_m2']=float(segment.solar_W_per_m2.sum())
            p['complete']=bool(p['completed_until_s']==p['target_end_time_s'] and p['completed_blocks']>0)
            p['rain_error_m']=p['integrated_precipitation_m']-p['source_precipitation_m'] if p['complete'] else None
            p['solar_error_Wh_m2']=p['integrated_solar_Wh_m2']-p['source_solar_Wh_m2'] if p['complete'] else None
            p['forcing_conservation_passed']=bool(p['complete'] and abs(p['rain_error_m'])<=1e-12 and abs(p['solar_error_Wh_m2'])<=1e-7)
            p['relative_balance_tolerance']=float(tolerance)
            p['balance_scale_m3']=p['external_inflow_volume_m3'] if p['external_inflow_volume_m3']>0 else (sum(p['initial_storage_m3']) if p['initial_storage_m3'] else None)
            p['balance_scale_rule']='Phase external inflow volume; phase initial total storage only if no positive external inflow'
            p['relative_balance_residual']=p['maximum_balance_residual_m3']/p['balance_scale_m3'] if p['balance_scale_m3'] else None
            p['balance_passed']=bool(p['complete'] and p['relative_balance_residual'] is not None and p['relative_balance_residual']<tolerance)
            p['net_nonoutlet_external_volume_m3']=p['external_net_volume_m3']+p['outlet_volume_m3']
            p['rhs_observations']=self.observers[key].summary()
            p['safety_exceeded_at_successful_rhs']=any(v['exceeded'] for v in p['rhs_observations']['safety_successful_rhs'])
            p['scope']='Completed control-block accounting and phase-specific RHS stages; partial-phase quantities are not complete event outcomes'
        return self.phases


def run_one(point,event,antecedent_hours,case,dt=15.):
    if case not in ('wetland','no_wetland'):raise ValueError('Choose the existing eight-state or six-state case')
    if antecedent_hours not in (24,72):raise ValueError('Only the approved 24 h/72 h conditioning durations are active')
    if event['event_id'] not in ('historical_20190825T0200Z','historical_20240804T1500Z'):
        raise ValueError('Only the two approved historical events are active')
    if point['sample_id'] not in ('mc_00016','mc_00024','mc_00018'):
        raise ValueError('Reuse only the three approved existing pilot vectors')
    original_plan_path=ROOT/'D/d2_analysis/parameter_uncertainty/results/study_plan.json'
    original_points={p['sample_id']:p for p in json.loads(original_plan_path.read_text())['mc']}
    approved=original_points[point['sample_id']]
    if point['physical_values']!=approved['physical_values'] or any(point[k]!=approved[k] for k in ('m_K_u','m_routing','m_R_r')):
        raise ValueError('The approved saved parameter vector was altered or resampled')
    forcing=copy.deepcopy(provider.record_for(event['event_id'],antecedent_hours))
    source=_resolve(forcing['forcing_file'])
    if frozen.sha(source)!=forcing['forcing_sha256']:raise ValueError('Prepared combined forcing bytes changed')
    if forcing['original_event_forcing_sha256']!=event['forcing_sha256']:raise ValueError('Original event forcing provenance changed')
    offset=float(antecedent_hours*3600);duration=float((antecedent_hours+EVENT_HOURS)*3600)
    if forcing['event_offset_s']!=offset or forcing['total_duration_s']!=duration:raise ValueError('Forcing manifest duration/offset mismatch')
    base,records=parameter_adapter._private_runner(point['physical_values'],dt)
    old_historical=base.__globals__['historical'];original_module=old_historical.adapter()
    module=SimpleNamespace(**vars(original_module))
    for name in ('read_forcing','hourly_input_segments','forcing_conservation'):setattr(module,name,getattr(provider,name))
    observed=[]
    def load_extended(*args,**kwargs):
        # Invoke the unchanged 144 h template loader first; do not change its
        # historical schedule-duration guard or any original cached object.
        fixed=old_historical.load_fixed(*args,**kwargs)
        fixed['model_config']=copy.deepcopy(fixed['model_config'])
        fixed['manifest']=copy.deepcopy(fixed['manifest'])
        fixed['official']=dict(fixed['official'])
        fixed['model_config']['numerics']['duration_s']=duration
        fixed['model_config']['geography']['event_start_utc']=forcing['start_utc']
        fixed['manifest'].update(duration_s=duration,conditioning_duration_s=offset,event_duration_s=EVENT_HOURS*3600.,
            initial_condition_note=CONDITIONING,schedule_description='Unchanged 120 s feedback policy, continuous across antecedent/event boundary')
        audit=BlockAudit(case,fixed,original_module,offset);observed.append(audit)
        fixed['official']['ode4u']=audit.wrap(fixed['official']['ode4u'])
        return fixed
    proxy=SimpleNamespace(adapter=lambda:module,load_fixed=load_extended,RHSObserver=old_historical.RHSObserver)
    namespace=dict(base.__globals__);namespace['historical']=proxy
    weather_record=dict(scenario_id=forcing['conditioning_id'],forcing_file=str(source),forcing_sha256=forcing['forcing_sha256'],seed_sequence_entropy=None)
    namespace['original_weather']=lambda:[copy.deepcopy(weather_record)]
    private=FunctionType(base.__code__,namespace,base.__name__,base.__defaults__,base.__closure__)
    private.__kwdefaults__=copy.copy(base.__kwdefaults__)
    result=private(case=case,scenario_id=forcing['conditioning_id'],design='best',dt=float(dt),hydraulics='D',Cw_US=2.8)
    if len(records)!=1 or len(observed)!=1:raise RuntimeError('Expected one private physical configuration and block observer')
    d=result['diagnostics'];record=copy.deepcopy(records[0]);audit=observed[0]
    frame,_=provider.read_forcing(source)
    if pd.Timestamp(d['event_start_utc'])!=pd.Timestamp(forcing['start_utc']):raise RuntimeError('Applied continuous run start differs from forcing')
    expected=record['initial_wetland_m3' if case=='wetland' else 'initial_no_wetland_m3']
    if d['initial_storage_m3']!=expected:raise RuntimeError('Accepted initial physical water was changed at antecedent start')
    if record['effective_physical_values']!=point['physical_values']:raise RuntimeError('Applied parameter vector differs from the approved saved vector')
    record['nominal_2022_template_provenance']={k:record[k] for k in ('deterministic_weather_id','prepared_forcing_path','prepared_forcing_sha256','prepared_forcing_manifest_sha256','adapter_sha256')}
    record.update(deterministic_weather_id=forcing['conditioning_id'],prepared_forcing_path=str(source.relative_to(ROOT)),
        prepared_forcing_sha256=forcing['forcing_sha256'],prepared_forcing_manifest_sha256=frozen.sha(HERE/'forcing/manifest.json'),
        adapter_sha256=frozen.sha(__file__),parameter_template_adapter_sha256=frozen.sha(parameter_adapter.__file__),
        conditioning=CONDITIONING,conditioning_start_utc=forcing['start_utc'],event_start_utc=forcing['event_start_utc'],
        event_end_utc=forcing['end_utc'],duration_s=duration,event_offset_s=offset,
        forcing_method='Separate unchanged hourly antecedent/event mappings; shifted grids preserve original event endpoint-hold inputs',
        weather_seed=None,weather_source_kind='Existing station-derived/interpolated weather and grid-estimated solar, with retained QC/trace flags')
    d['parameter_evaluation']=record;d['fixed_configuration']['parameter_evaluation']=copy.deepcopy(record)
    d['fixed_configuration']['conditional_on']=CONDITIONING
    tolerance=float(d['balance_relative_tolerance'] if 'balance_relative_tolerance' in d else audit.module._parameters(audit.fixed)[0].balance_rtol)
    phases=audit.finish(result,frame,antecedent_hours,tolerance)
    d['phase_diagnostics']=phases;d['event_boundary_handoff']=audit.boundary
    d['antecedent_experiment']=dict(event_id=event['event_id'],sample_id=point['sample_id'],antecedent_hours=antecedent_hours,
        event_offset_s=offset,event_duration_s=EVENT_HOURS*3600.,total_duration_s=duration,
        forcing=forcing,original_parameter_plan_sha256=frozen.sha(original_plan_path),continuous_controller_instance=True,no_state_or_controller_reset=True,
        maximum_completed_block_state_continuity_error_m3=audit.maximum_state_continuity_error_m3,
        maximum_completed_block_gate_continuity_error_m=audit.maximum_gate_continuity_error_m,
        failed_phase=None if d.get('full_integration_horizon_completed') else
            ('antecedent' if (d.get('first_failure_time_s') or d.get('completed_until_s',0.))<offset else 'event'),
        physical_peak_scope='Original runner peaks cover the entire extended run; event metrics must be derived from t>=event_offset_s, retaining that event-start state',
        residual_scope='Event residual equals full_run_balance(t)-full_run_balance(event_offset), with its own event inflow scale')
    if result['valid'] and not all(p['balance_passed'] and p['forcing_conservation_passed'] for p in phases.values()):
        d['classification']='numerical_or_state_failure';d['failure_reason']='phase_specific_balance_or_forcing_conservation_failure'
        d['failure_stage']='post_integration_phase_validation';d['full_valid_trajectory_available']=False
        d['physical_river_peak_m']=None;d['physical_river_peak_discharge_m3_s']=None
        d['physical_river_peak_time_s']=None;d['physical_river_peak_discharge_time_s']=None
        result['valid']=False;result['solver_status']='numerical_failure'
    result['block_audit']=audit.rows
    return result
