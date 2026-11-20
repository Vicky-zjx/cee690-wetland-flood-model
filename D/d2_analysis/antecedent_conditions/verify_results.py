"""Read-only audit of the 12 approved conditioning pairs and numerical repeat.

No ODE, parameter sampling, model fitting or source/output mutation. A fresh
controller instance is replayed only against saved accepted water states.
Print the report to stdout; redirect outside frozen result files if desired.
"""
from __future__ import annotations
from pathlib import Path
import hashlib
import json
import sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
import numpy as np
import pandas as pd
from D.feedback_control.policy import Controller
from D.feedback_control import run_policy as frozen
from D.d2_analysis.antecedent_conditions import preservation

LIMITS={'H_d1':101.2,'H_d2':101.,'H_w1':98.9,'H_w2':98.8}
CASES=('wetland','no_wetland')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check():
    plan_path=HERE/'study_plan.json';plan=json.loads(plan_path.read_text())
    for relative,digest in plan['source_sha256'].items():assert sha(ROOT/relative)==digest,relative
    old_plan_path=ROOT/'D/d2_analysis/parameter_uncertainty/results/study_plan.json'
    assert sha(old_plan_path)==plan['source_parameter_plan_sha256']
    original={p['sample_id']:p for p in json.loads(old_plan_path.read_text())['mc']}
    assert plan['parameter_sampling_calls']==plan['weather_generation_calls']==0
    assert plan['unique_physical_pairs']==12 and plan['antecedent_hours']==[24,72]
    for point in plan['points']:assert point==original[point['sample_id']]
    module=frozen.historical.adapter();fixed=frozen.historical.load_fixed(dt=15,design='best')
    params,initial8,gates=module._parameters(fixed)
    gate_max=np.asarray([g.gate_max_opening for g in gates]);sills=np.asarray([g.gate_sill_z for g in gates])
    audit_rows=[];pair_rows=[];controller_checks=0;largest_gate_replay_error=0.;largest_state_seam_error=0.
    largest_phase_net_error=largest_phase_residual_error=largest_block_sum_error=0.
    max_phase_relative=0.;min_safety=np.inf;min_ratio=np.inf;max_saved_derivative_mass_error=0.
    for event in plan['events']:
        for hours in plan['antecedent_hours']:
            for point in plan['points']:
                sid=point['sample_id'];eid=event['event_id'];directory=HERE/'results'/f'{eid}__{hours}h__{sid}'
                pilot=(eid==plan['pilot']['event_id'] and hours==plan['pilot']['antecedent_hours'] and sid==plan['pilot']['sample_id'])
                for dt in ([15.,7.5] if pilot else [15.]):
                    tag=f'{dt:g}'.replace('.','p')+'s';loaded={}
                    for case in CASES:
                        name=f'{sid}__{case}__{tag}';path=directory/(name+'.json')
                        assert path.is_file(),f'Missing planned result: {path}'
                        d=json.loads(path.read_text());valid=d.get('full_valid_trajectory_available',False)
                        rec=dict(event_id=eid,antecedent_hours=hours,sample_id=sid,case=case,dt_s=dt,
                            classification=d['classification'],valid=valid,unsafe=d.get('safety_violation'),
                            full_horizon_completed=d.get('full_integration_horizon_completed'),
                            failure_reason=d.get('failure_reason') or d.get('failure_message'),
                            diagnostic_sha256=sha(path),trajectory_checks_available=False)
                        audit_rows.append(rec)
                        assert d['antecedent_study_plan_sha256']==sha(plan_path)
                        if not valid:
                            assert rec['failure_reason'],'Invalid result needs a recorded cause'
                            continue
                        assert d['full_integration_horizon_completed'] and d['completed_until_s']==(hours+144)*3600.
                        record=d['parameter_evaluation'];assert record==d['fixed_configuration']['parameter_evaluation']
                        assert record['effective_physical_values']==point['physical_values']
                        assert record['policy_sha256']==sha(frozen.POLICY_PATH)
                        assert pd.Timestamp(d['event_start_utc'])+pd.Timedelta(hours=hours)==pd.Timestamp(event['start_utc'])
                        for relative,digest in d['fixed_configuration']['stage_d_feedback_source_sha256'].items():assert sha(ROOT/relative)==digest
                        start=initial8 if case=='wetland' else initial8[module.dry.SHARED]
                        assert np.array_equal(d['initial_storage_m3'],start)
                        assert record['initial_actual_opening_m']==[.0015,.0015]
                        assert record['controller_period_s']==120. and d['control_period_s']==120.
                        with np.load(directory/(name+'.npz'),allow_pickle=False) as saved:z={k:saved[k].copy() for k in saved.files}
                        t=z['time'];x=z['x'];u=z['u'];y=z['y'];balance=z['balance']
                        names=list(z['output_names']);states=list(z['state_names']);Y={n:i for i,n in enumerate(names)}
                        assert t[0]==0. and t[-1]==(hours+144)*3600. and np.all(np.diff(t)==dt)
                        assert np.array_equal(x[:,0],start) and np.array_equal(x[:,-1],d['last_completed_state_m3'])
                        assert all(np.isfinite(z[k]).all() for k in ('x','dx','y','u','balance'))
                        external=u[module.a2.U['P']]*params.A_rain+u[module.a2.U['Q_ext']]+y[Y['Q_gx']]-y[Y['Q_out']]-y[Y['E_tot']]
                        derivative_error=float(np.max(abs(z['dx'].sum(axis=0)-external)))
                        max_saved_derivative_mass_error=max(max_saved_derivative_mass_error,derivative_error)
                        assert derivative_error<1e-12 and np.max(abs(y[Y['external_net']]-external))<1e-12
                        assert x.min()>=0 and np.max(x[states.index('S_m')])<=params.C_m
                        boundary=hours*3600.;k=int(np.searchsorted(t,boundary));assert t[k]==boundary
                        prior_path=ROOT/'D/d2_analysis/historical_storm_robustness/results/pilot'/eid/(name+'.npz')
                        with np.load(prior_path,allow_pickle=False) as prior:
                            assert np.array_equal(t[k:]-boundary,prior['time'])
                            assert np.array_equal(u[2:,k:],prior['u'][2:]),'Original event weather inputs changed'
                        handoff=d['event_boundary_handoff']
                        assert all(handoff[q] for q in ('reached','state_exact','actual_opening_exact','latch_memory_exact','raw_actual_opening_exact'))
                        assert np.array_equal(x[:,k],handoff['state_before_m3']) and np.array_equal(x[:,k],handoff['state_after_m3'])
                        largest_state_seam_error=max(largest_state_seam_error,float(np.max(abs(np.asarray(handoff['state_before_m3'])-handoff['state_after_m3']))))
                        assert handoff['previous_block_end_s']==handoff['next_update_time_s']==boundary
                        assert handoff['latches_before']==handoff['latches_before_next_update']
                        ctl=pd.read_csv(directory/(name+'_controller.csv.gz'),float_precision='round_trip')
                        blocks=pd.read_csv(directory/(name+'_blocks.csv.gz'),float_precision='round_trip')
                        count=int(t[-1]/120.);assert len(ctl)==2*count and len(blocks)==count
                        assert np.array_equal(blocks.start_time_s,np.arange(count)*120.)
                        assert np.array_equal(blocks.end_time_s,(np.arange(count)+1)*120.)
                        assert np.array_equal(ctl.time_s.to_numpy().reshape(-1,2),np.repeat((np.arange(count)*120.)[:,None],2,axis=1))
                        assert np.array_equal(ctl.branch.to_numpy().reshape(-1,2),np.tile([1,2],(count,1)))
                        assert ctl.completed_block.all() and blocks.state_continuity_error_m3.eq(0).all() and blocks.actual_opening_continuity_error_m.eq(0).all()
                        controller=Controller.from_initial([.0015,.0015]);block_samples=int(120/dt)
                        numeric_fields={'start_opening_m':'start_opening_m','target_opening_m':'target_opening_m',
                            'planned_end_opening_m':'end_opening_m','rate_m_per_s':'rate_m_per_s','attained_end_opening_m':'end_opening_m'}
                        numeric={c:ctl[c].to_numpy(float).reshape(-1,2) for c in numeric_fields}
                        enabled=ctl.enabled.to_numpy().reshape(-1,2);changed=ctl.latch_changed.to_numpy().reshape(-1,2)
                        modes=ctl['mode'].to_numpy().reshape(-1,2)
                        for j in range(count):
                            index=j*block_samples;heads=x[[states.index('S_d1'),states.index('S_d2')],index]/params.A_d+params.z_d-sills
                            command=controller.update(heads,j*120.)
                            for column,key in numeric_fields.items():assert np.array_equal(numeric[column][j],command[key]),(eid,hours,sid,case,j,column)
                            assert np.array_equal(enabled[j],command['enabled']) and np.array_equal(changed[j],command['latch_changed'])
                            assert np.array_equal(modes[j],command['modes'])
                            segment=slice(index,index+block_samples+1);elapsed=(t[segment]-j*120.)/120.
                            expected_opening=np.asarray(command['start_opening_m'])[:,None]+(np.asarray(command['end_opening_m'])-command['start_opening_m'])[:,None]*elapsed
                            error=float(np.max(abs(u[:2,segment]*gate_max[:,None]-expected_opening)))
                            largest_gate_replay_error=max(largest_gate_replay_error,error);assert error<1e-14
                            controller_checks+=2
                        unsafe=False
                        for phase,selection in [('antecedent',slice(0,k+1)),('event',slice(k,None))]:
                            p=d['phase_diagnostics'][phase];bp=blocks[blocks.phase.eq(phase)]
                            assert p['complete'] and p['balance_passed'] and p['forcing_conservation_passed']
                            assert p['relative_balance_tolerance']==plan['numerics']['mass_balance_relative_tolerance']==1e-4
                            xx=x[:,selection];bb=balance[selection]-balance[selection][0]
                            assert np.array_equal(xx[:,0],p['initial_storage_m3']) and np.array_equal(xx[:,-1],p['final_storage_m3'])
                            net=float(xx[:,-1].sum()-xx[:,0].sum()-bb[-1])
                            err=abs(net-p['external_net_volume_m3']);largest_phase_net_error=max(largest_phase_net_error,err);assert err<1e-7
                            residual_error=abs(np.max(abs(bb))-p['maximum_balance_residual_m3'])
                            largest_phase_residual_error=max(largest_phase_residual_error,float(residual_error));assert residual_error<1e-7
                            assert abs(bb[-1]-p['final_balance_residual_m3'])<1e-7
                            for col,key in [('external_net_m3','external_net_volume_m3'),('outlet_m3','outlet_volume_m3'),
                                ('external_inflow_m3','external_inflow_volume_m3'),('precipitation_m','integrated_precipitation_m'),('solar_Wh_m2','integrated_solar_Wh_m2')]:
                                error=abs(float(bp[col].sum())-p[key]);largest_block_sum_error=max(largest_block_sum_error,error);assert error<1e-7
                            out=float(np.trapezoid(z['downstream_discharge'][selection],t[selection]));assert abs(out-p['outlet_volume_m3'])<1e-7
                            assert abs(p['net_nonoutlet_external_volume_m3']-(p['external_net_volume_m3']+p['outlet_volume_m3']))<1e-9
                            scale=p['external_inflow_volume_m3'] if p['external_inflow_volume_m3']>0 else sum(p['initial_storage_m3'])
                            assert scale==p['balance_scale_m3'] and p['relative_balance_residual']==p['maximum_balance_residual_m3']/scale
                            max_phase_relative=max(max_phase_relative,p['relative_balance_residual'])
                            observation=p['rhs_observations'];phase_unsafe=False
                            for safety in observation['safety_successful_rhs']:
                                variable=safety['variable'];assert safety['threshold_m']==LIMITS[variable]
                                peak=max(float(y[Y[variable],selection].max()),safety['peak_m'])
                                margin=LIMITS[variable]-peak;min_safety=min(min_safety,margin);phase_unsafe|=margin<0
                            assert phase_unsafe==p['safety_exceeded_at_successful_rhs'];unsafe|=phase_unsafe
                            for ratio in p['minimum_open_gate_H_over_a']:
                                if ratio is not None:min_ratio=min(min_ratio,ratio)
                            assert observation['rhs_failures']==0
                        assert unsafe==bool(d['safety_violation'])
                        all_net=float(blocks.external_net_m3.sum());all_out=float(blocks.outlet_m3.sum())
                        assert abs(all_net-d['completed_blocks_mass_diagnostics']['external_net_volume_m3'])<1e-7
                        assert abs(all_out-float(np.trapezoid(z['downstream_discharge'],t)))<1e-7
                        assert abs(float(x[:,-1].sum()-x[:,0].sum())-all_net-balance[-1])<1e-7
                        rec.update(trajectory_checks_available=True,parameter_vector_exact=True,event_weather_exact=True,
                            boundary_exact=True,controller_replay_branch_checks=2*count,full_balance_relative=d['balance_relative'],
                            antecedent_balance_relative=d['phase_diagnostics']['antecedent']['relative_balance_residual'],
                            event_balance_relative=d['phase_diagnostics']['event']['relative_balance_residual'])
                        loaded[case]=(d,z)
                    pair=dict(event_id=eid,sample_id=sid,antecedent_hours=hours,dt_s=dt,both_valid=len(loaded)==2)
                    if len(loaded)==2:
                        dw,zw=loaded['wetland'];dn,zn=loaded['no_wetland']
                        assert np.array_equal(zw['time'],zn['time']) and np.array_equal(zw['u'][2:],zn['u'][2:])
                        assert np.array_equal(zw['x'][module.dry.SHARED,0],zn['x'][:,0])
                        pair['maximum_within_pair_actual_gate_difference_m']=float(np.max(abs((zw['u'][:2]-zn['u'][:2])*gate_max[:,None])))
                        pair['paired_weather_exact']=True
                    pair_rows.append(pair)
    assert len(audit_rows)==26 and len(pair_rows)==13 and sum(r['dt_s']==15. for r in pair_rows)==12
    pilot=json.loads((HERE/'results/pilot_summary.json').read_text());assert pilot['passed']
    guard=preservation.verify();assert guard['passed'],guard
    return dict(status='pass',physical_ODE_runs_performed_by_audit=0,controller_only_replay_branch_checks=controller_checks,
        primary_physical_pairs=12,numerical_refinement_pairs=1,retained_model_records=len(audit_rows),
        valid_model_records=sum(r['valid'] for r in audit_rows),invalid_model_records=sum(not r['valid'] for r in audit_rows),
        unsafe_valid_model_records=sum(r['valid'] and bool(r['unsafe']) for r in audit_rows),
        maximum_gate_replay_error_m=largest_gate_replay_error,maximum_state_seam_error_m3=largest_state_seam_error,
        maximum_phase_net_identity_error_m3=largest_phase_net_error,maximum_phase_rebased_residual_difference_m3=largest_phase_residual_error,
        maximum_block_integral_sum_difference=largest_block_sum_error,maximum_phase_relative_mass_residual=max_phase_relative,
        minimum_safety_margin_m=min_safety,minimum_open_gate_H_over_a=min_ratio,
        maximum_saved_derivative_mass_identity_error_m3_s=max_saved_derivative_mass_error,
        all_original_parameter_vectors_reused=True,event_input_and_time_mapping_preserved=True,
        no_full_joined_discontinuous_flux_integration=True,
        study_plan_sha256=sha(plan_path),verification_source_sha256=sha(__file__),
        rows=audit_rows,pairs=pair_rows,numerical_pilot=pilot,preservation=guard)


if __name__=='__main__':print(json.dumps(check(),indent=2,allow_nan=False))
