"""Focused checks of stateful policy semantics and pre-run physical margins."""
from pathlib import Path
import json
import numpy as np
from math import sqrt
from .policy import Controller, CONTROL_PERIOD_S
from D.g5_diagnostics import runner
from D.gate_hydraulics_revision import run_study
from D.gate_hydraulics_revision.hydraulics import gate_flow, Assumptions, ApplicabilityError

HERE=Path(__file__).resolve().parent


def run():
    checks=[]
    def ck(name,value,**evidence):checks.append(dict(check=name,passed=bool(value),**evidence))
    c=Controller.from_initial([.0015,.0015]); cfg=c.config
    first=c.update([.3,.3],0)
    ck('Original actual initial opening retained',first['start_opening_m']==[.0015,.0015])
    ck('Initial command is not instantaneous actual motion',first['end_opening_m']==[.0615,.0615] and first['target_opening_m']==[.1,.1])
    ck('Fixed clock independent of ODE grid',CONTROL_PERIOD_S==120)
    try:c.update([.3,.3],60);rejected=False
    except ValueError:rejected=True
    ck('RHS-like off-clock update rejected',rejected)
    on=Controller.from_initial([.1,.1]);on.update([.25,.25],0)
    off=Controller.from_initial([.1,.1]);off.update([.2,.2],0)
    on_plan=on.update([.25,.25],120);off_plan=off.update([.25,.25],120)
    ck('Hysteresis changes command at identical intermediate head',on_plan['target_opening_m']==[.1,.1] and off_plan['target_opening_m']==[0.,0.])
    ck('Reopen occurs at upper threshold',off.update([.3,.3],240)['enabled']==[True,True])
    ck('Close occurs at lower threshold',on.update([.2,.2],240)['enabled']==[False,False])
    close=Controller.from_initial([.3,.3]);plans=[close.update([.02,.02],k*120) for k in range(5)]
    ck('No instantaneous hydraulic-avoidance clipping',plans[0]['end_opening_m']==[.24,.24])
    ck('Finite full stroke reaches exact shutoff',plans[-1]['end_opening_m']==[0.,0.])
    for j,plan in enumerate(plans):
        ck(f'Closing block {j} rate and endpoint bounds',max(abs(v) for v in plan['rate_m_per_s'])<=.0005+1e-15 and all(0<=v<=.3 for v in plan['end_opening_m']))
    previous=None
    for j,h in enumerate([.1,.15,.2,.25,.3,.4,.6,.8,1.2,2.]):
        ctr=Controller.from_initial([.0015,.0015])
        row=ctr.update([h,h],0)
        target=row['target_opening_m'][0]
        ck(f'Command bounds at head {h}',0<=target<=.3)
        if target>0:ck(f'Commanded energy-gate ratio at head {h}',h/target>=2-1e-12,ratio=h/target)
        if previous is not None:ck(f'Enabled target monotonic head {h}',target>=previous)
        previous=target
    fixed=runner.load_fixed();ad=runner.adapter();p,x0,gates=ad._parameters(fixed)
    state=x0.copy();initial_heads=state[[4,5]]/p.A_d+p.z_d-np.array([g.gate_sill_z for g in gates])
    ck('Nominal initial heads support explicit ON initialization',np.max(abs(initial_heads-.3))<1e-12)
    # ET bound is checked against actual saved forcings, not used by control.
    scenarios=run_study.weather_manifest();source_t,source_u=ad.wet.fixed_forcing(fixed['model_config'],p)
    demand=[ad.a2.evaporative_demand(float(t),float(r),float(w),p) for t,r,w in zip(source_u[1],source_u[2],source_u[3])]
    forcings=0
    for item in scenarios:
        if not item['original_liquid_domain_admissible']:continue
        frame,_=ad.read_forcing(item['forcing_file']);forcings+=1
        for t,r,w in zip(frame['temperature_C'],frame['solar_W_per_m2'],frame['wind_m_per_s']):
            demand.append(ad.a2.evaporative_demand(float(t),float(r),float(w),p))
    # Interpolation can change demand nonlinearly. Also bound it using the
    # componentwise maximum positive input values across the saved warm data.
    maxima=np.array([np.max(source_u[1]),np.max(source_u[2]),np.max(source_u[3])])
    minimum_temperature=float(np.min(source_u[1]))
    for item in scenarios:
        if item['original_liquid_domain_admissible']:
            frame,_=ad.read_forcing(item['forcing_file'])
            maxima=np.maximum(maxima,frame[['temperature_C','solar_W_per_m2','wind_m_per_s']].max().to_numpy())
            minimum_temperature=min(minimum_temperature,float(frame['temperature_C'].min()))
    # Each weighting fraction in the reduced ET expression is <=1. Bound both
    # positive terms separately; do not assume the combined expression is
    # monotone in temperature between sampled observations.
    es_max=610.8*np.exp(17.27*maxima[0]/(maxima[0]+237.3))
    envelope_demand=float(p.eta*maxima[1]/(p.rho_w*p.lambda_v)+
        (p.k_0+p.k_v*maxima[2])*(1-p.r_H)*es_max/(p.R_v*(minimum_temperature+273.15)*p.rho_w))
    ck('ET envelope exceeds sampled demand',envelope_demand>=max(demand)-1e-15)
    total_time=cfg['control_period_s']+cfg['maximum_opening_m']/cfg['closing_rate_limit_m_per_s']
    buffer=[]
    for i,gate in enumerate(gates):
        ideal=sqrt(gate.gravity)*(2/3)**1.5
        draw=(ideal*gate.gate_width*.2**1.5/p.A_d[i]+p.K_ETd[i]*envelope_demand)*total_time
        buffer.append(dict(branch=i+1,head_at_close_m=.2,unsupported_weir_head_m=.05,
            delay_plus_closure_s=total_time,ideal_release_plus_ET_drawdown_bound_m=float(draw),
            remaining_head_margin_m=float(.2-.05-draw)))
        ck(f'Branch {i+1} model-based low-head buffer including saved-weather ET',draw<.15,**buffer[-1])
        capacities=[]
        for cw in (2.6,2.8,3.):
            capacities.append(dict(Cw_US=cw,Q_m3_s=gate_flow(gate.gate_sill_z+1.2,97.,1.,gate,Assumptions(Cw_US=cw)).discharge))
        ck(f'Branch {i+1} finite full-open safe-head capacities',all(v['Q_m3_s']>0 for v in capacities),capacities=capacities)
    report=dict(passed=all(r['passed'] for r in checks),check_count=len(checks),checks=checks,
        policy=cfg,evaporation_screen=dict(nonfreezing_saved_scenarios=forcings,
            largest_sampled_demand_m_s=float(max(demand)),componentwise_input_envelope=maxima.tolist(),
            minimum_temperature_C=minimum_temperature,envelope_demand_m_s=envelope_demand,
            derivation='Bound each positive ET term using weighting fractions <=1; maximum radiation/wind/vapor pressure and minimum absolute temperature',
            interpretation='Model-based deterministic upper-envelope screen; not a physical field guarantee'),
        scope='Controller consistency and physical-margin screen; actual ODE safety and numerical sensitivity are tested separately')
    (HERE/'results').mkdir(exist_ok=True)
    (HERE/'results/policy_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(passed=report['passed'],checks=len(checks),failed=[r['check'] for r in checks if not r['passed']])))
    return report


if __name__=='__main__':raise SystemExit(0 if run()['passed'] else 1)
