#!/usr/bin/env python3
"""HW B.2 direct shooting. Missing approvals block project optimization.

A.2 remains untouched; only the two gate transfers are replaced in the B.2 RHS.
Run --preflight first. --plumbing-check runs the original synthetic A.2 model,
not a B.2 constrained solution. See b2/docs/02_model_and_optimization.md.
"""
from __future__ import annotations
import argparse
from collections import Counter, OrderedDict
import copy
import csv
from dataclasses import asdict
import hashlib
import itertools
import json
from pathlib import Path
from B.project_paths import resolve_path
import time as clock
import numpy as np
from A import wetland_flood_model as a2

ROOT = Path(__file__).resolve().parent
CONFIG = ROOT / 'b2_config.json'
CONSTRAINT_COUNT = 20  # 18 existing residuals plus two approved gate applicability constraints

class PreflightError(ValueError):
    pass

class TrajectoryFailure(RuntimeError):
    """Expected physical-domain or integration failure, not a programming error."""

class HydraulicApplicabilityFailure(TrajectoryFailure):
    """Stop at first encountered invalid RK stage; no fabricated continuation."""
    def __init__(self, diagnostic):
        self.diagnostic=diagnostic
        super().__init__('Approved gate head/opening constraint violated: '+json.dumps(diagnostic))

class RunBudgetExceeded(RuntimeError):
    pass


def serial(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, Path): return str(value)
    raise TypeError(type(value).__name__)


def save_json(path, data):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, default=serial, allow_nan=False)+'\n')


def write_table(name, rows):
    path = ROOT / 'tables' / name
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows: return
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def load_config(path=CONFIG):
    return json.loads(resolve_path(path).read_text())


def config_hash(config):
    return hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()


def base_parameters(config):
    params = a2.c()
    for name, record in config['parameter_overrides'].items():
        if name in {'gate1_values','gate2_values','rho','tau','duration','time_step'}:
            raise PreflightError('Decisions/time settings belong in their dedicated configuration fields: '+name)
        if not hasattr(params, name) or not isinstance(record, dict):
            raise PreflightError('Invalid parameter override: '+name)
        if not record.get('source_id') or not record.get('reason') or record.get('approved') is not True:
            raise PreflightError('Unapproved/unsourced parameter override: '+name)
        old = getattr(params, name)
        setattr(params, name, np.asarray(record['value'], dtype=float) if isinstance(old, np.ndarray) else record['value'])
    params.duration = float(config['numerics']['duration_s'])
    params.time_step = float(config['numerics']['dt_s'])
    return params


def gate_parameters(config):
    from B.gate_hydraulics import GateParameters
    out=[]
    for i in range(2):
        out.append(GateParameters(
            gravity=config['gates']['gravity'],
            gate_width=config['gates']['gate_width'][i],
            gate_max_opening=config['gates']['gate_max_opening'][i],
            gate_sill_z=config['gates']['gate_sill_z'][i],
            Cd_gate_free=config['gates']['Cd_gate_free'][i],
            Cd_gate_submerged=config['gates']['Cd_gate_submerged'][i],
            neglect_approach_velocity_head=config['gates']['neglect_approach_velocity_head']))
    return out


def reference_decision(config):
    return np.r_[config['controls']['gate1_values'], config['controls']['gate2_values'], config['reference']['covers']]


def decode(theta, config):
    theta=np.asarray(theta, dtype=float)
    n1=len(config['controls']['gate1_hours']); n2=len(config['controls']['gate2_hours'])
    if theta.shape != (n1+n2+2,) or not np.all(np.isfinite(theta)):
        raise ValueError('Decision vector must contain both gate grids followed by two covers.')
    if np.any(theta<0) or np.any(theta>1):
        raise TrajectoryFailure('Decision outside native [0,1] bounds; never clipped.')
    return theta[:n1],theta[n1:n1+n2],theta[-2:]


def preflight(config):
    """One grouped approval/gap table; no invented defaults for physical choices."""
    gaps=[]
    def gap(code,topic,missing,action):
        gaps.append(dict(id=code, topic=topic, status='blocked', missing=missing, required_action=action))
    if config['mode'] not in {'approved_illustrative','approved_local'}:
        gap('G1','Study unit / event / run label','Proposed observed-event-forced hypothetical0.66km2 installation is not approved.',
            'Approve the0.66km2 analytical installation and Sep9–15UTC2022 event; gauge46.879km2 is context only, Q_ext=0. See b2/docs/01_study_configuration.md.')
    if any(v is None for v in config['safety']['dam_m']+config['safety']['wetland_m']) or not config['safety']['approved'] or len(config['safety']['source_ids'])!=4 or not all(config['safety']['source_ids']):
        gap('G2','Four upstream safety levels','Proposed safety levels or their supporting definitions remain unapproved/incomplete.',
            'Approve dam[101.20,101.00]m and wetland[98.90,98.80]m in assumed modeldatum, derived from specifiedsettledtops minus0.35m; keepfixed regardless of feasibility.')
    if not config['categories']['approved'] or not all(config['categories']['admissible']):
        gap('G3','Vegetation categories and paired endpoints','The proposed lower/higher response classes are synthetic and unapproved, not species mappings.',
            'Approve paired type_a R=[2.8,3.0],K=[1.10,1.10]; type_b R=[4.0,4.5],K=[1.25,1.25], for bothbranches; exactspeciesproperties remain unsupported.')
    gate_fields=['gate_width','gate_max_opening','gate_sill_z','Cd_gate_free','Cd_gate_submerged']
    missing=[k for k in gate_fields if any(v is None for v in config['gates'][k])]
    if config['gates']['gravity'] is None: missing.append('gravity')
    if missing or not config['gates']['approved'] or config['gates']['neglect_approach_velocity_head'] is not True:
        gap('G4','Physical gate configuration','Missing/unapproved: '+', '.join(missing or ['geometry, coefficients or negligible-approach-head assumption']),
            'Approve widths[1,.8]m, travel[.3,.3]m, sills[100,99.8]m, Cd_free=.60,Cd_submerged=.80, gravity9.80665 and negligible approach-head assumption; hypothetical bed-mounted channel gates.')
    if config['gates']['limited_regimes_approved'] is not True:
        gap('G5','Unimplemented hydraulic regimes','Partial/weir opening and the weir-to-gate transition remain explicit unsupported regimes; submergence blending is implemented from USACE 7.0.',
            'Approve H_pool-z_sill >=1.25*a_max*g for each open/wet/forward-flow gate throughout integration; no clipping. Otherwise crest/transition method remains unresolved.')
    if not config['fixed_model']['approved'] or not config['forcing']['approved'] or not config['initial_conditions']['approved']:
        gap('G6','Fixed geometry, forcing, initial state and datum','Event forcing is prepared and fixed assumptions specified, but their use is not approved.',
            'Approve KMKE145.3mm quantifiedrain including explicittrace0 lowerbound; KMWC T/wind, r_H=.7084916612,p_0=98692.51172Pa; POWERhourstart/energy reconstruction; Q_ext=0,h_d=h_g=97.15m,B=0; fixedA2 x0 and remaining declared assumptions.')
    n=config['numerics']
    if config['optimizer']['feasibility_tol'] != n['feasibility_tol']:
        raise PreflightError('Optimizer and independent audit feasibility tolerances must match.')
    if not all(np.isfinite(v) for v in list(n['volume_scales_m3'])+[n['height_scale_m'],n['strict_margin_m'],n['feasibility_tol']]):
        raise PreflightError('Residual scales and numerical tolerances must be finite.')
    if n['height_scale_m']<=0 or min(n['volume_scales_m3'])<=0 or len(n['volume_scales_m3'])!=8:
        raise PreflightError('All fixed residual scales must be positive and dimensionally specified.')
    if n['strict_margin_m'] <= n['feasibility_tol']*n['height_scale_m']:
        raise PreflightError('Numerical strict-domain margin must exceed the feasibility tolerance in metres.')
    if n['dt_s']<=0 or n['refine_factor']<2 or int(n['refine_factor'])!=n['refine_factor']:
        raise PreflightError('Require positive dt and integer refinement factor >=2.')
    for branch in (1,2):
        times=np.asarray(config['controls'][f'gate{branch}_hours'])*3600
        if times[0]!=0 or times[-1]!=n['duration_s'] or np.any(np.diff(times)<=0):
            raise PreflightError('Gate knot times must increase and cover the fixed horizon exactly.')
        # The existing A.2 knot grid must be represented exactly on the ODE grid.
        if not np.allclose(times/n['dt_s'],np.round(times/n['dt_s']),rtol=0,atol=1e-10):
            raise PreflightError('This grid would approximate rather than retain the gate-knot interpolation.')
    decode(reference_decision(config),config)
    if not gaps:
        if not np.all(np.isfinite(config['safety']['dam_m']+config['safety']['wetland_m'])):
            raise PreflightError('Approved safety levels must be finite.')
        if not all(config['gates']['source_ids'].get(k) for k in gate_fields+['gravity','neglect_approach_velocity_head']):
            raise PreflightError('Every approved gate parameter/assumption requires a source ID.')
        params=base_parameters(config)
        gp=gate_parameters(config)
        for i,p in enumerate(gp):
            if p.gate_sill_z<params.z_d[i] or p.gate_sill_z>=params.z_dsp[i]:
                raise PreflightError('Gate sill must be at/above the pool bed and below the independent spill crest.')
        for i,labels in enumerate(config['categories']['admissible']):
            for label in labels:
                if label not in config['categories']['endpoints']:
                    raise PreflightError('Missing paired endpoints for '+label)
        if config['mode']=='approved_local':
            required=['polygon_path','polygon_sha256','outlet','crs','vertical_datum','event_start_utc']
            if not config['geography']['approved'] or any(not config['geography'].get(k) for k in required):
                raise PreflightError('Local run requires approved polygon/outlet/CRS/datum/time metadata.')
            if config['forcing']['kind']!='prepared_csv':
                raise PreflightError('A local-data run cannot use the A.2 synthetic forcing.')
            path=resolve_path(config['geography']['polygon_path'])
            if hashlib.sha256(path.read_bytes()).hexdigest()!=config['geography']['polygon_sha256']:
                raise PreflightError('Approved polygon checksum mismatch.')
        fixed_forcing(config,params)  # QC before any optimization
    return gaps


def fixed_forcing(config,params):
    """Return fixed forcing knots. Refinement interpolates the SAME knot values."""
    kind=config['forcing']['kind']
    if kind=='a2_synthetic_fixed_samples':
        return params.time.copy(), a2.make_inputs(params)[2:].copy()
    if kind!='prepared_csv': raise PreflightError('Unknown forcing kind.')
    path=resolve_path(config['forcing']['prepared_csv']); meta_path=resolve_path(config['forcing']['manifest_json'])
    meta=json.loads(meta_path.read_text())
    for field in ['source_ids','original_timezone','timestamp_convention','qc_summary','gap_policy','derivations','vertical_datum','Q_ext_accounting']:
        if not meta.get(field): raise PreflightError('Prepared forcing metadata missing '+field)
    if meta.get('timezone')!='UTC' or meta.get('approved') is not True:
        raise PreflightError('Prepared inputs must be approved and UTC.')
    if meta['vertical_datum']!=config['geography']['vertical_datum']:
        raise PreflightError('Boundary stages and geometry do not share a declared datum.')
    if hashlib.sha256(path.read_bytes()).hexdigest()!=meta['sha256']:
        raise PreflightError('Prepared forcing checksum mismatch.')
    data=np.genfromtxt(path,delimiter=',',names=True)
    names=['time_s']+a2.INPUT_NAMES[2:]
    if list(data.dtype.names)!=names: raise PreflightError('Prepared CSV columns must be '+','.join(names))
    t=np.asarray(data['time_s']); values=np.vstack([data[name] for name in names[1:]])
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(t)) or np.any(np.diff(t)<=0):
        raise PreflightError('Missing, nonfinite or unordered prepared forcing; never fill rainfall gaps with zero.')
    if t[0]!=0 or t[-1]!=params.duration: raise PreflightError('Forcing must cover exactly the configured horizon.')
    if not np.allclose(t/params.time_step,np.round(t/params.time_step),rtol=0,atol=1e-10):
        raise PreflightError('Forcing knots must be represented on the integration grid without dropping pulses.')
    if np.any(values[0]<0) or np.any(values[2:5]<0) or np.any(values[-1]<0) or np.any(values[-1]>1):
        raise PreflightError('Invalid precipitation/radiation/wind/inflow/blockage.')
    if abs(np.trapezoid(values[0],t)-meta['rain_depth_m'])>meta['rain_mass_tolerance_m']:
        raise PreflightError('Linearly interpolated rainfall does not preserve documented accumulated depth.')
    if abs(np.trapezoid(values[2],t)/3600-meta['solar_energy_Wh_m2'])>1e-6:
        raise PreflightError('Prepared radiation does not preserve documented solar energy.')
    return t,values


def make_b2_rhs(gates, regimes=None):
    from B.gate_hydraulics import sluice_gate_flow,UnsupportedGateRegime
    def rhs(t,x,u,params):
        # Original A.2 computes ALL other fluxes; zero only its two gate controls.
        u_without_gate=np.array(u,copy=True); u_without_gate[:2]=0.0
        dx,y=a2.system(t,x,u_without_gate,params)
        for name,bound in [('H_g',params.H_g_valid_max),('H_g',params.soil_bottom),('H_g',min(params.z_w)),('H_r',params.H_r_valid_max)]:
            if y[a2.Y[name]]>=bound: raise TrajectoryFailure('Applicability domain exceeded at t='+str(t)+': '+name)
        for i in range(2):
            upstream=float(y[a2.Y[f'H_d{i+1}']]); tail=float(y[a2.Y[f'H_w{i+1}']])
            try:
                flow=sluice_gate_flow(upstream,tail,u[i],gates[i])
            except UnsupportedGateRegime as exc:
                H=upstream-gates[i].gate_sill_z; opening=gates[i].gate_max_opening*float(u[i])
                raise HydraulicApplicabilityFailure(dict(branch=i+1,time_s=float(t),
                    water_depth_m=H,opening_height_m=opening,normalized_opening=float(u[i]),
                    residual_m=1.25*opening-H,upstream_stage_m=upstream,tailwater_stage_m=tail,
                    rule='1.25*a-H <= 0 when open, wet and forward-flow',regime=exc.regime,
                    detection='First encountered invalid ode4u RHS/RK stage; not a continuous-time crossing estimate')) from exc
            q=flow.discharge
            dx[4+i]-=q; dx[2+i]+=q
            y[a2.Y[f'Q_G{i+1}']]=q
            if regimes is not None: regimes[(i+1,flow.regime)]+=1
        return dx,y
    return rhs


def constraint_audit(time,x,y,u,theta,params,config):
    n=config['numerics']; hscale=n['height_scale_m']; tol=n['feasibility_tol']; rows=[]
    def add(label,physical,scale,units,threshold,source,kind):
        physical=np.asarray(physical); k=int(np.argmax(physical)); raw=float(physical[k]); scaled=raw/scale
        rows.append(dict(label=label,kind=kind,threshold=threshold,source_id=source,worst_time_s=float(time[k]),
                         margin=-raw,units=units,scale=scale,scaled_residual=scaled,physical_tolerance=tol*scale,
                         strict_pass=raw<=0,pass_at_tolerance=scaled<=tol))
    for j,name in enumerate(a2.STATE_NAMES):
        add('nonnegative_'+name,-x[j],n['volume_scales_m3'][j],'m3',0,'B1','B1 path')
    def parameter_source(name):
        return config['parameter_overrides'].get(name,{}).get('source_id','A2_SOURCE')
    add('soil_capacity',x[1]-params.C_m,n['volume_scales_m3'][1],'m3',params.C_m,parameter_source('C_m'),'B1 path')
    for i,name in enumerate(['H_d1','H_d2','H_w1','H_w2']):
        bound=(config['safety']['dam_m']+config['safety']['wetland_m'])[i]
        if bound is None: raise PreflightError('Missing safety threshold: '+name)
        add(name+'_safe',y[a2.Y[name]]-bound,hscale,'m',bound,config['safety']['source_ids'][i],'B1 upstream safety')
    for label,name,bound,parameter in [('groundwater_geometry','H_g',params.H_g_valid_max,'H_g_valid_max'),('groundwater_soil','H_g',params.soil_bottom,'soil_bottom'),('groundwater_wetland_bed','H_g',min(params.z_w),'z_w'),('river_geometry','H_r',params.H_r_valid_max,'H_r_valid_max')]:
        add(label,y[a2.Y[name]]-(bound-n['strict_margin_m']),hscale,'m',bound,parameter_source(parameter),'model applicability; numerical margin only')
    add('integration_status',np.full(len(time),-1.),1.,'1',0,'B2_CONFIG','numerical validity')
    for i in range(2):
        H=y[a2.Y[f'H_d{i+1}']]-config['gates']['gate_sill_z'][i]
        opening=config['gates']['gate_max_opening'][i]*u[i]
        active=(u[i]>0)&(H>0)&(y[a2.Y[f'H_d{i+1}']]>y[a2.Y[f'H_w{i+1}']])
        residual=np.where(active,1.25*opening-H,0.)
        add(f'gate{i+1}_head_opening',residual,hscale,'m',0,'B2_CONFIG',
            'additional approved hydraulic-applicability path constraint; inactive closed/dry/adverse cases have residual 0')
    residuals=np.array([r['scaled_residual'] for r in rows])
    assert len(residuals)==CONSTRAINT_COUNT
    # Native bounds checked independently after reconstruction; not penalty-weighted.
    for i in range(2):
        add(f'gate{i+1}_lower',-u[i],1.,'1',0,'B1','native/reconstructed bound')
        add(f'gate{i+1}_upper',u[i]-1,1.,'1',1,'B1','native/reconstructed bound')
        add(f'cover{i+1}_lower',np.full(len(time),-theta[-2+i]),1.,'1',0,'B1','native bound')
        add(f'cover{i+1}_upper',np.full(len(time),theta[-2+i]-1),1.,'1',1,'B1','native bound')
    return residuals,rows


def simulate(theta,pair,config,*,refine=1,forcing=None):
    params=base_parameters(config); params.time_step/=refine
    g1,g2,covers=decode(theta,config)
    params.rho=covers.copy(); params.tau=list(pair)
    endpoints=config['categories']['endpoints']
    for label in pair:
        params.R_veg[label]=np.array(endpoints[label]['R_veg'],dtype=float)
        params.K_veg[label]=np.array(endpoints[label]['K_veg'],dtype=float)
    gates=gate_parameters(config)
    # The NEW explicit sill is used for retained geometry-order checks only.
    params.z_G=np.array([g.gate_sill_z for g in gates])
    t= params.time; source_t,source_u=forcing if forcing is not None else fixed_forcing(config,base_parameters(config))
    u=np.vstack([np.interp(t,np.array(config['controls']['gate1_hours'])*3600,g1),
                 np.interp(t,np.array(config['controls']['gate2_hours'])*3600,g2),
                 *[np.interp(t,source_t,row) for row in source_u]])
    initial= np.asarray(config['initial_conditions']['x0_m3'],dtype=float)
    if initial.shape!=(8,) or not np.all(np.isfinite(initial)): raise PreflightError('Eight fixed finite initial storages required.')
    params.u_depth_init=initial[0]/params.A_u
    params.soil_fraction_init=initial[1]/params.C_m
    params.w_depth_init=initial[2:4]/params.A_w
    params.d_depth_init=initial[4:6]/params.A_d
    params.H_g_init=a2.stage(initial[6],params.S_y*params.A_g,params.z_g)
    params.r_depth_init=initial[7]/params.A_r
    a2.validate(params,u)  # identical fixed x0; preserved geometry/order checks
    if np.any(initial<0) or initial[1]>params.C_m: raise PreflightError('Fixed initial state outside storage domain.')
    a2.check_groundwater_levels(a2.stage(initial[6],params.S_y*params.A_g,params.z_g),params,'B.2 fixed x0')
    regimes=Counter(); ode=a2.load_ode4u()
    t,x,dx,y=ode(make_b2_rhs(gates,regimes),t,initial.copy(),u,params)
    if not all(np.all(np.isfinite(a)) for a in [x,dx,y]): raise TrajectoryFailure('Official ode4u returned incomplete/nonfinite arrays.')
    balance=x.sum(axis=0)-x[:,0].sum()-a2.cumulative_trapezoid(y[a2.Y['external_net']],t)
    inflow=np.trapezoid(y[a2.Y['P_total']]+u[a2.U['Q_ext']]+y[a2.Y['Q_gx_in']],t)
    scale=inflow if inflow>0 else x[:,0].sum()
    relative=float(np.max(np.abs(balance))/scale) if scale>0 else (0. if not np.any(balance) else float('inf'))
    if relative>=params.balance_rtol: raise TrajectoryFailure('Water-balance numerical check failed.')
    residuals,audit=constraint_audit(t,x,y,u,theta,params,config)
    k=int(np.argmax(y[a2.Y['H_r']]))
    metrics=dict(peak_stage_m=float(y[a2.Y['H_r'],k]),peak_time_s=float(t[k]),peak_outlet_m3s=float(max(y[a2.Y['Q_out']])),
                 outlet_volume_m3=float(np.trapezoid(y[a2.Y['Q_out']],t)),end_storage_m3=x[:,-1].tolist(),
                 balance_max_m3=float(max(abs(balance))),balance_relative=relative,
                 peak_near_horizon=bool(t[k]>=.95*t[-1]),river_rising_at_end=bool(y[a2.Y['H_r'],-1]>y[a2.Y['H_r'],-2]),
                 flood_duration_h=None if config['reporting']['flood_threshold_m'] is None else a2.flood_duration(t,y[a2.Y['H_r']],config['reporting']['flood_threshold_m'])/3600,
                 regime_rhs_counts={str(key):value for key,value in regimes.items()})
    return dict(time=t,x=x,dx=dx,y=y,u=u,residuals=residuals,audit=audit,metrics=metrics,theta=np.array(theta),pair=list(pair),balance=balance)


def is_feasible(result,config):
    return bool(all(r['pass_at_tolerance'] for r in result['audit']))


class Evaluator:
    def __init__(self,config,pair,forcing,global_budget):
        self.config=copy.deepcopy(config); self.pair=pair; self.forcing=forcing
        self.records=[]; self.cache=OrderedDict(); self.best=None; self.budget=global_budget
    def __call__(self,theta,constants=None):
        from B.gate_hydraulics import UnsupportedGateRegime
        key=(config_hash(self.config),tuple(self.pair),np.asarray(theta,dtype=float).tobytes())
        cache_hit=key in self.cache
        if cache_hit: objective,residuals,result,error=self.cache[key]
        else:
            if self.budget['evaluations']>=self.config['numerics']['total_ode_budget']-3 or clock.monotonic()-self.budget['start']>=self.config['numerics']['total_seconds']:
                raise RunBudgetExceeded('Shared ODE/time budget reached; last in-flight evaluation may exceed wall time.')
            self.budget['evaluations']+=1
            try:
                result=simulate(theta,self.pair,self.config,forcing=self.forcing)
                objective=(result['metrics']['peak_stage_m']-self.config['numerics']['objective_origin_m'])/self.config['numerics']['height_scale_m']
                residuals=result['residuals']; error=None
            except (FloatingPointError,TrajectoryFailure,UnsupportedGateRegime) as exc:
                objective=self.config['numerics']['failure_objective']; residuals=np.full(CONSTRAINT_COUNT,self.config['numerics']['failure_residual'])
                result=None; error=type(exc).__name__+': '+str(exc)
            self.cache[key]=(objective,residuals,result,error)
            if len(self.cache)>8: self.cache.popitem(last=False)
        if result is not None and is_feasible(result,self.config) and (self.best is None or result['metrics']['peak_stage_m']<self.best['metrics']['peak_stage_m']): self.best=result
        self.records.append(dict(evaluation=len(self.records)+1,category_pair='/'.join(self.pair),cached=cache_hit,
            raw_peak_stage_m=None if result is None else result['metrics']['peak_stage_m'],conditioned_objective=float(objective),
            max_scaled_violation=float(max(0.,max(residuals))),feasible=result is not None and is_feasible(result,self.config),
            best_feasible_peak_m=None if self.best is None else self.best['metrics']['peak_stage_m'],failure=error or '',hydraulically_invalid=bool(error and ('HydraulicApplicabilityFailure' in error or 'UnsupportedGateRegime' in error)),
            decision_json=json.dumps(np.asarray(theta).tolist())))
        return float(objective),np.array(residuals,copy=True)


def export_gate_rating(reference,config):
    """Approved-parameter hydraulic diagnostic before any optimization search."""
    import matplotlib.pyplot as plt
    plotdir=ROOT/'plots'; plotdir.mkdir(parents=True,exist_ok=True)
    from B.gate_hydraulics import sluice_gate_flow,UnsupportedGateRegime
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for i,(ax,gate) in enumerate(zip(axes,gate_parameters(config))):
        tail=float(np.median(reference['y'][a2.Y[f'H_w{i+1}']]))
        levels=np.linspace(min(reference['y'][a2.Y[f'H_d{i+1}']]),max(reference['y'][a2.Y[f'H_d{i+1}']]),150)
        for opening in [0,.25,.5,1]:
            flows=[]
            for level in levels:
                try: flows.append(sluice_gate_flow(level,tail,opening,gate).discharge)
                except UnsupportedGateRegime: flows.append(np.nan)
            ax.plot(levels-gate.gate_sill_z,flows,label=f'g={opening}')
        ax.set(xlabel='Pool depth above explicit sill (m)',ylabel='Gate discharge (m³/s)',title=f'Gate {i+1}; tailwater={tail:.3f} m\nGaps mark unsupported regimes')
        ax.legend()
    fig.tight_layout();fig.savefig(plotdir/'B2_0_gate_rating.png',dpi=160);plt.close(fig)


def export_comparison(reference,best,config,records):
    import matplotlib.pyplot as plt
    plotdir=ROOT/'plots'; plotdir.mkdir(parents=True,exist_ok=True)
    fig,ax=plt.subplots(figsize=(9,4))
    for run,label in [(reference,'Reference'),(best,'Best feasible found')]:
        ax.plot(run['time']/3600,run['y'][a2.Y['H_r']],label=label)
        m=run['metrics'];ax.plot(m['peak_time_s']/3600,m['peak_stage_m'],'o')
        ax.annotate(f"{label}: {m['peak_stage_m']:.6f} m at {m['peak_time_s']/3600:.2f} h",
                    (m['peak_time_s']/3600,m['peak_stage_m']),xytext=(5,10 if label=='Reference' else -20),textcoords='offset points',fontsize=8)
    reduction=reference['metrics']['peak_stage_m']-best['metrics']['peak_stage_m']
    ax.set(xlabel='Elapsed time (h)',ylabel='River stage (m, configured datum)',title=f'Peak reduction {reduction:.6f} m; refined-grid feasible, illustrative unless locally supported')
    ax.legend();fig.tight_layout();fig.savefig(plotdir/'B2_1_river_stage.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(9,6),sharex=True)
    for i,ax in enumerate(axes):
        for run,label in [(reference,'Reference'),(best,'Best feasible found')]:
            ax.plot(run['time']/3600,run['u'][i],label=label)
        knots=config['controls'][f'gate{i+1}_hours']
        values=decode(best['theta'],config)[i]
        ax.scatter(knots,values,s=14)
        ax.set(ylabel=f'g{i+1} (fraction)',ylim=(-.02,1.02));ax.axhline(0,color='gray',lw=.7);ax.axhline(1,color='gray',lw=.7)
        ax.legend()
    axes[-1].set_xlabel('Elapsed time (h)');fig.tight_layout();fig.savefig(plotdir/'B2_2_controls.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(10,7))
    limits=config['safety']['dam_m']+config['safety']['wetland_m']
    for i,(name,ax) in enumerate(zip(['H_d1','H_d2','H_w1','H_w2'],axes.flat)):
        for run,label in [(reference,'Reference'),(best,'Best feasible found')]: ax.plot(run['time']/3600,run['y'][a2.Y[name]],label=label)
        ax.axhline(limits[i],color='red',ls='--',label='Configured safety level')
        ax.set(title=name,xlabel='Elapsed time (h)',ylabel='Stage (m)');ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(plotdir/'B2_3_upstream_safety.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(9,6),sharex=True)
    for pair in sorted(set(r['category_pair'] for r in records)):
        subset=[r for r in records if r['category_pair']==pair]; ind=[r['global_evaluation'] for r in subset]
        axes[0].plot(ind,[np.nan if r['raw_peak_stage_m'] is None else r['raw_peak_stage_m'] for r in subset],'.',label=pair)
        axes[0].plot(ind,[np.nan if r['best_feasible_peak_m'] is None else r['best_feasible_peak_m'] for r in subset],lw=1)
        axes[1].plot(ind,[r['max_scaled_violation'] for r in subset],'.')
    axes[0].set_ylabel('Raw peak stage (m)');axes[0].legend(title='Category pair; line=best feasible')
    axes[1].axhline(config['numerics']['feasibility_tol'],color='red',ls='--',label='Feasibility tolerance')
    axes[1].set(xlabel='Objective/constraint callback evaluations (including cache hits)',ylabel='Max scaled violation');axes[1].set_yscale('symlog',linthresh=1e-6);axes[1].legend()
    fig.tight_layout();fig.savefig(plotdir/'B2_4_optimization_history.png',dpi=160);plt.close(fig)
    rows=[]
    for i in range(2):
        rows.append(dict(decision=f'tau_{i+1}',time_h='',reference=reference['pair'][i],optimized=best['pair'][i]))
        rows.append(dict(decision=f'rho_{i+1}',time_h='',reference=reference['theta'][-2+i],optimized=best['theta'][-2+i]))
        for hour,old,new in zip(config['controls'][f'gate{i+1}_hours'],decode(reference['theta'],config)[i],decode(best['theta'],config)[i]):
            rows.append(dict(decision=f'g{i+1}',time_h=hour,reference=old,optimized=new))
    write_table('decision_values.csv',rows)
    write_table('performance.csv',[dict(metric=k,reference=json.dumps(v),best_feasible=json.dumps(best['metrics'][k])) for k,v in reference['metrics'].items()]+[dict(metric='absolute_peak_reduction_m',reference='0',best_feasible=reduction)])
    write_table('constraint_audit.csv',[dict(run=label,**row) for run,label in [(reference,'reference_refined'),(best,'best_refined')] for row in run['audit']])


def run_project(config):
    """First exploratory case: baseline check, explicit reference search, official OPT."""
    from B.official_opt import run_opt,load_official
    from B.gate_hydraulics import UnsupportedGateRegime
    from B.exploratory_reporting import export_case,save_run
    gaps=preflight(config)
    if gaps: raise PreflightError('Project optimization blocked; run --preflight.')
    out=ROOT/'results';out.mkdir(parents=True,exist_ok=True)
    save_json(out/'run_config.json',config)
    official=load_official()
    save_json(out/'data_manifest.json',dict(configuration_sha256=config_hash(config),forcing=config['forcing'],
        geography=config['geography'],parameter_overrides=config['parameter_overrides'],
        official_revision=official['revision'],official_version=official['version'],
        wetland_to_gate_sill_clearance_m=(np.asarray(config['gates']['gate_sill_z'])-np.asarray(config['safety']['wetland_m'])),
        tailwater_interpretation='Gate tailwater equals associated wetland level; safety-feasible wetlands lie below both sills, excluding wetland-induced submergence.'))
    forcing=fixed_forcing(config,base_parameters(config)); original=reference_decision(config)
    budget=dict(evaluations=0,start=clock.monotonic()); reference=None; ref_error=''
    pair_ref=config['reference']['categories']; label='Original schedule under B.2 hydraulics and event forcing'
    def integrate(theta,pair,refine=1):
        if budget['evaluations']>=config['numerics']['total_ode_budget'] or clock.monotonic()-budget['start']>=config['numerics']['total_seconds']:
            raise RunBudgetExceeded('Shared integration/time budget exhausted.')
        budget['evaluations']+=1
        return simulate(theta,pair,config,refine=refine,forcing=forcing)
    start=clock.monotonic(); replacement=False; trials=[]
    try:
        reference=integrate(original,pair_ref)
        original_check=dict(status='valid',label=label,metrics=reference['metrics'],audit=reference['audit'])
    except HydraulicApplicabilityFailure as exc:
        ref_error=str(exc)
        original_check=dict(status='hydraulically_invalid',label=label,diagnostic=exc.diagnostic,physical_peak_m=None)
    except (FloatingPointError,TrajectoryFailure,UnsupportedGateRegime) as exc:
        ref_error=str(exc);original_check=dict(status='simulation_failed',label=label,failure=ref_error,physical_peak_m=None)
    save_json(out/'original_baseline_check.json',original_check)
    save_json(out/'pilot.json',dict(evaluations=1,seconds=clock.monotonic()-start,
        reference_feasible=False if reference is None else is_feasible(reference,config),failure=ref_error))
    if original_check['status']=='hydraulically_invalid':
        # New open-loop trials, never a state-dependent clip or altered geometry.
        for factor in config['reference_search']['uniform_scales']:
            theta=original.copy();theta[:-2]*=factor
            try:
                trial=integrate(theta,pair_ref)
                row=dict(scale=factor,status='hydraulically_and_numerically_valid',safety_and_other_constraints_feasible=is_feasible(trial,config),
                    physical_peak_m=trial['metrics']['peak_stage_m'],failure='',diagnostic=None)
                if reference is None or factor>selected_factor:
                    reference=trial;selected_factor=factor
            except (FloatingPointError,TrajectoryFailure,UnsupportedGateRegime) as exc:
                row=dict(scale=factor,status='rejected',safety_and_other_constraints_feasible=False,
                    physical_peak_m=None,failure=str(exc),diagnostic=getattr(exc,'diagnostic',None))
            trials.append(row)
            print('Reference search:',row['scale'],row['status'],flush=True)
        replacement=reference is not None
        if replacement:label=f'Replacement reference ({selected_factor:g} × original openings; feasibility search)'
    save_json(out/'reference_search.json',dict(method='Declared uniform open-loop scales, unchanged knot times/covers/classes; choose largest tested hydraulically valid scale. Not an optimized original baseline.',trials=trials,
        selected_scale=selected_factor if replacement else None,label=label))
    if reference is None:
        save_json(out/'run_status.json',dict(status='reference_simulation_failed',reason=ref_error,reported_feasible=False,ode_evaluations=budget['evaluations']))
        return
    theta=reference['theta'];save_run(reference,'reference_coarse')
    save_json(out/'reference_assessment.json',dict(label=label,replacement=replacement,metrics=reference['metrics'],
        constraints=reference['audit'],feasible=is_feasible(reference,config),original=original_check))
    print('Reference accepted:',label,'physical peak',reference['metrics']['peak_stage_m'],'feasible',is_feasible(reference,config),flush=True)
    export_gate_rating(reference,config)
    export_case(reference,None,config,[],label,'No optimized candidate yet')
    records=[];summaries=[];best=reference if is_feasible(reference,config) else None
    pairs=list(itertools.product(*config['categories']['admissible']))
    stop=False
    for pair in pairs:
        evaluator=Evaluator(config,pair,forcing,budget)
        try:evaluator(theta) # preserve exact physical reference; NMS may move its interior start
        except RunBudgetExceeded:break
        for phase in ('full_box','reference_neighborhood'):
            settings=dict(config['optimizer']);settings['max_evals']=config['search'][phase+'_evaluations']
            upper=np.ones_like(theta)
            if phase=='reference_neighborhood':
                # Explicit NUMERICAL search subdomain, not a changed physical constraint.
                # It avoids losing a small-opening feasible seed to NMS's interior-start clamp.
                upper[:-2]=np.minimum(1.,np.maximum(config['search']['minimum_local_upper'],
                    config['search']['local_gate_upper_multiplier']*theta[:-2]))
            outcome={};status='not_started'
            try:
                outcome=run_opt(evaluator,theta,np.zeros_like(theta),upper,settings)
                status=outcome['status']
            except RunBudgetExceeded as exc:status='shared_budget: '+str(exc);stop=True
            slug='_'.join(pair)+'_'+phase
            save_json(out/(slug+'_optimizer.json'),dict(category_pair=list(pair),phase=phase,settings=settings,
                lower=np.zeros_like(theta),upper=upper,requested_start=theta,termination=status,
                objective_fields='solver.raw_objective is conditioned (Hpeak-97m)/1m for valid candidates, or a numerical failure sentinel; solver.penalized_objective includes constraint penalties. Physical peaks are in optimization_evaluations.csv.',solver=outcome))
            summaries.append(dict(category_pair='/'.join(pair),phase=phase,termination=status,
                callback_budget=settings['max_evals'],actual_callbacks=len(outcome.get('evaluations',[])),
                official_returned=outcome.get('official_returned',False),iterations=outcome.get('n_iterations'),
                converged=outcome.get('converged',False),incumbent_feasible=evaluator.best is not None,
                phase_feasible_callbacks=sum(r['feasible'] for r in outcome.get('evaluations',[])),
                phase_invalid_callbacks=sum(bool(np.all(np.asarray(r['constraints'])==config['numerics']['failure_residual'])) for r in outcome.get('evaluations',[])),
                incumbent_best_physical_peak_m=None if evaluator.best is None else evaluator.best['metrics']['peak_stage_m']))
            print('OPT:',slug,status,'best physical peak',summaries[-1]['incumbent_best_physical_peak_m'],flush=True)
            if stop:break
        for row in evaluator.records:records.append(dict(global_evaluation=len(records)+1,**row))
        if evaluator.best is not None and (best is None or evaluator.best['metrics']['peak_stage_m']<best['metrics']['peak_stage_m']):best=evaluator.best
        write_table('optimization_evaluations.csv',records);write_table('category_run_summary.csv',summaries)
        if stop:break
    if best is None:
        export_case(reference,None,config,records,label,'No feasible solution found')
        save_json(out/'run_status.json',dict(status='no_feasible_solution_found',reference_label=label,
            reference_metrics=reference['metrics'],reference_feasible=is_feasible(reference,config),
            ode_evaluations=budget['evaluations'],reported_feasible=False,note='Finite exploratory search; does not prove infeasibility.'))
        return
    save_run(best,'best_coarse')
    try:
        fresh=integrate(best['theta'],best['pair'])
        fine=integrate(best['theta'],best['pair'],config['numerics']['refine_factor'])
        fine_ref=integrate(theta,pair_ref,config['numerics']['refine_factor'])
        if not is_feasible(fresh,config) or not is_feasible(fine,config):raise TrajectoryFailure('Fresh/refined candidate violates constraints.')
        sensitivity=abs(fine['metrics']['peak_stage_m']-fresh['metrics']['peak_stage_m'])
        reference_sensitivity=abs(fine_ref['metrics']['peak_stage_m']-reference['metrics']['peak_stage_m'])
        if max(sensitivity,reference_sensitivity)>config['numerics']['peak_grid_tolerance_m']:raise TrajectoryFailure('Reference or candidate peak-stage sensitivity exceeds declared grid check tolerance.')
    except (FloatingPointError,TrajectoryFailure,UnsupportedGateRegime,RunBudgetExceeded) as exc:
        export_case(reference,best,config,records,label,'Coarse candidate; refinement NOT verified')
        save_json(out/'run_status.json',dict(status='refined_verification_failed',reason=str(exc),reported_feasible=False,ode_evaluations=budget['evaluations']))
        return
    save_run(fine_ref,'reference');save_run(fine,'best_feasible')
    export_case(fine_ref,fine,config,records,label,'Best feasible found (bounded official OPT search)')
    save_json(out/'run_status.json',dict(status='best_feasible_solution_found_under_stated_configuration',reported_feasible=True,
        reference_label=label,reference_replaced=replacement,original_baseline=original_check,
        peak_grid_difference_m=sensitivity,reference_peak_grid_difference_m=reference_sensitivity,
        ode_evaluations=budget['evaluations'],seconds=clock.monotonic()-budget['start'],perfect_foresight=True,
        continuous_time_guarantee=False,reference_feasible=is_feasible(fine_ref,config),
        reference_metrics=fine_ref['metrics'],best_metrics=fine['metrics'],best_pair=best['pair'],
        post_horizon_check_performed=False,termination=summaries,
        note='First exploratory regional-weather-driven hypothetical installation. Neither local calibration nor safety certification. Physical peak values are never penalized scores.'))


def plumbing_check():
    """Time one unchanged A.2 ODE evaluation; do not create or overwrite A.2 outputs."""
    params=a2.c();u=a2.make_inputs(params);a2.validate(params,u);ode=a2.load_ode4u();start=clock.monotonic()
    t,x,dx,y=ode(a2.system,params.time,params.x_init,u,params)
    summary,_=a2.baseline_diagnostics(params,t,x,dx,u,y,ode)
    elapsed=clock.monotonic()-start
    summary.update(label='Original synthetic A.2 plumbing check only; NOT revised B.2 baseline or optimization',seconds=elapsed,
                   local_observations_used=False,b2_safe_thresholds_assumed=False)
    save_json(ROOT/'results/a2_plumbing_check.json',summary)
    write_table('execution_status.csv',[dict(check='A.2 official ode4u plumbing',status='executed',value=summary['Q_peak_m3s'],units='m3/s',interpretation=summary['label'])])
    print(json.dumps(summary,indent=2))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=CONFIG)
    action=parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--preflight',action='store_true');action.add_argument('--plumbing-check',action='store_true');action.add_argument('--run',action='store_true')
    args=parser.parse_args()
    if args.plumbing_check: plumbing_check();return
    config=load_config(args.config)
    if args.preflight:
        gaps=preflight(config)
        write_table('approval_gaps.csv',gaps or [dict(id='ready',topic='Configuration',status='ready',missing='',required_action='Run bounded pilot and official OPT')])
        save_json(ROOT/'results/preflight.json',dict(status='blocked' if gaps else 'ready',gap_count=len(gaps),configuration_sha256=config_hash(config),local_observations_used=False))
        for row in gaps: print(row['id']+' '+row['topic']+': '+row['missing'])
        if gaps: raise SystemExit(2)
        print('Preflight ready.');return
    try: run_project(config)
    except PreflightError as exc: parser.exit(2,str(exc)+'\n')

if __name__=='__main__': main()
