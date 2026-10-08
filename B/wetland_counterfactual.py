"""Fixed-design six-state bypass comparison; no optimization or historical edits.

Run: python3 -B -m B.wetland_counterfactual --run [--report]
The former wetland areas remain rainfall footprints, not dynamic compartments.
"""
from __future__ import annotations

import argparse
from collections import Counter
import copy
import hashlib
import json
from pathlib import Path
from B.project_paths import matches_historical_source
import numpy as np
from A import wetland_flood_model as a2
from B import optimize_wetland as old
from B import no_spillway_model as wet
from B.gate_hydraulics import sluice_gate_flow, UnsupportedGateRegime
from B.official_opt import load_official
from B.no_spillway_study import load_trajectory
from B.wetland_only_study import load_study

ROOT = old.ROOT
OUT = ROOT/'results/wetland_counterfactual'
SAVED = ROOT/'results/wetland_only'
STATE_NAMES = ['S_u', 'S_m', 'S_d1', 'S_d2', 'S_g', 'S_r']
SHARED = [a2.STATE_NAMES.index(name) for name in STATE_NAMES]
OUTPUT_NAMES = ['H_d1', 'H_d2', 'H_g', 'H_r', 'Q_G1', 'Q_G2', 'Q_inf',
                'Q_mg', 'Q_gr', 'Q_gx', 'Q_out', 'E_tot', 'S_tot', 'E_star',
                'Q_gx_in', 'Q_gx_out', 'Q_ud1', 'Q_ud2', 'Q_ur', 'E_m',
                'E_d1', 'E_d2', 'E_r', 'P_total', 'P_former_wetland', 'external_net']
X = {name: i for i, name in enumerate(STATE_NAMES)}
Y = {name: i for i, name in enumerate(OUTPUT_NAMES)}
LABELS = {'bypass': 'No-wetland dynamic bypass', 'reference': 'Reference wetlands',
          'best': 'Best-found wetland design',
          'reference_zero': 'Reference wetlands; zero initial wetland storage',
          'best_zero': 'Best-found wetland design; zero initial wetland storage'}


def sha(path):
    return hashlib.sha256(old.resolve_path(path).read_bytes()).hexdigest()


def inputs(config, gates, times, forcing):
    """Same gate interpolation and SAME 120-second weather knots at every step."""
    gates = np.asarray(gates, dtype=float)
    n1 = len(config['controls']['gate1_hours'])
    n2 = len(config['controls']['gate2_hours'])
    if gates.shape != (n1+n2,) or not np.all(np.isfinite(gates)):
        raise ValueError('Expected the full fixed gate schedule, with no covers.')
    if np.any((gates < 0) | (gates > 1)):
        raise ValueError('Gate openings outside [0,1]; not clipped.')
    source_t, source_u = forcing
    return np.vstack([np.interp(times, np.array(config['controls']['gate1_hours'])*3600, gates[:n1]),
                      np.interp(times, np.array(config['controls']['gate2_hours'])*3600, gates[n1:]),
                      *[np.interp(times, source_t, row) for row in source_u]])


def make_rhs(gates, diagnostics=None):
    """Genuine six-state equations; no eight-state RHS, wetland law or dummy state."""
    def rhs(t, x, u, p):
        x = np.asarray(x)
        if x.shape != (6,) or not np.all(np.isfinite(x)) or np.any(x < 0) or x[X['S_m']] > p.C_m:
            raise FloatingPointError(f'Inadmissible six-state storage at t={t:g} s: {x}')
        Su, Sm, Sd1, Sd2, Sg, Sr = x
        g1, g2, P, Ta, Rs, wind, Qext, hd, hg, B = u
        Sd = np.array([Sd1, Sd2])
        Hd = a2.stage(Sd, p.A_d, p.z_d)
        Hg = a2.stage(Sg, p.S_y*p.A_g, p.z_g)
        Hr = a2.stage(Sr, p.A_r, p.z_r)
        # No wetland-bed check: no wetland compartment. Retain unchanged shared
        # groundwater ceiling, which numerically equals the former lowest bed.
        for name, value, bound in [('H_g', Hg, p.H_g_valid_max),
                                   ('H_g', Hg, p.soil_bottom), ('H_r', Hr, p.H_r_valid_max)]:
            if value >= bound:
                raise wet.TrajectoryFailure(f'Applicability domain exceeded at t={t:g}: {name}')
        Qud = np.minimum(p.Q_udcap, p.k_ud*Su)
        Qur = min(p.Q_urcap, p.k_ur*Su)
        Qinf = p.A_u*p.K_u*a2.availability(Su, p.V_u)*a2.positive_part(1-Sm/p.C_m)
        Qmg = p.k_mg*a2.positive_part(Sm-p.S_fc)
        QG = np.empty(2)
        for i, opening in enumerate((g1, g2)):
            try:
                flow = sluice_gate_flow(float(Hd[i]), float(Hr), opening, gates[i])
            except UnsupportedGateRegime as exc:
                depth = float(Hd[i])-gates[i].gate_sill_z
                height = gates[i].gate_max_opening*float(opening)
                raise wet.HydraulicApplicabilityFailure(dict(branch=i+1, time_s=float(t),
                    water_depth_m=depth, opening_height_m=height, normalized_opening=float(opening),
                    residual_m=1.25*height-depth, upstream_stage_m=float(Hd[i]), tailwater_stage_m=float(Hr),
                    rule='1.25*a-H <= 0 when open, wet and forward-flow', regime=exc.regime,
                    detection='First encountered invalid official ode4u RHS/RK stage; not a continuous crossing')) from exc
            QG[i] = flow.discharge
            if diagnostics is not None:
                diagnostics['regimes'][(i+1, flow.regime)] += 1
                diagnostics['minimum_tailwater_clearance_m'][i] = min(
                    diagnostics['minimum_tailwater_clearance_m'][i], gates[i].gate_sill_z-Hr)
                if opening > 0 and flow.energy_head > 0 and flow.differential_head > 0:
                    diagnostics['minimum_G5_margin_m'][i] = min(
                        diagnostics['minimum_G5_margin_m'][i], flow.energy_head-1.25*flow.opening_height)
        Qgr = p.C_gr*a2.availability(Sg, p.V_g)*a2.positive_part(Hg-Hr)
        Qgx_in = p.C_gx*a2.positive_part(hg-Hg)
        Qgx_out = p.C_gx*a2.availability(Sg, p.V_g)*a2.positive_part(Hg-hg)
        Qgx = Qgx_in-Qgx_out
        Qout = (1-B)*a2.head_flow(Hr, hd, p.z_ro, p.R_r, p.Q_rcap)
        E_star = a2.evaporative_demand(Ta, Rs, wind, p)
        Em = p.A_u*p.K_ETm*a2.availability(Sm, p.V_m)*E_star
        Ed = p.A_d*p.K_ETd*a2.availability(Sd, p.V_d)*E_star
        Er = p.A_r*p.K_ETr*a2.availability(Sr, p.V_r)*E_star
        Etot = Em+sum(Ed)+Er
        Pformer = P*sum(p.A_w)  # Rain area only; A_r and A_g remain unchanged.
        dx = np.array([P*p.A_u-Qinf-sum(Qud)-Qur, Qinf-Qmg-Em,
                       P*p.A_d[0]+Qud[0]-QG[0]-Ed[0],
                       P*p.A_d[1]+Qud[1]-QG[1]-Ed[1],
                       Qmg+Qgx-Qgr,
                       P*p.A_r+Pformer+Qur+sum(QG)+Qgr+Qext-Qout-Er])
        net = P*p.A_rain+Qext+Qgx-Qout-Etot
        y = np.array([*Hd, Hg, Hr, *QG, Qinf, Qmg, Qgr, Qgx, Qout, Etot, sum(x), E_star,
                      Qgx_in, Qgx_out, *Qud, Qur, Em, *Ed, Er, P*p.A_rain, Pformer, net])
        if diagnostics is not None:
            diagnostics['rhs_evaluations'] += 1
            diagnostics['derivative_balance_max_m3s'] = max(
                diagnostics['derivative_balance_max_m3s'], abs(float(sum(dx)-net)))
        return dx, y
    return rhs


def validate_six(p, u, initial, gates):
    if u.shape != (10, len(p.time)) or not np.all(np.isfinite(u)):
        raise ValueError('Expected complete finite common inputs.')
    for name in ('g1', 'g2', 'B'):
        if np.any((u[a2.U[name]] < 0) | (u[a2.U[name]] > 1)):
            raise ValueError(name+' outside [0,1]')
    for name in ('P', 'R_s', 'v', 'Q_ext'):
        if np.any(u[a2.U[name]] < 0): raise ValueError(name+' must be nonnegative')
    if np.any(u[a2.U['T_a']] <= 0): raise ValueError('Liquid-water assumption requires positive temperature')
    positive = ['A_u','A_d','A_w','A_r','A_g','S_y','C_m','V_u','V_m','V_d','V_g','V_r',
                'Q_udcap','Q_urcap','k_ud','k_ur','k_mg','R_r','Q_rcap','lambda_v','c_p','epsilon','rho_w','R_v','p_0']
    nonnegative = ['K_u','C_gr','C_gx','K_ETm','K_ETd','K_ETr','k_0','k_v']
    for name in positive+nonnegative:
        values = np.asarray(getattr(p, name))
        if not np.all(np.isfinite(values)) or np.any(values <= 0 if name in positive else values < 0):
            raise ValueError('Invalid retained parameter: '+name)
    for name in ('z_d','z_r','z_g','z_ro','H_g_valid_max','soil_bottom','H_r_valid_max'):
        if not np.all(np.isfinite(getattr(p, name))): raise ValueError('Nonfinite '+name)
    if not (0 <= p.S_fc <= p.C_m and 0 < p.S_y <= 1 and 0 <= p.r_H <= 1 and 0 <= p.eta <= 1):
        raise ValueError('Invalid retained parameter domain')
    if any(gate.gate_sill_z < p.z_d[i] for i, gate in enumerate(gates)) or p.z_ro < p.z_r:
        raise ValueError('Retained outlet lies below donor bed')
    make_rhs(gates)(0., initial, u[:, 0], p)


def audit_six(t, x, y, u, p, config):
    n = config['numerics']; tol = n['feasibility_tol']; rows = []
    def add(label, values, scale, units, threshold, kind):
        k = int(np.argmax(values)); raw = float(values[k])
        rows.append(dict(label=label, kind=kind, threshold=float(threshold), worst_time_s=float(t[k]),
                         margin=-raw, units=units, scale=float(scale), scaled_residual=raw/scale,
                         physical_tolerance=tol*scale, strict_pass=bool(raw <= 0), pass_at_tolerance=bool(raw/scale <= tol)))
    for j, name in enumerate(STATE_NAMES):
        add('nonnegative_'+name, -x[j], n['volume_scales_m3'][SHARED[j]], 'm3', 0, 'shared storage')
    add('soil_capacity', x[X['S_m']]-p.C_m, n['volume_scales_m3'][1], 'm3', p.C_m, 'shared capacity')
    for i in range(2):
        bound = config['safety']['dam_m'][i]
        add(f'H_d{i+1}_safe', y[Y[f'H_d{i+1}']]-bound, n['height_scale_m'], 'm', bound, 'unchanged dam safety')
    for label, name, bound in [('groundwater_geometry','H_g',p.H_g_valid_max),
                               ('groundwater_soil','H_g',p.soil_bottom), ('river_geometry','H_r',p.H_r_valid_max)]:
        add(label, y[Y[name]]-(bound-n['strict_margin_m']), n['height_scale_m'], 'm', bound, 'retained applicability margin')
    add('integration_status', np.full(len(t), -1.), 1., '1', 0, 'numerical validity')
    for i in range(2):
        H = y[Y[f'H_d{i+1}']]-config['gates']['gate_sill_z'][i]
        height = config['gates']['gate_max_opening'][i]*u[i]
        active = (u[i]>0)&(H>0)&(y[Y[f'H_d{i+1}']]>y[Y['H_r']])
        add(f'gate{i+1}_head_opening', np.where(active, 1.25*height-H, 0.), n['height_scale_m'], 'm', 0, 'G5 applicability only')
        add(f'gate{i+1}_lower', -u[i], 1., '1', 0, 'prescribed input bound')
        add(f'gate{i+1}_upper', u[i]-1, 1., '1', 1, 'prescribed input bound')
    return rows


def simulate_bypass(config, fixed_gates, *, dt=120, forcing=None):
    p = wet.base_parameters(config); p.time_step = dt
    forcing = forcing if forcing is not None else wet.fixed_forcing(config, wet.base_parameters(config))
    t = p.time; u = inputs(config, fixed_gates, t, forcing)
    original = np.asarray(config['initial_conditions']['x0_m3'], dtype=float)
    if original.shape != (8,): raise ValueError('Expected approved eight-state source initial conditions')
    initial = original[SHARED].copy()
    gates = wet.gate_parameters(config)
    validate_six(p, u, initial, gates)
    diag = dict(regimes=Counter(), minimum_tailwater_clearance_m=[float('inf')]*2,
                minimum_G5_margin_m=[float('inf')]*2, rhs_evaluations=0, derivative_balance_max_m3s=0.)
    ode = load_official()['ode4u']
    t, x, dx, y = ode(make_rhs(gates, diag), t, initial, u, p)
    if x.shape != (6,len(t)) or not all(np.all(np.isfinite(a)) for a in (x,dx,y)):
        raise wet.TrajectoryFailure('Official ode4u returned incomplete six-state arrays')
    balance = x.sum(axis=0)-x[:,0].sum()-a2.cumulative_trapezoid(y[Y['external_net']], t)
    inflow = np.trapezoid(y[Y['P_total']]+u[a2.U['Q_ext']]+y[Y['Q_gx_in']], t)
    relative = float(max(abs(balance))/(inflow if inflow > 0 else sum(initial)))
    if relative >= p.balance_rtol: raise wet.TrajectoryFailure('Integrated six-state water balance failed')
    rows = audit_six(t, x, y, u, p, config)
    k = int(np.argmax(y[Y['H_r']]))
    diag['regimes'] = {str(k):v for k,v in diag['regimes'].items()}
    metrics = dict(peak_stage_m=float(y[Y['H_r'],k]), peak_time_s=float(t[k]),
                   balance_max_m3=float(max(abs(balance))), balance_relative=relative,
                   derivative_balance_max_m3s=diag['derivative_balance_max_m3s'],
                   peak_near_horizon=bool(t[k] >= .95*t[-1]),
                   river_rising_at_end=bool(y[Y['H_r'],-1] > y[Y['H_r'],-2]),
                   integration_stage_diagnostics=diag)
    return dict(time=t, x=x, dx=dx, y=y, u=u, balance=balance, theta=np.asarray(fixed_gates),
                pair=[], audit=rows, residuals=np.array([r['scaled_residual'] for r in rows]), metrics=metrics,
                state_names=STATE_NAMES, output_names=OUTPUT_NAMES, model='six_state_no_wetland_no_dam_spillway')


def save_trajectory(path, run):
    arrays = {k:run[k] for k in ('time','x','dx','y','u','theta','residuals','balance') if k in run}
    metadata = {k:v for k,v in run.items() if k not in arrays}
    np.savez_compressed(path, **arrays, metadata_json=np.asarray(json.dumps(metadata, default=old.serial)))


def column(run, name):
    names = run.get('output_names', a2.OUTPUT_NAMES)
    return run['y'][names.index(name)]


def verify_saved(run, config, expected, forcing, dt):
    """Check exact saved inputs/design/x0 and all accepted-node fluxes/balances.

    Unchanged source hashes and saved execution records support compatible reuse;
    fresh 120s reproduction additionally checks the full integrator trajectory.
    """
    p = wet.base_parameters(config); p.time_step = dt
    if run['pair'] != expected['pair'] or not np.array_equal(run['theta'], expected['theta']):
        raise RuntimeError('Saved fixed design differs from full-precision source')
    if not np.array_equal(run['time'], p.time) or not np.array_equal(run['x'][:,0], config['initial_conditions']['x0_m3']):
        raise RuntimeError('Saved timestep or initial water mismatch')
    u = inputs(config, expected['theta'][:-2], p.time, forcing)
    if not np.array_equal(run['u'], u): raise RuntimeError('Saved weather/gate inputs differ')
    p.rho = expected['theta'][-2:].copy(); p.tau = list(expected['pair'])
    for pair in p.tau:
        endpoints = config['categories']['endpoints'][pair]
        p.R_veg[pair] = np.array(endpoints['R_veg']); p.K_veg[pair] = np.array(endpoints['K_veg'])
    rhs = wet.make_b2_rhs(wet.gate_parameters(config)); max_y = 0.; max_balance = 0.
    for k,t in enumerate(run['time']):
        dx,y = rhs(t,run['x'][:,k],u[:,k],p)
        max_y = max(max_y,float(max(abs(y-run['y'][:,k]))))
        max_balance = max(max_balance,abs(float(sum(dx)-y[a2.Y['external_net']])))
    if max_y > 1e-10: raise RuntimeError('Saved fluxes cannot be reproduced')
    residuals,audit = wet.constraint_audit(run['time'],run['x'],run['y'],u,expected['theta'],p,config)
    if not np.array_equal(residuals,run['residuals']) or audit != run['audit']:
        raise RuntimeError('Saved constraint audit differs')
    balance = run['x'].sum(axis=0)-sum(run['x'][:,0])-a2.cumulative_trapezoid(column(run,'external_net'),run['time'])
    if not np.allclose(balance,run['balance'],atol=1e-9,rtol=0): raise RuntimeError('Saved balance differs')
    return dict(accepted_node_output_max_error=max_y, derivative_balance_max_m3s=max_balance,
                exact_time_inputs_initial_design=True, accepted_nodes_checked=len(run['time']))


def run_comparison():
    if (OUT/'comparison.json').exists(): raise RuntimeError('Completed counterfactual results exist; preserving them')
    OUT.mkdir(parents=True,exist_ok=True)
    config,spec,gates = load_study(); official = load_official()
    forcing = wet.fixed_forcing(config,wet.base_parameters(config))
    reference = load_trajectory(SAVED/'reference.npz'); best = load_trajectory(SAVED/'batch_best.npz')
    if reference['pair'] != spec['reference_categories'] or not np.array_equal(reference['theta'][-2:],spec['reference_covers']):
        raise RuntimeError('Unexpected saved reference design')
    for design in (reference,best):
        if not np.array_equal(design['theta'][:-2],gates): raise RuntimeError('Saved gate schedule mismatch')
    plan = old.load_config(SAVED/'study_plan.json')
    if plan['physical_model_configuration'] != config or plan['specification'] != spec:
        raise RuntimeError('Historical study configuration mismatch')
    before = old.load_config(SAVED/'preservation_before.json')
    # Prior snapshot precedes wetland-only execution. Check the scientific source
    # and forcing files that govern all reused trajectories against that snapshot.
    paths = ['wetland_flood_model.py','optimize_wetland.py','b2/no_spillway_model.py',
             'b2/gate_hydraulics.py','b2/official_opt.py','b2/no_spillway_config.json']
    if 'files' in before: before = before['files']
    for path in paths:
        record = before.get(path)
        digest = record.get('sha256') if isinstance(record,dict) else record
        if not matches_historical_source(path, digest): raise RuntimeError('Historical scientific source changed: '+path)
    manifest = old.load_config(ROOT/config['forcing']['manifest_json'])
    for source in manifest['source_files']:
        if sha(ROOT/source['path']) != source['sha256']: raise RuntimeError('Forcing provenance source changed')
    lock_paths = paths+[str((SAVED/name).relative_to(ROOT)) for name in ['reference.npz','batch_best.npz','study_plan.json','fixed_design_timesteps.json']]
    lock_paths += [config['forcing']['prepared_csv'],config['forcing']['manifest_json'],'b2/wetland_only_study.py','b2/wetland_only_config.json',
                   'work/multivarious/multivarious/ode/ode4u.py']
    lock_paths += [str((SAVED/f'{label}_{dt}s.npz').relative_to(ROOT)) for label in ['fixed_reference','wetland_best'] for dt in [120,60,30]]
    lock = dict(configuration=config, wetland_study=spec, reference_pair=reference['pair'], reference_theta=reference['theta'],
                best_pair=best['pair'],best_theta=best['theta'],fixed_gate_knots=gates,source_sha256={p:sha(ROOT/p) for p in lock_paths},
                official_revision=official['revision'],official_version=official['version'],
                state_names_six=STATE_NAMES,former_wetland_rain_to_river=True,initial_water_transfer=False,
                timesteps_s=[120,60,30],optimizer_run=False)
    old.save_json(OUT/'configuration_and_provenance.json',lock)
    records=[]; reproduction=[]
    # Fresh whole-trajectory checks at the original step; no optimization.
    for name,design in [('reference',reference),('best',best)]:
        check = wet.simulate(design['theta'],design['pair'],config,forcing=forcing)
        saved = load_trajectory(SAVED/f'{"fixed_reference" if name=="reference" else "wetland_best"}_120s.npz')
        errors = {key:float(np.max(np.abs(check[key]-saved[key]))) for key in ['time','x','u','y','balance']}
        if any(value > 1e-9 for value in errors.values()): raise RuntimeError('Saved wetland run reproduction mismatch: '+str(errors))
        reproduction.append(dict(case=name,dt_s=120,errors=errors,completed=True))
        print('Reproduced saved',name,'120s',errors,flush=True)
    for dt in [120,60,30]:
        for name in LABELS:
            design = reference if name.startswith('reference') else best
            row = dict(case=name,label=LABELS[name],dt_s=dt,reused=False,valid=False,feasible=False)
            try:
                if name in ('reference','best'):
                    path = SAVED/f'{"fixed_reference" if name=="reference" else "wetland_best"}_{dt}s.npz'
                    run = load_trajectory(path)
                    row.update(reused=True,reuse_verification=verify_saved(run,config,design,forcing,dt))
                else:
                    path = OUT/f'{name}_{dt}s.npz'
                    if path.exists(): raise RuntimeError('Existing new trajectory would be overwritten: '+str(path))
                    if name == 'bypass': run = simulate_bypass(config,gates,dt=dt,forcing=forcing)
                    else:
                        diagnostic_config=copy.deepcopy(config)
                        diagnostic_config['initial_conditions']['x0_m3'][2:4]=[0.,0.]
                        run=wet.simulate(design['theta'],design['pair'],diagnostic_config,refine=int(120/dt),forcing=forcing)
                    save_trajectory(path,run)
                row.update(valid=True,feasible=all(r['pass_at_tolerance'] for r in run['audit']),
                           status='valid_feasible' if all(r['pass_at_tolerance'] for r in run['audit']) else 'valid_constraint_violation',
                           trajectory=str(path.relative_to(ROOT)),metrics=run['metrics'],audit=run['audit'],
                           state_count=run['x'].shape[0],pair=run['pair'],covers=None if name=='bypass' else run['theta'][-2:],
                           rain_depth_m=float(np.trapezoid(run['u'][a2.U['P']],run['time'])),
                           solar_energy_Wh_m2=float(np.trapezoid(run['u'][a2.U['R_s']],run['time'])/3600))
                if abs(row['rain_depth_m']-manifest['rain_depth_m'])>manifest['rain_mass_tolerance_m'] or abs(row['solar_energy_Wh_m2']-manifest['solar_energy_Wh_m2'])>1e-6:
                    raise RuntimeError('Forcing integral changed at refined step')
            except (wet.TrajectoryFailure,FloatingPointError) as exc:
                row.update(status='invalid_trajectory',failure=str(exc),diagnostic=getattr(exc,'diagnostic',None))
            records.append(row)
            old.save_json(OUT/'progress.json',dict(records=records))
            print(name,dt,row['status'],'peak',row.get('metrics',{}).get('peak_stage_m'),flush=True)
    for path,digest in lock['source_sha256'].items():
        if sha(ROOT/path) != digest: raise RuntimeError('A frozen source changed during this run')
    old.save_json(OUT/'comparison.json',dict(records=records,reproduction_checks=reproduction,
        new_ode_attempts=2+sum(not r['reused'] for r in records),reused_trajectories=sum(r['reused'] for r in records),
        optimizer_runs=0,completion='All requested fixed comparisons attempted; inspect per-run validity and feasibility'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',action='store_true'); parser.add_argument('--report',action='store_true')
    args=parser.parse_args()
    if args.run: run_comparison()
    if args.report:
        from B.wetland_counterfactual_reporting import report
        report()

if __name__=='__main__': main()
