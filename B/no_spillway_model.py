"""Revised B.2 model with no dam-pool spillway and retained wetland overflow.

A.2 and the first exploratory B.2 implementation remain historical, unchanged.
This module reuses their constitutive helpers, fixed forcing and gate law, but
writes the revised eight balances directly. It never calls A.2 system/validate.
The two historical Q_dsp output positions contain literal compatibility zeros;
they are not computed fluxes, parameters, storage caps, or feasibility limits.
"""
from __future__ import annotations

from collections import Counter
import numpy as np
from A import wetland_flood_model as a2
from B import optimize_wetland as previous
from B.gate_hydraulics import sluice_gate_flow, UnsupportedGateRegime

PreflightError = previous.PreflightError
TrajectoryFailure = previous.TrajectoryFailure
HydraulicApplicabilityFailure = previous.HydraulicApplicabilityFailure
CONSTRAINT_COUNT = previous.CONSTRAINT_COUNT
constraint_audit = previous.constraint_audit
is_feasible = previous.is_feasible
fixed_forcing = previous.fixed_forcing
gate_parameters = previous.gate_parameters
decode = previous.decode
reference_decision = previous.reference_decision
REMOVED_DAM_SPILLWAY_PARAMETERS = frozenset({'z_dsp', 'c_dsp', 'C_d'})
MODEL_VARIANT = 'B.2 no dam spillway; wetland overflow retained'


def base_parameters(config):
    """Fresh retained parameters, with old dam-spillway attributes removed.

    A.2's old C_d property consequently has no backing crest attribute. Nothing
    in this module or the reused constraint audit reads that historical property.
    Supplying removed parameters as overrides is an explicit configuration error.
    """
    params = a2.c()
    del params.z_dsp
    del params.c_dsp
    for name, record in config['parameter_overrides'].items():
        if name in REMOVED_DAM_SPILLWAY_PARAMETERS:
            raise PreflightError('Dam-spillway parameter is not part of the revised model: '+name)
        if name in {'gate1_values', 'gate2_values', 'rho', 'tau', 'duration', 'time_step'}:
            raise PreflightError('Decisions/time settings require their dedicated configuration fields: '+name)
        if not hasattr(params, name) or not isinstance(record, dict):
            raise PreflightError('Invalid parameter override: '+name)
        if not record.get('source_id') or not record.get('reason') or record.get('approved') is not True:
            raise PreflightError('Unapproved/unsourced parameter override: '+name)
        old = getattr(params, name)
        setattr(params, name, np.asarray(record['value'], dtype=float) if isinstance(old, np.ndarray) else record['value'])
    params.duration = float(config['numerics']['duration_s'])
    params.time_step = float(config['numerics']['dt_s'])
    return params


def make_b2_rhs(gates, regimes=None):
    """Use Q_G as the only dam-to-wetland transfer; never clip pool storage."""
    def rhs(t, x, u, params):
        if not np.all(np.isfinite(x)) or np.any(x < 0) or x[1] > params.C_m:
            raise FloatingPointError(f'Inadmissible state at t={t:g} s: {x}. Reduce timestep.')
        S_u, S_m, S_w1, S_w2, S_d1, S_d2, S_g, S_r = x
        g1, g2, P, T_a, R_s, wind, Q_ext, h_d, h_g, B = u
        Sw, Sd = np.array([S_w1, S_w2]), np.array([S_d1, S_d2])
        Hw = a2.stage(Sw, params.A_w, params.z_w)
        Hd = a2.stage(Sd, params.A_d, params.z_d)
        Hg = a2.stage(S_g, params.S_y*params.A_g, params.z_g)
        Hr = a2.stage(S_r, params.A_r, params.z_r)
        for name, level, bound in [('H_g', Hg, params.H_g_valid_max),
                                    ('H_g', Hg, params.soil_bottom),
                                    ('H_g', Hg, min(params.z_w)),
                                    ('H_r', Hr, params.H_r_valid_max)]:
            if level >= bound:
                raise TrajectoryFailure('Applicability domain exceeded at t='+str(t)+': '+name)
        Rw, KETw = params.design()

        # Every transfer is in m3/s and occurs with opposite signs at its ends.
        Qud = np.minimum(params.Q_udcap, params.k_ud*S_u)
        Qur = min(params.Q_urcap, params.k_ur*S_u)
        Qinf = params.A_u*params.K_u*a2.availability(S_u, params.V_u)*a2.positive_part(1-S_m/params.C_m)
        Qmg = params.k_mg*a2.positive_part(S_m-params.S_fc)
        QG = np.empty(2)
        for i, opening_fraction in enumerate((g1, g2)):
            upstream, tail = float(Hd[i]), float(Hw[i])
            try:
                flow = sluice_gate_flow(upstream, tail, opening_fraction, gates[i])
            except UnsupportedGateRegime as exc:
                H = upstream-gates[i].gate_sill_z
                opening = gates[i].gate_max_opening*float(opening_fraction)
                raise HydraulicApplicabilityFailure(dict(
                    branch=i+1, time_s=float(t), water_depth_m=H,
                    opening_height_m=opening, normalized_opening=float(opening_fraction),
                    residual_m=1.25*opening-H, upstream_stage_m=upstream,
                    tailwater_stage_m=tail, rule='1.25*a-H <= 0 when open, wet and forward-flow',
                    regime=exc.regime,
                    detection='First encountered invalid ode4u RHS/RK stage; not a continuous-time crossing estimate'
                )) from exc
            QG[i] = flow.discharge
            if regimes is not None:
                regimes[(i+1, flow.regime)] += 1
        Qwr = a2.head_flow(Hw, Hr, params.z_wo, Rw, params.Q_wrcap)
        Qout = (1-B)*a2.head_flow(Hr, h_d, params.z_ro, params.R_r, params.Q_rcap)
        # Wetland -> river overflow is retained; no dam-pool overflow is computed.
        Qwsp = params.c_wsp*a2.positive_part(Hw-np.maximum(params.z_wsp, Hr))**1.5
        Qwg = params.C_wg*a2.availability(Sw, params.V_w)*a2.positive_part(Hw-Hg)
        Qgr = params.C_gr*a2.availability(S_g, params.V_g)*a2.positive_part(Hg-Hr)
        Qgx_in = params.C_gx*a2.positive_part(h_g-Hg)
        Qgx_out = params.C_gx*a2.availability(S_g, params.V_g)*a2.positive_part(Hg-h_g)
        Qgx = Qgx_in-Qgx_out
        E_star = a2.evaporative_demand(T_a, R_s, wind, params)
        Em = params.A_u*params.K_ETm*a2.availability(S_m, params.V_m)*E_star
        Ew = params.A_w*KETw*a2.availability(Sw, params.V_w)*E_star
        Ed = params.A_d*params.K_ETd*a2.availability(Sd, params.V_d)*E_star
        Er = params.A_r*params.K_ETr*a2.availability(S_r, params.V_r)*E_star
        Etot = Em+sum(Ew)+sum(Ed)+Er

        dSu = P*params.A_u-Qinf-Qud[0]-Qud[1]-Qur
        dSm = Qinf-Qmg-Em
        dSw1 = P*params.A_w[0]+QG[0]-Qwr[0]-Qwsp[0]-Qwg[0]-Ew[0]
        dSw2 = P*params.A_w[1]+QG[1]-Qwr[1]-Qwsp[1]-Qwg[1]-Ew[1]
        dSd1 = P*params.A_d[0]+Qud[0]-QG[0]-Ed[0]
        dSd2 = P*params.A_d[1]+Qud[1]-QG[1]-Ed[1]
        dSg = Qmg+Qwg[0]+Qwg[1]+Qgx-Qgr
        dSr = P*params.A_r+Qur+sum(Qwr)+sum(Qwsp)+Qgr+Q_ext-Qout-Er
        dxdt = np.array([dSu, dSm, dSw1, dSw2, dSd1, dSd2, dSg, dSr])
        external_net = P*params.A_rain+Q_ext+Qgx-Qout-Etot
        # Literal zeros preserve historical output-column positions only.
        y = np.array([*Hw, *Hd, Hg, Hr, *QG, *Qwr, 0., 0., *Qwsp, Qinf, Qmg,
                      *Qwg, Qgr, Qgx, Qout, Etot, sum(x), E_star, Qgx_in, Qgx_out,
                      *Qud, Qur, Em, *Ew, *Ed, Er, P*params.A_rain, external_net])
        return dxdt, y
    return rhs


def validate(params, u, gates):
    """Retained input/parameter checks without any dam-spillway dependency."""
    if u.shape != (10, len(params.time)) or not np.all(np.isfinite(u)):
        raise ValueError('Expected finite u with shape (10, len(params.time)).')
    if params.time_step <= 0 or params.duration <= 0 or not np.all(np.diff(params.time) > 0):
        raise ValueError('Time must be strictly increasing.')
    for name in ('g1', 'g2', 'B'):
        if np.any((u[a2.U[name]] < 0) | (u[a2.U[name]] > 1)):
            raise ValueError(name+' must lie in [0,1].')
    for name in ('P', 'R_s', 'v', 'Q_ext'):
        if np.any(u[a2.U[name]] < 0):
            raise ValueError(name+' must be nonnegative.')
    if np.any(u[a2.U['T_a']] <= 0):
        raise ValueError('This liquid-water demonstration requires T_a > 0 degC.')
    positive = ['A_u', 'A_d', 'A_w', 'A_r', 'A_g', 'S_y', 'C_m', 'V_u', 'V_m',
                'V_w', 'V_d', 'V_g', 'V_r', 'Q_udcap', 'Q_urcap', 'k_ud', 'k_ur',
                'k_mg', 'R_bare', 'Q_wrcap', 'R_r', 'Q_rcap', 'c_wsp', 'lambda_v',
                'c_p', 'epsilon', 'rho_w', 'R_v', 'p_0']
    nonnegative = ['K_u', 'C_wg', 'C_gr', 'C_gx', 'K_bare', 'K_ETm', 'K_ETd',
                   'K_ETr', 'k_0', 'k_v']
    for name in positive+nonnegative:
        values = np.asarray(getattr(params, name))
        if not np.all(np.isfinite(values)) or np.any(values <= 0 if name in positive else values < 0):
            raise ValueError('Invalid parameter: '+name)
    for name in ('z_d', 'z_w', 'z_r', 'z_g', 'z_wo', 'z_wsp', 'z_ro',
                 'H_g_valid_max', 'soil_bottom', 'H_r_valid_max', 'H_g_init'):
        if not np.all(np.isfinite(getattr(params, name))):
            raise ValueError('Nonfinite elevation: '+name)
    if not (0 <= params.S_fc <= params.C_m and 0 < params.S_y <= 1 and
            0 <= params.r_H <= 1 and 0 <= params.eta <= 1 and
            np.all((params.rho >= 0) & (params.rho <= 1))):
        raise ValueError('Invalid soil, yield, humidity, radiation, or design domain.')
    for category in set(params.tau):
        if category not in params.R_veg or category not in params.K_veg:
            raise ValueError('Unknown categorical vegetation type.')
        if (not np.all(np.isfinite(params.R_veg[category])) or
                not np.all(np.isfinite(params.K_veg[category])) or
                np.any(params.R_veg[category] <= 0) or np.any(params.K_veg[category] < 0)):
            raise ValueError('Invalid design endpoints.')
    for i, gate in enumerate(gates):
        if gate.gate_sill_z < params.z_d[i]:
            raise ValueError(f'Branch {i+1}: gate sill lies below the donor pool bed.')
        if params.z_wsp[i] <= params.z_wo[i]:
            raise ValueError(f'Branch {i+1}: wetland overflow crest must exceed its drainage sill.')
    if (np.any(params.z_wo < params.z_w) or params.z_ro < params.z_r or
            np.any(params.z_wsp <= params.z_w)):
        raise ValueError('Retained outlet/overflow threshold lies below required donor bed.')
    a2.check_groundwater_levels(params.H_g_init, params, 'Initial condition')
    make_b2_rhs(gates)(0., params.x_init, u[:, 0], params)


def preflight(config):
    """Check the approved revised case without consulting an old spill crest."""
    gaps = []
    checks = [
        ('G1', config['mode'] in {'approved_illustrative', 'approved_local'}, 'Study case'),
        ('G2', config['safety']['approved'] and all(v is not None for v in config['safety']['dam_m']+config['safety']['wetland_m']), 'Four fixed safety levels'),
        ('G3', config['categories']['approved'] and all(config['categories']['admissible']), 'Synthetic response classes'),
        ('G4', config['gates']['approved'], 'Fixed gates'),
        ('G5', config['gates']['limited_regimes_approved'] is True, 'Additional gate applicability constraint'),
        ('G6', config['fixed_model']['approved'] and config['forcing']['approved'] and config['initial_conditions']['approved'], 'Fixed model, forcing and initial state'),
    ]
    for code, approved, topic in checks:
        if not approved:
            gaps.append(dict(id=code, topic=topic, status='blocked', missing='Required case approval/value missing.', required_action='Supply explicit approved configuration.'))
    n = config['numerics']
    if config['optimizer']['feasibility_tol'] != n['feasibility_tol']:
        raise PreflightError('Optimizer and independent audit feasibility tolerances must match.')
    values = list(n['volume_scales_m3'])+[n['height_scale_m'], n['strict_margin_m'], n['feasibility_tol']]
    if not all(np.isfinite(v) for v in values) or n['feasibility_tol'] < 0:
        raise PreflightError('Residual scales and numerical tolerances must be finite and admissible.')
    if len(n['volume_scales_m3']) != 8 or min(n['volume_scales_m3']) <= 0 or n['height_scale_m'] <= 0:
        raise PreflightError('All fixed residual scales must be positive and dimensionally specified.')
    if n['strict_margin_m'] <= n['feasibility_tol']*n['height_scale_m']:
        raise PreflightError('Numerical strict-domain margin must exceed feasibility tolerance in metres.')
    if n['dt_s'] <= 0 or n['refine_factor'] < 2 or int(n['refine_factor']) != n['refine_factor']:
        raise PreflightError('Require positive dt and integer refinement factor >=2.')
    for branch in (1, 2):
        times = np.asarray(config['controls'][f'gate{branch}_hours'])*3600
        if times[0] != 0 or times[-1] != n['duration_s'] or np.any(np.diff(times) <= 0):
            raise PreflightError('Gate knot times must increase and cover the fixed horizon exactly.')
        if not np.allclose(times/n['dt_s'], np.round(times/n['dt_s']), rtol=0, atol=1e-10):
            raise PreflightError('Gate knots must be represented exactly on the ODE grid.')
    decode(reference_decision(config), config)
    if gaps:
        return gaps
    levels = config['safety']['dam_m']+config['safety']['wetland_m']
    if len(levels) != 4 or not np.all(np.isfinite(levels)):
        raise PreflightError('Exactly four finite fixed safety levels required.')
    if len(config['safety']['source_ids']) != 4 or not all(config['safety']['source_ids']):
        raise PreflightError('Each safety level requires its declared source/assumption ID.')
    fields = ['gate_width', 'gate_max_opening', 'gate_sill_z', 'Cd_gate_free', 'Cd_gate_submerged', 'gravity', 'neglect_approach_velocity_head']
    if not all(config['gates']['source_ids'].get(k) for k in fields):
        raise PreflightError('Every gate parameter/assumption requires a source ID.')
    params = base_parameters(config)
    for i, gate in enumerate(gate_parameters(config)):
        if gate.gate_sill_z < params.z_d[i]:
            raise PreflightError('Gate sill must be at or above the donor pool bed.')
    for labels in config['categories']['admissible']:
        if any(label not in config['categories']['endpoints'] for label in labels):
            raise PreflightError('Missing paired synthetic-class endpoints.')
    if config['mode'] == 'approved_local':
        # This revision is not a new local calibration path.
        required = ['polygon_path', 'polygon_sha256', 'outlet', 'crs', 'vertical_datum', 'event_start_utc']
        if not config['geography']['approved'] or any(not config['geography'].get(k) for k in required):
            raise PreflightError('Local run requires approved polygon/outlet/CRS/datum/time metadata.')
        if config['forcing']['kind'] != 'prepared_csv':
            raise PreflightError('A local-data run cannot use synthetic forcing.')
        import hashlib
        if hashlib.sha256((previous.resolve_path(config['geography']['polygon_path'])).read_bytes()).hexdigest() != config['geography']['polygon_sha256']:
            raise PreflightError('Approved polygon checksum mismatch.')
    fixed_forcing(config, params)
    return []


def simulate(theta, pair, config, *, refine=1, forcing=None):
    """Official ode4u, unchanged forcing/x0; only the two spillway links removed."""
    params = base_parameters(config)
    params.time_step /= refine
    g1, g2, covers = decode(theta, config)
    params.rho = covers.copy()
    params.tau = list(pair)
    endpoints = config['categories']['endpoints']
    for label in pair:
        params.R_veg[label] = np.array(endpoints[label]['R_veg'], dtype=float)
        params.K_veg[label] = np.array(endpoints[label]['K_veg'], dtype=float)
    gates = gate_parameters(config)
    t = params.time
    source_t, source_u = forcing if forcing is not None else fixed_forcing(config, base_parameters(config))
    u = np.vstack([np.interp(t, np.array(config['controls']['gate1_hours'])*3600, g1),
                   np.interp(t, np.array(config['controls']['gate2_hours'])*3600, g2),
                   *[np.interp(t, source_t, row) for row in source_u]])
    initial = np.asarray(config['initial_conditions']['x0_m3'], dtype=float)
    if initial.shape != (8,) or not np.all(np.isfinite(initial)):
        raise PreflightError('Eight fixed finite initial storages required.')
    params.u_depth_init = initial[0]/params.A_u
    params.soil_fraction_init = initial[1]/params.C_m
    params.w_depth_init = initial[2:4]/params.A_w
    params.d_depth_init = initial[4:6]/params.A_d
    params.H_g_init = a2.stage(initial[6], params.S_y*params.A_g, params.z_g)
    params.r_depth_init = initial[7]/params.A_r
    validate(params, u, gates)
    if np.any(initial < 0) or initial[1] > params.C_m:
        raise PreflightError('Fixed initial state outside storage domain.')
    regimes = Counter()
    ode = a2.load_ode4u()
    t, x, dx, y = ode(make_b2_rhs(gates, regimes), t, initial.copy(), u, params)
    if not all(np.all(np.isfinite(arr)) for arr in [x, dx, y]):
        raise TrajectoryFailure('Official ode4u returned incomplete/nonfinite arrays.')
    balance = x.sum(axis=0)-x[:, 0].sum()-a2.cumulative_trapezoid(y[a2.Y['external_net']], t)
    inflow = np.trapezoid(y[a2.Y['P_total']]+u[a2.U['Q_ext']]+y[a2.Y['Q_gx_in']], t)
    scale = inflow if inflow > 0 else x[:, 0].sum()
    relative = float(np.max(np.abs(balance))/scale) if scale > 0 else (0. if not np.any(balance) else float('inf'))
    if relative >= params.balance_rtol:
        raise TrajectoryFailure('Water-balance numerical check failed.')
    residuals, audit = constraint_audit(t, x, y, u, theta, params, config)
    k = int(np.argmax(y[a2.Y['H_r']]))
    metrics = dict(peak_stage_m=float(y[a2.Y['H_r'], k]), peak_time_s=float(t[k]),
        peak_outlet_m3s=float(max(y[a2.Y['Q_out']])), outlet_volume_m3=float(np.trapezoid(y[a2.Y['Q_out']], t)),
        end_storage_m3=x[:, -1].tolist(), balance_max_m3=float(max(abs(balance))), balance_relative=relative,
        derivative_balance_max_m3s=float(np.max(np.abs(dx.sum(axis=0)-y[a2.Y['external_net']]))),
        peak_near_horizon=bool(t[k] >= .95*t[-1]),
        river_rising_at_end=bool(y[a2.Y['H_r'], -1] > y[a2.Y['H_r'], -2]),
        flood_duration_h=None if config['reporting']['flood_threshold_m'] is None else a2.flood_duration(t, y[a2.Y['H_r']], config['reporting']['flood_threshold_m'])/3600,
        regime_rhs_counts={str(key): value for key, value in regimes.items()}, model_variant=MODEL_VARIANT,
        historical_output_slots='Q_dsp1/Q_dsp2 are literal compatibility zeros; no dam-spillway law exists in this RHS.')
    return dict(time=t, x=x, dx=dx, y=y, u=u, residuals=residuals, audit=audit,
                metrics=metrics, theta=np.array(theta), pair=list(pair), balance=balance)
