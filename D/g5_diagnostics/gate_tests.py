#!/usr/bin/env python3
"""Read-only Stage B hydraulic regression checks; no low-head closure is invented.

Run with python -B g5_gate_tests.py --project-root PATH --output-dir PATH.
Only the requested output directory is written. RHS identity checks are not
integration-convergence evidence, safety certification, or complete trajectories.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from dataclasses import replace

sys.dont_write_bytecode = True
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/g5-test-matplotlib')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root', type=Path, default=Path('/Users/zhangjiaxuan/Documents/Codex/Modeling'))
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    root = args.project_root.resolve()
    sys.path.insert(0, str(root))
    import numpy as np
    from A import wetland_flood_model as a2
    from B import no_spillway_model as wet
    from B import wetland_counterfactual as dry
    from B.gate_hydraulics import sluice_gate_flow, UnsupportedGateRegime, GateInputError
    from B.wetland_only_study import load_study
    from B.no_spillway_study import load_trajectory

    config, spec, fixed_gates = load_study()
    design = load_trajectory(root / 'B/results/wetland_only/wetland_best_120s.npz')
    p = wet.base_parameters(config)
    p.rho = np.asarray(design['theta'][-2:]).copy()
    p.tau = list(design['pair'])
    for label in p.tau:
        p.R_veg[label] = np.asarray(config['categories']['endpoints'][label]['R_veg'])
        p.K_veg[label] = np.asarray(config['categories']['endpoints'][label]['K_veg'])
    gates = wet.gate_parameters(config)
    x0 = np.asarray(config['initial_conditions']['x0_m3'], dtype=float)
    historical = load_trajectory(root / 'B/results/no_spillway_revision/reference.npz')
    u0 = historical['u'][:, 0].copy()
    records, coverage = [], []

    def check(name, group, callback):
        try:
            detail = callback() or {}
            records.append(dict(name=name, group=group, regression_pass=True, **detail))
        except Exception as exc:
            records.append(dict(name=name, group=group, regression_pass=False,
                                error_type=type(exc).__name__, error=str(exc)))

    def close(a, b, atol=1e-12, rtol=1e-10):
        assert np.allclose(a, b, atol=atol, rtol=rtol), (a, b)

    # SI geometric values are unchanged; zero sill merely changes vertical datum.
    # A binary-exact 0.125 m opening avoids a rounding ambiguity exactly at H/a=1.
    for branch, original in enumerate(gates, 1):
        gp = replace(original, gate_sill_z=0.)
        opening = .125 / gp.gate_max_opening
        a = gp.gate_max_opening * opening
        for ratio in (.5, 1., 1.125, 1.25, 1.5):
            name = f'gate{branch}_H_over_a_{ratio:g}'
            def test(gp=gp, opening=opening, a=a, ratio=ratio, branch=branch):
                row = dict(branch=branch, H_over_a=ratio, H_m=ratio*a, a_m=a)
                try:
                    flow = sluice_gate_flow(ratio*a, -1., opening, gp)
                except UnsupportedGateRegime as exc:
                    row.update(implemented=False, regime=exc.regime, discharge_m3_s=None)
                    coverage.append(row)
                    assert ratio < 1.25
                    expected = 'weir_or_partially_wetted_opening' if ratio <= 1 else 'weir_gate_transition'
                    assert exc.regime == expected
                    return dict(expected_rejection=True, physical_coverage_pass=False, **row)
                row.update(implemented=True, regime=flow.regime, discharge_m3_s=flow.discharge)
                coverage.append(row)
                assert ratio >= 1.25 and flow.regime == 'free' and flow.discharge > 0
                close(flow.discharge, gp.Cd_gate_free*gp.gate_width*a*math.sqrt(2*gp.gravity*ratio*a))
                return dict(expected_rejection=False, physical_coverage_pass=True, **row)
            check(name, 'head_opening_coverage', test)

        def zeros(gp=gp):
            cases = [(1., -1., 0., 'closed'), (0., -1., .5, 'dry'),
                     (-.1, -1., .5, 'dry'), (1., 1., .5, 'adverse_head'), (1., 2., .5, 'adverse_head')]
            for up, down, opening, expected in cases:
                f = sluice_gate_flow(up, down, opening, gp)
                assert f.regime == expected and f.discharge == 0
            return dict(cases=len(cases), scope='One-way zero routes; low-head wet limit remains unsupported')
        check(f'gate{branch}_closed_dry_adverse', 'gate_limits', zeros)

        def scaling(gp=gp, original=original):
            f = sluice_gate_flow(2., 0., .25, gp)
            width = sluice_gate_flow(2., 0., .25, replace(gp, gate_width=2*gp.gate_width))
            half = sluice_gate_flow(2., 0., .125, gp)
            close(width.discharge, 2*f.discharge)
            close(half.discharge, .5*f.discharge)
            shifted = sluice_gate_flow(original.gate_sill_z+2., original.gate_sill_z, .25, original)
            close(shifted.discharge, f.discharge)
            # Unit covariance: every length and gravity's length unit scaled 100.
            s = 100.
            metric = replace(gp, gravity=gp.gravity*s, gate_width=gp.gate_width*s,
                             gate_max_opening=gp.gate_max_opening*s)
            converted = sluice_gate_flow(2*s, 0., .25, metric)
            close(converted.discharge, f.discharge*s**3, atol=1e-7)
            flows = [sluice_gate_flow(2., tail, .25, gp).discharge for tail in (1.6, 1.9, 1.99, 2.)]
            assert all(v >= 0 for v in flows) and all(a>b for a,b in zip(flows[:-1], flows[1:]))
            return dict(unit='m3/s', length_unit_scaling_factor=s, discharge_scaling_factor=converted.discharge/f.discharge,
                        approaches_zero_at_equal_head=True)
        check(f'gate{branch}_dimensions_scaling_limits', 'gate_units', scaling)

        def transitions(gp=gp):
            # H=100 m is an arithmetic fixture, not a study trajectory.
            H, opening, eps = 100., .25, 1e-10
            samples = [(s, sluice_gate_flow(H, H*s, opening, gp)) for s in (.67-eps,.67,.67+eps,.8-eps,.8,.8+eps)]
            q = [f.discharge for s,f in samples]
            ratio = q[2]/q[1]
            close(ratio, math.sqrt(.99), atol=2e-9)
            close(q[3], q[4], atol=2e-7)
            assert all(v>0 for v in q)
            return dict(lower_limit_Q_m3_s=q[1], lower_right_Q_m3_s=q[2],
                        lower_relative_jump=ratio-1, lower_boundary_continuous=False,
                        upper_boundary_continuous=True, upper_absolute_difference_m3_s=abs(q[3]-q[4]),
                        known_limitation='Published rounded 0.67 causes sqrt(0.99)-1 jump; not smoothed',
                        samples=[dict(submergence=s, discharge_m3_s=f.discharge, regime=f.regime) for s,f in samples])
        check(f'gate{branch}_submergence_boundaries', 'continuity', transitions)

        def bad_inputs(gp=gp):
            for opening in (-.01, 1.01, float('nan')):
                try: sluice_gate_flow(2., 0., opening, gp)
                except GateInputError: pass
                else: raise AssertionError('Invalid opening was not rejected')
            return dict(rejections=3)
        check(f'gate{branch}_invalid_inputs_not_clipped', 'input_validation', bad_inputs)

    theta = np.r_[fixed_gates, p.rho]
    for case, rhs, initial, output_map in [('wetland',wet.make_b2_rhs(gates),x0.copy(),a2.Y),
                                         ('no_wetland',dry.make_rhs(gates),x0[dry.SHARED].copy(),dry.Y)]:
        def mass_identity(case=case, rhs=rhs, initial=initial, om=output_map):
            samples = [(initial.copy(), u0.copy(), 'unchanged_initial_state_and_input')]
            for groundwater_boundary in (97.,97.5):
                x, u = initial.copy(), u0.copy()
                names = a2.STATE_NAMES if case=='wetland' else dry.STATE_NAMES
                idx = {v:i for i,v in enumerate(names)}
                x[idx['S_u']] *= 2
                x[idx['S_m']] *= .9
                x[idx['S_d1']] *= 1.1
                x[idx['S_d2']] *= 1.1
                u[a2.U['P']] = .003/3600
                u[a2.U['R_s']] = 450.
                u[a2.U['Q_ext']] = .27
                u[a2.U['h_g']] = groundwater_boundary
                samples.append((x,u,f'perturbed_boundary_hg_{groundwater_boundary}'))
            residuals=[]
            for x,u,label in samples:
                before=x.copy(); dx,y=rhs(0.,x,u,p)
                independent = u[a2.U['P']]*p.A_rain+u[a2.U['Q_ext']]+y[om['Q_gx']]-y[om['Q_out']]-y[om['E_tot']]
                close(dx.sum(),independent,atol=2e-13,rtol=0)
                close(y[om['external_net']],independent,atol=2e-13,rtol=0)
                assert np.array_equal(before,x)
                residuals.append(dict(fixture=label, derivative_balance_m3_s=float(dx.sum()-independent)))
            dry_u=u0.copy();dry_u[a2.U['P']]=0.
            rain_u=dry_u.copy();rain_u[a2.U['P']]=.001/3600
            d0,_=rhs(0.,initial,dry_u,p); d1,_=rhs(0.,initial,rain_u,p)
            expected = (np.array([p.A_u,0.,*p.A_w,*p.A_d,0.,p.A_r]) if case=='wetland'
                        else np.array([p.A_u,0.,*p.A_d,0.,p.A_r+sum(p.A_w)]))*(.001/3600)
            close(d1-d0,expected,atol=2e-13,rtol=0)
            return dict(samples=residuals, rainfall_1mm_per_hour_total_m3_s=float((d1-d0).sum()),
                        expected_1mm_per_hour_total_m3_s=float(p.A_rain*.001/3600),
                        scope='Instantaneous full RHS conservation identity, not integrated numerical mass balance')
        check(f'{case}_full_ODE_mass_identity_and_rain_units', 'ODE_mass', mass_identity)

        def safety(case=case,rhs=rhs,initial=initial,om=output_map):
            x,u=initial.copy(),u0.copy()
            names = a2.STATE_NAMES if case=='wetland' else dry.STATE_NAMES
            idx={v:i for i,v in enumerate(names)}
            for i in range(2):
                x[idx[f'S_d{i+1}']]=(config['safety']['dam_m'][i]+.02-p.z_d[i])*p.A_d[i]
                if case=='wetland':
                    x[idx[f'S_w{i+1}']]=(config['safety']['wetland_m'][i]+.02-p.z_w[i])*p.A_w[i]
            before=x.copy();dx,y=rhs(0.,x,u,p)
            assert np.array_equal(before,x)
            kwargs=(np.array([0.]),x[:,None],y[:,None],u[:,None])
            rows=(wet.constraint_audit(*kwargs,theta,p,config)[1] if case=='wetland' else dry.audit_six(*kwargs,p,config))
            safety_rows=[r for r in rows if r['label'].endswith('_safe')]
            assert len(safety_rows)==(4 if case=='wetland' else 2)
            assert all(not bool(r['pass_at_tolerance']) and r['margin']<0 for r in safety_rows)
            if case=='wetland':
                assert y[a2.Y['Q_wsp1']]>0 and y[a2.Y['Q_wsp2']]>0
                assert y[a2.Y['Q_dsp1']]==0 and y[a2.Y['Q_dsp2']]==0
            return dict(RHS_completed=True, safety_feasible=False, storage_unchanged_by_RHS=True,
                        safety_margins_m={r['label']:float(r['margin']) for r in safety_rows},
                        scope='Artificial safety-exceeding states distinguish safety outcomes from hydraulic validity')
        check(f'{case}_safety_exceedance_is_not_storage_clipping', 'safety_classification', safety)

    def native(v):
        if isinstance(v,np.generic):return v.item()
        if isinstance(v,np.ndarray):return v.tolist()
        raise TypeError(type(v).__name__)
    sources=['B/gate_hydraulics.py','B/no_spillway_model.py','B/wetland_counterfactual.py',
             'B/no_spillway_config.json','B/wetland_only_config.json']
    summary=dict(regression_checks=len(records),regression_passed=sum(bool(r['regression_pass']) for r in records),
                 regression_failed=sum(not r['regression_pass'] for r in records),
                 requested_low_head_physical_coverage_pass=False,
                 hydraulic_G5_resolved=False,
                 explanation='Expected low-head rejections are passing regression tests, not supported low-head physics. Known 0.67 jump remains.',
                 original_physical_configuration=True, trajectory_integration_performed=False,
                 sources_sha256={s:hashlib.sha256((root/s).read_bytes()).hexdigest() for s in sources},
                 initial_state_wetland_m3=x0.tolist(),initial_state_no_wetland_m3=x0[dry.SHARED].tolist(),
                 records=records)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    (args.output_dir/'gate_regression_results.json').write_text(json.dumps(summary,indent=2,default=native)+'\n')
    for filename,rows in [('gate_regression_results.csv',records),('gate_head_opening_coverage.csv',coverage)]:
        keys=sorted(set().union(*(r.keys() for r in rows)))
        with (args.output_dir/filename).open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=keys);writer.writeheader()
            for row in rows:
                writer.writerow({k:json.dumps(v,default=native) if isinstance(v,(dict,list)) else v for k,v in row.items()})
    print(json.dumps({k:summary[k] for k in ('regression_checks','regression_passed','regression_failed','requested_low_head_physical_coverage_pass','hydraulic_G5_resolved')},indent=2))
    return 0 if summary['regression_failed']==0 else 1

if __name__=='__main__':
    raise SystemExit(main())
