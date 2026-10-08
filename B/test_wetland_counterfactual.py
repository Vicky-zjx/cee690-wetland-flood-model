"""Focused six-state checks; artificial fixtures are not study observations.

Short/stress experiments here never alter approved files or saved results.
"""
from __future__ import annotations
import copy
import unittest
from unittest.mock import patch
import numpy as np
from A import wetland_flood_model as a2
from B import optimize_wetland as old
from B import no_spillway_model as wet
from B import wetland_counterfactual as case


class SixStateCounterfactualTests(unittest.TestCase):
    def setUp(self):
        self.config = old.load_config(old.ROOT/'no_spillway_config.json')
        self.params = wet.base_parameters(self.config)
        self.gates = wet.gate_parameters(self.config)
        self.x = np.asarray(self.config['initial_conditions']['x0_m3'])[case.SHARED].copy()
        self.u = np.array([.1, .1, 0., 20., 100., 2., 0., 97.15, 97.15, 0.])
        self.rhs = case.make_rhs(self.gates)

    def short_fixture(self):
        """TEST-ONLY four-minute constant-weather experiment with closed gates."""
        cfg = copy.deepcopy(self.config)
        cfg['numerics']['duration_s'] = 240.
        for branch in (1, 2):
            cfg['controls'][f'gate{branch}_hours'] = [0., 240./3600.]
            cfg['controls'][f'gate{branch}_values'] = [0., 0.]
        t = np.array([0., 120., 240.])
        return cfg, np.zeros(4), (t, np.repeat(self.u[2:, None], 3, axis=1))

    def test_six_states_have_no_dummy_wetland_or_spill_quantities(self):
        self.assertEqual(case.STATE_NAMES, ['S_u', 'S_m', 'S_d1', 'S_d2', 'S_g', 'S_r'])
        self.assertEqual(case.SHARED, [0, 1, 4, 5, 6, 7])
        dx, y = self.rhs(0., self.x, self.u, self.params)
        self.assertEqual(dx.shape, (6,))
        self.assertEqual(y.shape, (len(case.OUTPUT_NAMES),))
        for name in ('S_w1', 'S_w2', 'H_w1', 'H_w2', 'Q_wr1', 'Q_wsp1', 'Q_wg1', 'E_w1', 'Q_dsp1'):
            self.assertNotIn(name, case.STATE_NAMES+case.OUTPUT_NAMES)
        with self.assertRaises(FloatingPointError):
            self.rhs(0., np.asarray(self.config['initial_conditions']['x0_m3']), self.u, self.params)
        self.assertEqual(y[case.Y['S_tot']], self.x.sum())
        self.assertEqual(y[case.Y['H_r']], self.params.z_r+self.x[case.X['S_r']]/self.params.A_r)

    def test_internal_transfers_cancel_for_both_groundwater_exchange_directions(self):
        for exterior_head, expected_sign in ((97.0, -1), (97.5, 1)):
            with self.subTest(exterior_head=exterior_head):
                u = self.u.copy()
                u[a2.U['P']] = 3e-6
                u[a2.U['Q_ext']] = .27
                u[a2.U['R_s']] = 450.
                u[a2.U['h_g']] = exterior_head
                dx, y = self.rhs(456., self.x, u, self.params)
                self.assertEqual(np.sign(y[case.Y['Q_gx']]), expected_sign)
                total_et = sum(y[case.Y[name]] for name in ('E_m', 'E_d1', 'E_d2', 'E_r'))
                net = (u[a2.U['P']]*self.params.A_rain+u[a2.U['Q_ext']]
                       +y[case.Y['Q_gx']]-y[case.Y['Q_out']]-total_et)
                self.assertAlmostEqual(float(dx.sum()), float(net), places=13)
                self.assertAlmostEqual(float(y[case.Y['external_net']]), float(net), places=13)
                self.assertAlmostEqual(float(y[case.Y['E_tot']]), float(total_et), places=13)
                river = (u[a2.U['P']]*(self.params.A_r+sum(self.params.A_w))
                         +sum(y[case.Y[name]] for name in ('Q_ur', 'Q_G1', 'Q_G2', 'Q_gr'))
                         +u[a2.U['Q_ext']]-y[case.Y['Q_out']]-y[case.Y['E_r']])
                self.assertAlmostEqual(float(dx[case.X['S_r']]), float(river), places=13)

    def test_rainfall_units_preserve_footprint_without_inflating_river_geometry(self):
        depth_rate = .001/3600.  # TEST-ONLY 1 mm/hour, represented in m/s.
        dry = self.u.copy(); dry[a2.U['P']] = 0.
        rainy = dry.copy(); rainy[a2.U['P']] = depth_rate
        dx0, _ = self.rhs(0., self.x, dry, self.params)
        dx1, y1 = self.rhs(0., self.x, rainy, self.params)
        areas = np.array([self.params.A_u, 0., *self.params.A_d, 0.,
                          self.params.A_r+sum(self.params.A_w)])
        np.testing.assert_allclose(dx1-dx0, depth_rate*areas, rtol=0., atol=1e-14)
        self.assertEqual(float(areas.sum()), self.params.A_rain)
        self.assertAlmostEqual(float((dx1-dx0).sum()), depth_rate*660000., places=13)
        self.assertAlmostEqual(y1[case.Y['P_former_wetland']], depth_rate*sum(self.params.A_w), places=13)
        # Artificially changing only former-wetland footprint changes rain
        # routing, never stage-storage area, evaporation area or aquifer area.
        altered = copy.deepcopy(self.params); altered.A_w *= 2
        _, changed = self.rhs(0., self.x, rainy, altered)
        for name in ('H_r', 'E_r', 'Q_out'):
            self.assertEqual(changed[case.Y[name]], y1[case.Y[name]])
        self.assertEqual(altered.A_r, self.params.A_r)
        self.assertEqual(altered.A_g, self.params.A_g)

    def test_gate_uses_internal_river_tailwater_and_responds_to_it(self):
        # TEST-ONLY domain expansion distinguishes internal Hr from external hd
        # with a submerged stress fixture, not an approved event configuration.
        p = copy.deepcopy(self.params); p.H_r_valid_max = 102.
        x = self.x.copy(); x[case.X['S_r']] = (100.25-p.z_r)*p.A_r
        u = self.u.copy(); u[a2.U['h_d']] = 96.5
        with patch.object(case, 'sluice_gate_flow', wraps=case.sluice_gate_flow) as gate:
            _, y = self.rhs(0., x, u, p)
        self.assertEqual(gate.call_count, 2)
        for call in gate.call_args_list:
            self.assertEqual(call.args[1], y[case.Y['H_r']])
            self.assertNotEqual(call.args[1], u[a2.U['h_d']])
        expected = case.sluice_gate_flow(y[case.Y['H_d1']], y[case.Y['H_r']], u[0], self.gates[0])
        self.assertEqual(expected.regime, 'submerged')
        self.assertEqual(y[case.Y['Q_G1']], expected.discharge)
        self.assertEqual(y[case.Y['Q_G2']], 0.)  # Adverse internal head; no reverse flow.
        _, lower_river = self.rhs(0., self.x, u, p)
        self.assertGreater(lower_river[case.Y['Q_G1']], y[case.Y['Q_G1']])
        self.assertGreater(lower_river[case.Y['Q_G2']], y[case.Y['Q_G2']])

    def test_no_wetland_design_or_flux_dependency_remains(self):
        expected_dx, expected_y = self.rhs(0., self.x, self.u, self.params)
        p = copy.deepcopy(self.params)
        # A_w survives as a rain footprint; dynamic wetland attributes do not.
        for name in ('z_w', 'z_wo', 'z_wsp', 'c_wsp', 'C_wg', 'V_w', 'Q_wrcap',
                     'R_bare', 'K_bare', 'R_veg', 'K_veg', 'rho', 'tau', 'w_depth_init'):
            if hasattr(p, name): delattr(p, name)
        with patch.object(p, 'design', side_effect=AssertionError('No wetland design evaluation')), \
             patch.object(a2, 'system', side_effect=AssertionError('No eight-state A.2 RHS')), \
             patch.object(wet, 'make_b2_rhs', side_effect=AssertionError('No eight-state B.2 RHS')):
            dx, y = self.rhs(0., self.x, self.u, p)
            case.validate_six(p, np.repeat(self.u[:, None], len(p.time), axis=1), self.x, self.gates)
        np.testing.assert_array_equal(dx, expected_dx)
        np.testing.assert_array_equal(y, expected_y)

    def test_common_flux_definitions_match_accepted_wetland_model(self):
        eight = np.asarray(self.config['initial_conditions']['x0_m3'])
        _, expected = wet.make_b2_rhs(self.gates)(0., eight, self.u, self.params)
        _, actual = self.rhs(0., self.x, self.u, self.params)
        for name in ('Q_inf', 'Q_mg', 'Q_gr', 'Q_gx', 'Q_out', 'Q_ud1', 'Q_ud2',
                     'Q_ur', 'E_m', 'E_d1', 'E_d2', 'E_r', 'H_d1', 'H_d2', 'H_g', 'H_r'):
            self.assertEqual(actual[case.Y[name]], expected[a2.Y[name]], name)

    def test_official_rk_stage_g5_rejection_has_complete_diagnostic(self):
        cfg, schedule, forcing = self.short_fixture()
        cfg['initial_conditions']['x0_m3'][0] = 0.
        cfg['initial_conditions']['x0_m3'][4] = .03755*self.params.A_d[0]
        schedule[:2] = .1  # a=.03 m; initial head slightly above 1.25*a.
        p = wet.base_parameters(cfg)
        x = np.asarray(cfg['initial_conditions']['x0_m3'])[case.SHARED]
        u = self.u.copy(); u[:2] = [.1, 0.]
        dx, _ = self.rhs(0., x, u, p)  # Initial gate state is admissible.
        stage_depth = ((x[case.X['S_d1']]+60.*dx[case.X['S_d1']])/p.A_d[0]
                       +p.z_d[0]-self.gates[0].gate_sill_z)
        self.assertLess(stage_depth, 1.25*.03)
        with self.assertRaises(wet.HydraulicApplicabilityFailure) as raised:
            case.simulate_bypass(cfg, schedule, dt=120, forcing=forcing)
        d = raised.exception.diagnostic
        self.assertEqual(d['branch'], 1)
        self.assertEqual(d['time_s'], 60.)  # Official RK half-step, not a saved node.
        self.assertAlmostEqual(d['water_depth_m'], stage_depth, places=12)
        self.assertAlmostEqual(d['opening_height_m'], .03, places=14)
        self.assertAlmostEqual(d['residual_m'], .0375-stage_depth, places=12)
        self.assertGreater(d['residual_m'], 0.)
        self.assertIn('tailwater_stage_m', d)
        self.assertIn('RK stage', d['detection'])

    def test_official_run_preserves_shared_initial_water_and_balance(self):
        cfg, schedule, forcing = self.short_fixture()
        run = case.simulate_bypass(cfg, schedule, dt=120, forcing=forcing)
        full_initial = np.asarray(cfg['initial_conditions']['x0_m3'])
        self.assertEqual(run['x'].shape, (6, 3))
        np.testing.assert_array_equal(run['x'][:, 0], full_initial[case.SHARED])
        self.assertAlmostEqual(full_initial.sum()-run['x'][:, 0].sum(), full_initial[2:4].sum(), places=8)
        self.assertEqual(run['x'][case.X['S_r'], 0], full_initial[7])
        self.assertLess(run['metrics']['derivative_balance_max_m3s'], 1e-12)
        self.assertLess(run['metrics']['balance_relative'], self.params.balance_rtol)
        self.assertTrue(all(row['strict_pass'] for row in run['audit']))
        self.assertTrue(all('S_w' not in row['label'] and 'H_w' not in row['label'] for row in run['audit']))
        expected = (run['x'].sum(axis=0)-run['x'][:, 0].sum()
                    -a2.cumulative_trapezoid(run['y'][case.Y['external_net']], run['time']))
        np.testing.assert_array_equal(run['balance'], expected)

    def test_high_dams_remain_unclipped_valid_but_safety_infeasible(self):
        cfg, schedule, forcing = self.short_fixture()
        cfg['initial_conditions']['x0_m3'][4:6] = (self.params.A_d*2.).tolist()
        forcing[1][0] = 1e-5  # TEST-ONLY rain drives already-unsafe pools upward.
        run = case.simulate_bypass(cfg, schedule, dt=120, forcing=forcing)
        self.assertTrue(np.all(run['x'][2:4, -1] > run['x'][2:4, 0]))
        for branch in (1, 2):
            row = next(r for r in run['audit'] if r['label'] == f'H_d{branch}_safe')
            self.assertFalse(row['pass_at_tolerance'])
            self.assertGreater(row['scaled_residual'], 0.)
            self.assertEqual(row['worst_time_s'], 240.)
        self.assertLess(run['metrics']['derivative_balance_max_m3s'], 1e-12)
        self.assertTrue(next(r for r in run['audit'] if r['label'] == 'integration_status')['strict_pass'])

    def test_final_point_audit_uses_correct_state_mapping_and_shared_bounds(self):
        t = np.array([0., 120.])
        x = np.repeat(self.x[:, None], 2, axis=1)
        x[case.X['S_d1'], -1] = self.params.A_d[0]*(self.config['safety']['dam_m'][0]-self.params.z_d[0]+.01)
        u = np.repeat(self.u[:, None], 2, axis=1)
        y = np.column_stack([self.rhs(t[k], x[:, k], u[:, k], self.params)[1] for k in range(2)])
        rows = case.audit_six(t, x, y, u, self.params, self.config)
        dam = next(r for r in rows if r['label'] == 'H_d1_safe')
        self.assertEqual(dam['worst_time_s'], 120.)
        self.assertAlmostEqual(dam['margin'], -.01, places=12)
        self.assertFalse(dam['strict_pass'])
        self.assertEqual(next(r for r in rows if r['label'] == 'nonnegative_S_d1')['scale'],
                         self.config['numerics']['volume_scales_m3'][4])
        labels = {r['label'] for r in rows}
        self.assertTrue({'groundwater_geometry', 'groundwater_soil', 'river_geometry',
                         'gate1_head_opening', 'gate2_head_opening', 'soil_capacity'} <= labels)
        self.assertFalse(any('wetland' in label or 'S_w' in label or 'H_w' in label for label in labels))


if __name__ == '__main__':
    unittest.main()
