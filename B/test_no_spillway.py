"""Targeted revised-model checks; artificial fixtures are not study observations."""
from __future__ import annotations

import copy
import unittest
from unittest.mock import patch
import numpy as np
from A import wetland_flood_model as a2
from B import optimize_wetland as old
from B import no_spillway_model as revised


class NoDamSpillwayTests(unittest.TestCase):
    def setUp(self):
        self.cfg = old.load_config()
        self.params = revised.base_parameters(self.cfg)
        self.gates = revised.gate_parameters(self.cfg)
        self.x = self.params.x_init.copy()
        self.u = a2.make_inputs(self.params)[:, 0].copy()
        self.u[:2] = .1
        self.rhs = revised.make_b2_rhs(self.gates)

    def test_removed_parameters_have_no_active_factory_defaults_or_overrides(self):
        for name in revised.REMOVED_DAM_SPILLWAY_PARAMETERS:
            self.assertFalse(hasattr(self.params, name))
            cfg = copy.deepcopy(self.cfg)
            cfg['parameter_overrides'][name] = dict(value=[999., 999.], approved=True,
                source_id='TEST_ONLY', reason='Check removed configuration is rejected.')
            with self.assertRaisesRegex(revised.PreflightError, 'not part of the revised model'):
                revised.base_parameters(cfg)

    def test_other_fluxes_match_previous_model_and_only_spill_links_are_removed(self):
        # High pool levels activate the historical spill law in the comparator.
        # New dynamics must never evaluate that law or depend on its parameters.
        x = self.x.copy()
        x[4:6] = self.params.A_d*np.array([1.35, 1.40])
        oldparams = old.base_parameters(self.cfg)
        oldparams.rho = self.params.rho.copy()
        olddx, oldy = old.make_b2_rhs(self.gates)(0., x, self.u, oldparams)
        spill = oldy[[a2.Y['Q_dsp1'], a2.Y['Q_dsp2']]]
        self.assertTrue(np.all(spill > 0.))
        with patch.object(a2, 'system', side_effect=AssertionError('Historical RHS must not run')):
            dx, y = self.rhs(0., x, self.u, self.params)
        kept = [j for j in range(len(a2.OUTPUT_NAMES)) if j not in (a2.Y['Q_dsp1'], a2.Y['Q_dsp2'])]
        np.testing.assert_array_equal(y[kept], oldy[kept])
        np.testing.assert_array_equal(y[[a2.Y['Q_dsp1'], a2.Y['Q_dsp2']]], [0., 0.])
        expected = olddx.copy()
        expected[2:4] -= spill
        expected[4:6] += spill
        np.testing.assert_allclose(dx, expected, rtol=0., atol=1e-14)

    def test_transfer_conservation_with_both_directions_of_boundary_exchange(self):
        x = self.x.copy()
        x[2:4] = self.params.A_w*np.array([.9, 1.0])
        x[4:6] = self.params.A_d*np.array([1.4, 1.6])
        for exterior_head in [97.0, 97.5]:
            u = self.u.copy()
            u[a2.U['P']] = 3e-6
            u[a2.U['h_g']] = exterior_head
            u[a2.U['Q_ext']] = .27
            u[a2.U['R_s']] = 450.
            dx, y = self.rhs(456., x, u, self.params)
            expected = (u[a2.U['P']]*self.params.A_rain + u[a2.U['Q_ext']] +
                        y[a2.Y['Q_gx']] - y[a2.Y['Q_out']] - y[a2.Y['E_tot']])
            self.assertAlmostEqual(float(dx.sum()), float(expected), places=13)
            self.assertAlmostEqual(float(y[a2.Y['external_net']]), float(expected), places=13)

    def test_rainfall_depth_rate_times_disjoint_area_gives_storage_rate(self):
        # 1 mm/hour in m/s: every recipient gets P*A m3/s; soil/aquifer are not
        # second rainfall footprints. The test isolates the units algebraically.
        delta = .001/3600.
        u = self.u.copy()
        u[a2.U['P']] = 0.
        dx0, y0 = self.rhs(0., self.x, u, self.params)
        u[a2.U['P']] = delta
        dx1, y1 = self.rhs(0., self.x, u, self.params)
        areas = np.array([self.params.A_u, 0., *self.params.A_w,
                          *self.params.A_d, 0., self.params.A_r])
        np.testing.assert_allclose(dx1-dx0, delta*areas, rtol=0., atol=1e-14)
        self.assertAlmostEqual(float((dx1-dx0).sum()), delta*self.params.A_rain, places=13)
        self.assertAlmostEqual(float(y1[a2.Y['P_total']]-y0[a2.Y['P_total']]), delta*660000., places=13)

    def test_high_dam_storage_is_not_capped_clipped_or_bypassed(self):
        x = self.x.copy()
        x[4:6] = self.params.A_d*2.0  # Above historical crests and fixed safety levels.
        before = x.copy()
        u = self.u.copy()
        u[:2] = 0.
        u[a2.U['P']] = 1e-5
        dx, y = self.rhs(0., x, u, self.params)
        np.testing.assert_array_equal(x, before)
        np.testing.assert_allclose(y[[a2.Y['H_d1'], a2.Y['H_d2']]], self.params.z_d+2.)
        self.assertTrue(np.all(dx[4:6] > 0.))
        np.testing.assert_array_equal(y[[a2.Y['Q_dsp1'], a2.Y['Q_dsp2'], a2.Y['Q_G1'], a2.Y['Q_G2']]], np.zeros(4))
        self.params.d_depth_init = np.array([2., 2.])
        full_u = np.repeat(u[:, None], len(self.params.time), axis=1)
        revised.validate(self.params, full_u, self.gates)  # No spillcrest/safety clipping check.
        t = np.array([0., 120.])
        theta = old.reference_decision(self.cfg)
        theta[:-2] = 0.
        _, rows = revised.constraint_audit(t, np.repeat(x[:, None], 2, axis=1),
            np.repeat(y[:, None], 2, axis=1), np.repeat(u[:, None], 2, axis=1), theta, self.params, self.cfg)
        damrows = [r for r in rows if r['label'] in ('H_d1_safe', 'H_d2_safe')]
        self.assertTrue(all(r['scaled_residual'] > 0. for r in damrows))
        self.assertFalse(revised.is_feasible({'audit': rows}, self.cfg))

    def test_wetland_overflow_is_retained_and_enters_river_exactly_once(self):
        x = self.x.copy()
        x[2:4] = self.params.A_w*np.array([.9, 1.0])
        u = self.u.copy()
        u[:2] = 0.
        dx, y = self.rhs(0., x, u, self.params)
        flow = y[[a2.Y['Q_wsp1'], a2.Y['Q_wsp2']]]
        self.assertTrue(np.all(flow > 0.))
        changed = copy.deepcopy(self.params)
        changed.c_wsp *= 2.
        dx2, y2 = self.rhs(0., x, u, changed)
        delta = np.zeros(8)
        delta[2:4] = -flow
        delta[7] = flow.sum()
        np.testing.assert_allclose(dx2-dx, delta, rtol=0., atol=1e-14)
        self.assertEqual(y2[a2.Y['external_net']], y[a2.Y['external_net']])

    def test_poisoned_spill_parameters_do_not_affect_rhs_validation_or_feasibility(self):
        # Even maliciously reintroduced historical attributes must remain unread.
        dx, y = self.rhs(0., self.x, self.u, self.params)
        poisoned = copy.deepcopy(self.params)
        poisoned.z_dsp = np.array([np.nan, -np.inf])
        poisoned.c_dsp = np.array([np.inf, -999.])
        full_u = np.repeat(self.u[:, None], len(self.params.time), axis=1)
        with patch.object(a2, 'validate', side_effect=AssertionError('Historical validation must not run')):
            revised.validate(poisoned, full_u, self.gates)
        dx2, y2 = self.rhs(0., self.x, self.u, poisoned)
        np.testing.assert_array_equal(dx2, dx)
        np.testing.assert_array_equal(y2, y)
        args = (np.array([0., 120.]), np.repeat(self.x[:, None], 2, axis=1),
                np.repeat(y[:, None], 2, axis=1), np.repeat(self.u[:, None], 2, axis=1),
                old.reference_decision(self.cfg))
        first = revised.constraint_audit(*args, self.params, self.cfg)
        second = revised.constraint_audit(*args, poisoned, self.cfg)
        np.testing.assert_array_equal(first[0], second[0])
        self.assertEqual(first[1], second[1])

    def test_gate_applicability_guard_remains_a_diagnostic_rejection(self):
        x = self.x.copy()
        x[4] = self.params.A_d[0]*.1
        u = self.u.copy()
        u[0] = .5
        with self.assertRaises(revised.HydraulicApplicabilityFailure) as raised:
            self.rhs(123., x, u, self.params)
        diagnostic = raised.exception.diagnostic
        self.assertEqual(diagnostic['branch'], 1)
        self.assertEqual(diagnostic['time_s'], 123.)
        self.assertAlmostEqual(diagnostic['residual_m'], .0875)

    def test_preflight_never_calls_historical_preflight_or_reads_old_spillcrest(self):
        with patch.object(old, 'preflight', side_effect=AssertionError('Historical preflight must not run')), \
             patch.object(revised, 'fixed_forcing', return_value=None):
            self.assertEqual(revised.preflight(self.cfg), [])

    def test_short_official_integration_preserves_mass_and_uncapped_dam_growth(self):
        # TEST-ONLY four-minute horizon, zero gates, constant rain/meteorology.
        cfg = copy.deepcopy(self.cfg)
        cfg['numerics']['duration_s'] = 240.
        cfg['controls']['gate1_hours'] = cfg['controls']['gate2_hours'] = [0., 240./3600.]
        cfg['controls']['gate1_values'] = cfg['controls']['gate2_values'] = [0., 0.]
        cfg['initial_conditions']['x0_m3'][4:6] = (self.params.A_d*2.).tolist()
        theta = revised.reference_decision(cfg)
        source_t = np.array([0., 120., 240.])
        one = np.array([1e-5, 20., 100., 2., 0., 97.15, 97.15, 0.])
        forcing = (source_t, np.repeat(one[:, None], 3, axis=1))
        with patch.object(a2, 'system', side_effect=AssertionError('No A.2 RHS delegation')), \
             patch.object(a2, 'validate', side_effect=AssertionError('No A.2 validation delegation')):
            run = revised.simulate(theta, ('type_a', 'type_a'), cfg, forcing=forcing)
        self.assertTrue(np.all(run['x'][4:6, -1] > run['x'][4:6, 0]))
        np.testing.assert_array_equal(run['y'][[a2.Y['Q_dsp1'], a2.Y['Q_dsp2']]], np.zeros((2, 3)))
        self.assertLess(run['metrics']['derivative_balance_max_m3s'], 1e-12)
        self.assertLess(run['metrics']['balance_relative'], self.params.balance_rtol)
        self.assertFalse(revised.is_feasible(run, cfg))  # Valid integration, fixed safety violated.


if __name__ == '__main__':
    unittest.main()
