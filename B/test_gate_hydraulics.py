"""Numerical fixtures ONLY: these are not B.2 geometry or calibrated coefficients.

Run from the project root: python3 -m unittest B.test_gate_hydraulics -v
No solver, optimization, scientific plot, or A.2 output is involved.
"""
from dataclasses import replace
import math
import unittest

from B.gate_hydraulics import (GateConfigurationError, GateInputError,
                              GateParameters, UnsupportedGateRegime,
                              sluice_gate_flow)


def fixture():
    # Deliberate arithmetic fixture: gravity=10 makes selected expected values exact.
    return GateParameters(gravity=10, gate_width=2, gate_max_opening=1,
                          gate_sill_z=0, Cd_gate_free=0.6,
                          Cd_gate_submerged=0.8,
                          neglect_approach_velocity_head=True)


class GateHydraulicsTests(unittest.TestCase):
    def setUp(self):
        self.p = fixture()

    def test_explicit_missing_configuration_blocks_even_before_simulation(self):
        with self.assertRaisesRegex(GateConfigurationError, "gravity"):
            GateParameters()
        for field in ("gravity", "gate_width", "gate_max_opening", "gate_sill_z",
                      "Cd_gate_free", "Cd_gate_submerged"):
            with self.subTest(field=field):
                with self.assertRaisesRegex(GateConfigurationError, field):
                    replace(self.p, **{field: None})

    def test_invalid_fixed_values_and_unapproved_head_assumption(self):
        for field in ("gravity", "gate_width", "gate_max_opening",
                      "Cd_gate_free", "Cd_gate_submerged"):
            with self.subTest(field=field):
                with self.assertRaises(GateConfigurationError):
                    replace(self.p, **{field: 0})
        for value in (None, False):
            with self.assertRaisesRegex(GateConfigurationError, "approach"):
                replace(self.p, neglect_approach_velocity_head=value)
        with self.assertRaises(GateConfigurationError):
            replace(self.p, gate_sill_z=math.nan)

    def test_closed_dry_and_adverse_head_are_explicit_zero_routes(self):
        cases = [(5, 0, 0, "closed"), (0, -1, .4, "dry"),
                 (-1, -2, .4, "dry"), (5, 5, .4, "adverse_head"),
                 (5, 6, .4, "adverse_head")]
        for up, down, opening, regime in cases:
            with self.subTest(regime=regime):
                flow = sluice_gate_flow(up, down, opening, self.p)
                self.assertEqual(flow.discharge, 0)
                self.assertEqual(flow.regime, regime)

    def test_free_fixture_uses_sill_head_and_no_second_opening_multiplier(self):
        # H=5, a=.4, width=2, sqrt(2*10*5)=10 => Q=.6*2*.4*10=4.8.
        flow = sluice_gate_flow(5, 0, .4, self.p)
        self.assertEqual(flow.regime, "free")
        self.assertAlmostEqual(flow.discharge, 4.8)
        self.assertAlmostEqual(flow.opening_height, .4)
        half = sluice_gate_flow(5, 0, .2, self.p)
        self.assertAlmostEqual(half.discharge, 2.4)

    def test_maximum_physical_opening_is_in_the_map(self):
        flow = sluice_gate_flow(5, 0, .4, replace(self.p, gate_max_opening=.5))
        self.assertAlmostEqual(flow.opening_height, .2)
        self.assertAlmostEqual(flow.discharge, 2.4)

    def test_submerged_fixture_uses_differential_head_and_its_own_coefficient(self):
        flow = sluice_gate_flow(5, 4, .4, self.p)
        self.assertEqual(flow.regime, "submerged")
        self.assertAlmostEqual(flow.discharge, 2.862167011199731)
        changed_free = sluice_gate_flow(5, 4, .4, replace(self.p, Cd_gate_free=.2))
        self.assertEqual(flow.discharge, changed_free.discharge)
        changed_sub = sluice_gate_flow(5, 4, .4, replace(self.p, Cd_gate_submerged=.4))
        self.assertAlmostEqual(changed_sub.discharge, flow.discharge / 2)

    def test_vertical_datum_shift_leaves_flow_and_heads_unchanged(self):
        for tail in (0, 4):
            original = sluice_gate_flow(5, tail, .4, self.p)
            shifted = sluice_gate_flow(105, tail+100, .4,
                                      replace(self.p, gate_sill_z=100))
            self.assertEqual(original, shifted)

    def test_length_unit_covariance_demonstrates_volume_per_time_dimensions(self):
        # Consistently converting every length and gravity length unit by 100
        # must multiply discharge by 100^3; this is not a physical scenario.
        original = sluice_gate_flow(5, 4, .4, self.p)
        scaled = replace(self.p, gravity=1000, gate_width=200,
                         gate_max_opening=100, gate_sill_z=0)
        other = sluice_gate_flow(500, 400, .4, scaled)
        self.assertAlmostEqual(other.discharge/original.discharge, 1_000_000)

    def test_no_old_empirical_capacity_cap(self):
        low = sluice_gate_flow(5, 0, .4, self.p)
        high = sluice_gate_flow(20, 0, .4, self.p)
        self.assertAlmostEqual(high.discharge, 2 * low.discharge)
        self.assertGreater(low.discharge, 1)  # It does not inherit A.2's 1 m3/s cap.

    def test_partial_opening_and_low_head_transition_are_not_full_orifices(self):
        for upstream, expected in ((.3, "weir_or_partially_wetted_opening"),
                                   (.4, "weir_or_partially_wetted_opening"),
                                   (.45, "weir_gate_transition")):
            with self.subTest(upstream=upstream):
                with self.assertRaises(UnsupportedGateRegime) as raised:
                    sluice_gate_flow(upstream, 0, .4, self.p)
                self.assertEqual(raised.exception.regime, expected)
        self.assertEqual(sluice_gate_flow(.5, 0, .4, self.p).regime, "free")

    def test_submergence_endpoints_and_documented_transition(self):
        self.assertEqual(sluice_gate_flow(100, 67, .4, self.p).regime, "free")
        for tail in (67.00001, 70, 75, 79.99999):
            with self.subTest(tail=tail):
                flow = sluice_gate_flow(100, tail, .4, self.p)
                self.assertEqual(flow.regime, "transition")
                self.assertGreater(flow.discharge, 0)
        self.assertEqual(sluice_gate_flow(100, 80, .4, self.p).regime, "submerged")

    def test_transition_midpoint_fixture_uses_documented_three_dH_head(self):
        # s=.735 -> equal weights; width*a=.8; dH=26.5; gravity=10.
        # Half of .48*sqrt(1590)+.64*sqrt(530) = 16.93690621520256.
        flow = sluice_gate_flow(100, 73.5, .4, self.p)
        self.assertEqual(flow.regime, "transition")
        self.assertAlmostEqual(flow.discharge, 16.93690621520256)

    def test_documented_lower_transition_rounding_is_visible_not_smoothed(self):
        free = sluice_gate_flow(100, 67, .4, self.p).discharge
        just_inside = sluice_gate_flow(100, 67.00000001, .4, self.p).discharge
        # The published .67 rather than exact 2/3 gives sqrt(.99) limiting ratio.
        self.assertAlmostEqual(just_inside/free, math.sqrt(.99), places=8)
        self.assertLess(just_inside, free)

    def test_transition_approaches_submerged_endpoint(self):
        just_inside = sluice_gate_flow(100, 79.99999999, .4, self.p).discharge
        endpoint = sluice_gate_flow(100, 80, .4, self.p).discharge
        self.assertAlmostEqual(just_inside, endpoint, places=7)
        flows = [sluice_gate_flow(100, tail, .4, self.p).discharge
                 for tail in (67, 68, 72, 76, 80)]
        self.assertTrue(all(a > b for a, b in zip(flows[:-1], flows[1:])))

    def test_submerged_head_can_approach_zero_without_negative_discharge(self):
        flows = [sluice_gate_flow(5, tail, .4, self.p).discharge
                 for tail in (4, 4.9, 4.99, 5)]
        self.assertTrue(all(a > b for a, b in zip(flows[:-1], flows[1:])))
        self.assertEqual(flows[-1], 0)

    def test_inputs_rejected_instead_of_clipped(self):
        for opening in (-.01, 1.01, math.nan, math.inf, None):
            with self.subTest(opening=opening):
                with self.assertRaises(GateInputError):
                    sluice_gate_flow(5, 0, opening, self.p)
        for upstream, downstream in ((math.nan, 0), (5, math.inf)):
            with self.assertRaises(GateInputError):
                sluice_gate_flow(upstream, downstream, .4, self.p)


class B2GateTransferTests(unittest.TestCase):
    """RHS-only fixture checks; never an approved project configuration or run."""
    @classmethod
    def setUpClass(cls):
        # Only wrapper tests need the normal project's NumPy/Matplotlib imports.
        # No ode4u or OPT routine is loaded or executed.
        import numpy as np
        from A import wetland_flood_model as a2
        from B.optimize_wetland import make_b2_rhs, TrajectoryFailure
        cls.np, cls.a2 = np, a2
        cls.make_b2_rhs = staticmethod(make_b2_rhs)
        cls.trajectory_failure = TrajectoryFailure

    def setUp(self):
        self.params = self.a2.c()
        self.x = self.params.x_init.copy()
        self.u = self.a2.make_inputs(self.params)[:, 0].copy()
        self.u[:2] = [.4, .6]
        # Deliberately invented UNIT-TEST fixtures, NOT project gate settings.
        # The .3 m initial pool depth fully wets these small fixture openings.
        self.gates = (
            replace(fixture(), gate_max_opening=.2, gate_sill_z=100.0),
            replace(fixture(), gate_width=1.5, gate_max_opening=.2,
                    gate_sill_z=99.8, Cd_gate_free=.65),
        )

    def zero_gate_reference(self, state=None):
        closed = self.u.copy()
        closed[:2] = 0
        return self.a2.system(0, self.x if state is None else state,
                              closed, self.params)

    def test_each_gate_is_evaluated_once_and_conserved_without_other_flux_changes(self):
        from collections import Counter
        from unittest.mock import patch
        dx0, y0 = self.zero_gate_reference()
        expected_q = self.np.array([
            sluice_gate_flow(y0[self.a2.Y[f'H_d{i+1}']],
                             y0[self.a2.Y[f'H_w{i+1}']], self.u[i], gate).discharge
            for i, gate in enumerate(self.gates)])
        original_x, original_u = self.x.copy(), self.u.copy()
        regimes = Counter()
        with patch('B.gate_hydraulics.sluice_gate_flow',
                   wraps=sluice_gate_flow) as actual_gate:
            rhs = self.make_b2_rhs(self.gates, regimes)
            dx, y = rhs(0, self.x, self.u, self.params)
            self.assertEqual(actual_gate.call_count, 2)
            self.assertEqual([call.args[2] for call in actual_gate.call_args_list],
                             self.u[:2].tolist())
        expected_delta = self.np.zeros(8)
        expected_delta[2:4] = expected_q  # Receiving wetlands: +m3/s.
        expected_delta[4:6] = -expected_q  # Donor pools: identical -m3/s.
        self.np.testing.assert_allclose(dx-dx0, expected_delta, rtol=0, atol=1e-14)
        self.np.testing.assert_array_equal(y[[self.a2.Y['Q_G1'], self.a2.Y['Q_G2']]],
                                           expected_q)
        other_outputs = self.np.ones(len(y), dtype=bool)
        other_outputs[[self.a2.Y['Q_G1'], self.a2.Y['Q_G2']]] = False
        self.np.testing.assert_array_equal(y[other_outputs], y0[other_outputs])
        self.assertAlmostEqual(float(sum(dx)), float(sum(dx0)), places=13)
        self.assertAlmostEqual(float(sum(dx)), y[self.a2.Y['external_net']], places=13)
        self.assertEqual(regimes, Counter({(1, 'free'): 1, (2, 'free'): 1}))
        self.np.testing.assert_array_equal(self.x, original_x)
        self.np.testing.assert_array_equal(self.u, original_u)

    def test_replaced_gate_flux_has_no_old_resistance_or_capacity_dependence(self):
        from copy import deepcopy
        rhs = self.make_b2_rhs(self.gates)
        dx, y = rhs(0, self.x, self.u, self.params)
        changed = deepcopy(self.params)
        changed.R_G = self.np.array([9876., .0001])
        changed.Q_Gcap = self.np.array([.00001, 10000.])
        changed_dx, changed_y = rhs(0, self.x, self.u, changed)
        self.np.testing.assert_array_equal(changed_dx, dx)
        self.np.testing.assert_array_equal(changed_y, y)

    def test_closed_gates_leave_both_independent_spillways_active(self):
        spilling = self.x.copy()
        spilling[4:6] = self.params.A_d * (self.params.z_dsp-self.params.z_d
                                         + self.np.array([.2, .3]))
        closed = self.u.copy()
        closed[:2] = 0
        dx0, y0 = self.zero_gate_reference(spilling)
        dx, y = self.make_b2_rhs(self.gates)(0, spilling, closed, self.params)
        self.assertGreater(y[self.a2.Y['Q_dsp1']], 0)
        self.assertGreater(y[self.a2.Y['Q_dsp2']], 0)
        self.assertEqual(y[self.a2.Y['Q_G1']], 0)
        self.assertEqual(y[self.a2.Y['Q_G2']], 0)
        self.np.testing.assert_array_equal(dx, dx0)
        self.np.testing.assert_array_equal(y, y0)

    def test_invalid_intermediate_storage_is_rejected_not_repaired(self):
        rhs = self.make_b2_rhs(self.gates)
        for index, value in ((0, -1), (1, self.params.C_m+1), (6, math.nan)):
            with self.subTest(index=index, value=value):
                bad = self.x.copy()
                bad[index] = value
                with self.assertRaises(FloatingPointError):
                    rhs(123, bad, self.u, self.params)
                if math.isnan(value):
                    self.assertTrue(math.isnan(bad[index]))
                else:
                    self.assertEqual(bad[index], value)

    def test_intermediate_groundwater_applicability_failure_is_preserved(self):
        from copy import deepcopy
        changed = deepcopy(self.params)
        changed.soil_bottom = changed.H_g_init-.01
        self.assertLess(changed.H_g_init, changed.H_g_valid_max)
        with self.assertRaises(self.trajectory_failure):
            self.make_b2_rhs(self.gates)(123, self.x, self.u, changed)

    def test_unsupported_gate_state_has_project_applicability_diagnostic(self):
        from B.optimize_wetland import HydraulicApplicabilityFailure
        # At g1=.4 this larger TEST opening is .4 m, above the .3 m pool depth.
        gates = (replace(self.gates[0], gate_max_opening=1), self.gates[1])
        with self.assertRaises(HydraulicApplicabilityFailure) as raised:
            self.make_b2_rhs(gates)(0, self.x, self.u, self.params)
        self.assertIsInstance(raised.exception.__cause__, UnsupportedGateRegime)
        self.assertEqual(raised.exception.__cause__.regime, 'weir_or_partially_wetted_opening')
        self.assertEqual(raised.exception.diagnostic['branch'],1)
        self.assertEqual(raised.exception.diagnostic['time_s'],0.)
        self.assertAlmostEqual(raised.exception.diagnostic['residual_m'],.2)


if __name__ == "__main__":
    unittest.main()
