"""Focused contracts for fixed-gate wetland design search; no project ODE runs."""
from __future__ import annotations

import copy
import unittest
from unittest.mock import patch
import numpy as np
from B import optimize_wetland as old
from B import no_spillway_model as model
from B import wetland_only_study as study


class WetlandOnlyTests(unittest.TestCase):
    def setUp(self):
        self.config = old.load_config()
        self.gates = np.linspace(.02, .05, 24)
        self.gates.setflags(write=False)
        self.covers = np.array([.3, .8])

    def fixture(self):
        """Artificial complete-audit arrays, never observations or study results."""
        params = model.base_parameters(self.config)
        time = np.array([0., 120.])
        x = np.repeat(params.x_init[:, None], 2, axis=1)
        y = np.zeros((len(old.a2.OUTPUT_NAMES), 2))
        for label, level in [('H_d1', 100.3), ('H_d2', 100.1),
                             ('H_w1', 98.25), ('H_w2', 98.15),
                             ('H_g', 97.2), ('H_r', 97.35)]:
            y[old.a2.Y[label]] = level
        u = np.zeros((len(old.a2.INPUT_NAMES), 2))
        u[:2] = [[.03, .04], [.04, .05]]
        theta = study.assemble_design(self.covers, self.gates)
        residuals, audit = model.constraint_audit(time, x, y, u, theta, params, self.config)
        return dict(theta=theta, audit=audit, residuals=residuals, pair=['type_a', 'type_b'],
                    metrics={'peak_stage_m':97.35}, time=time, x=x, y=y, u=u)

    def test_only_two_covers_change_and_all_twenty_four_gate_knots_are_exact(self):
        gates_before = self.gates.copy()
        covers_before = self.covers.copy()
        for covers in (self.covers, np.array([0., 1.]), np.array([1., 0.])):
            full = study.assemble_design(covers, self.gates)
            self.assertEqual(full.shape, (26,))
            np.testing.assert_array_equal(full[:-2], gates_before)
            np.testing.assert_array_equal(full[-2:], covers)
            full[0] = .9  # New output does not expose mutable frozen gate storage.
            np.testing.assert_array_equal(self.gates, gates_before)
        np.testing.assert_array_equal(self.covers, covers_before)
        self.assertFalse(self.gates.flags.writeable)

    def test_extra_gate_decisions_and_invalid_covers_are_rejected_without_clipping(self):
        for covers in (np.zeros(26), np.zeros(1), np.zeros((1, 2)), [np.nan, .5], [np.inf, .5]):
            with self.assertRaises(ValueError):
                study.assemble_design(covers, self.gates)
        for covers in (np.array([-.001, .5]), np.array([.5, 1.001])):
            before = covers.copy()
            with self.assertRaises(old.TrajectoryFailure):
                study.assemble_design(covers, self.gates)
            np.testing.assert_array_equal(covers, before)
        with self.assertRaises(ValueError):
            study.assemble_design(self.covers, np.zeros(23))

    def test_optimizer_has_eighteen_residuals_and_g5_is_separate_validation(self):
        run = self.fixture()
        original = copy.deepcopy(run['audit'])
        view = study.optimization_view(run)
        self.assertEqual(view['residuals'].shape, (18,))
        self.assertEqual(len(view['active_audit']), 18)
        self.assertEqual({r['label'] for r in view['applicability_audit']}, study.G5_LABELS)
        self.assertFalse(study.G5_LABELS & {r['label'] for r in view['active_audit']})
        self.assertEqual([r['label'] for r in view['active_audit']],
                         [r['label'] for r in run['audit'][:18]])
        self.assertEqual(view['active_audit'][-1]['label'], 'integration_status')
        np.testing.assert_array_equal(view['residuals'], run['residuals'][:18])
        self.assertEqual(len(view['fixed_gate_bound_audit']), 4)
        self.assertEqual(len(view['native_cover_audit']), 4)
        self.assertTrue(view['feasible'])
        self.assertEqual(run['audit'], original)

    def test_g5_violation_is_model_rejection_even_within_optimizer_tolerance(self):
        run = self.fixture()
        row = next(r for r in run['audit'] if r['label'] == 'gate1_head_opening')
        row.update(scaled_residual=5e-7, strict_pass=False, pass_at_tolerance=True)
        with self.assertRaisesRegex(old.TrajectoryFailure, 'model applicability'):
            study.optimization_view(run)

    def test_safety_violation_stays_an_optimizer_constraint_not_a_model_failure(self):
        run = self.fixture()
        row = next(r for r in run['audit'] if r['label'] == 'H_d1_safe')
        row.update(scaled_residual=.2, margin=-.2, strict_pass=False, pass_at_tolerance=False)
        view = study.optimization_view(run)
        self.assertFalse(view['feasible'])
        self.assertEqual(view['residuals'][9], .2)
        self.assertTrue(all(r['strict_pass'] for r in view['applicability_audit']))

    def test_invalid_evaluation_uses_only_generic_integration_failure_sentinel(self):
        before = copy.deepcopy(self.config)
        f, g = study.invalid_response(self.config)
        self.assertEqual(f, self.config['numerics']['failure_objective'])
        self.assertEqual(g.shape, (18,))
        np.testing.assert_array_equal(g[:-1], np.zeros(17))
        self.assertEqual(g[-1], self.config['numerics']['failure_residual'])
        self.assertGreater(g[-1], self.config['numerics']['feasibility_tol'])
        g[0] = 999.
        _, again = study.invalid_response(self.config)
        self.assertEqual(again[0], 0.)
        self.assertEqual(self.config, before)

    def test_simulation_passes_frozen_gates_categories_forcing_and_refinement_unchanged(self):
        pair = ('type_b', 'type_a')
        forcing = (np.array([0., 120.]), np.zeros((8, 2)))
        stub = self.fixture()
        with patch.object(model, 'simulate', return_value=stub) as simulate:
            result, view = study.simulate_covers(self.covers, pair, self.config,
                self.gates, forcing=forcing, refine=2)
        args, kwargs = simulate.call_args
        np.testing.assert_array_equal(args[0][:-2], self.gates)
        np.testing.assert_array_equal(args[0][-2:], self.covers)
        self.assertEqual(args[1], pair)
        self.assertIs(args[2], self.config)
        self.assertIs(kwargs['forcing'], forcing)
        self.assertEqual(kwargs['refine'], 2)
        self.assertIs(result, stub)
        self.assertTrue(view['feasible'])

    def test_simulation_preserves_existing_gate_guard_exception_and_diagnostic(self):
        diagnostic = dict(branch=2, time_s=123., water_depth_m=.1,
                          opening_height_m=.15, residual_m=.0875)
        failure = model.HydraulicApplicabilityFailure(diagnostic)
        with patch.object(model, 'simulate', side_effect=failure):
            with self.assertRaises(model.HydraulicApplicabilityFailure) as raised:
                study.simulate_covers(self.covers, ('type_a', 'type_b'), self.config, self.gates)
        self.assertIs(raised.exception, failure)
        self.assertEqual(raised.exception.diagnostic, diagnostic)

    def test_returned_trajectory_cannot_silently_change_fixed_gate_schedule(self):
        stub = self.fixture()
        stub['theta'][0] += .001
        with patch.object(model, 'simulate', return_value=stub):
            with self.assertRaises(AssertionError):
                study.simulate_covers(self.covers, ('type_a', 'type_b'), self.config, self.gates)


    def test_analysis_returns_scaled_physical_objective_and_independent_ordered_residuals(self):
        # Fixture-only interface contract; no ODE or study observation is produced.
        stub = self.fixture()
        constants = study.AnalysisConstants(self.config, self.gates, ('type_a', 'type_b'))
        with patch.object(model, 'simulate', return_value=stub):
            f, g = study.analysis(self.covers, constants)
        numerics = self.config['numerics']
        self.assertIsInstance(f, float)
        self.assertEqual(f, (stub['metrics']['peak_stage_m'] - numerics['objective_origin_m'])
                         / numerics['height_scale_m'])
        np.testing.assert_array_equal(g, stub['residuals'][:18])
        evaluation = constants.last_evaluation
        self.assertIs(evaluation['run'], stub)
        self.assertTrue(evaluation['valid'])
        self.assertTrue(evaluation['feasible'])
        self.assertEqual(evaluation['physical_peak_m'], stub['metrics']['peak_stage_m'])
        self.assertIsNone(evaluation['failure'])
        self.assertIsNone(evaluation['diagnostic'])
        self.assertFalse(study.G5_LABELS & {r['label'] for r in evaluation['view']['active_audit']})
        before = evaluation['view']['residuals'].copy()
        g[:] = 999.
        np.testing.assert_array_equal(evaluation['view']['residuals'], before)

    def test_analysis_keeps_valid_unsafe_trajectory_and_physical_objective(self):
        stub = self.fixture()
        row = next(r for r in stub['audit'] if r['label'] == 'H_d1_safe')
        row.update(scaled_residual=.2, margin=-.2, strict_pass=False, pass_at_tolerance=False)
        constants = study.AnalysisConstants(self.config, self.gates, ('type_a', 'type_b'))
        with patch.object(model, 'simulate', return_value=stub):
            f, g = study.analysis(self.covers, constants)
        numerics = self.config['numerics']
        self.assertEqual(f, (stub['metrics']['peak_stage_m'] - numerics['objective_origin_m'])
                         / numerics['height_scale_m'])
        self.assertEqual(g[9], .2)
        self.assertEqual(g[-1], -1.)
        self.assertTrue(constants.last_evaluation['valid'])
        self.assertFalse(constants.last_evaluation['feasible'])
        self.assertIs(constants.last_evaluation['run'], stub)
        self.assertIsNone(constants.last_evaluation['failure'])

    def test_analysis_hydraulic_rejection_replaces_previous_results_and_retains_diagnostic(self):
        constants = study.AnalysisConstants(self.config, self.gates, ('type_a', 'type_b'))
        diagnostic = dict(branch=2, time_s=123., water_depth_m=.1,
                          opening_height_m=.15, residual_m=.0875)
        failure = model.HydraulicApplicabilityFailure(diagnostic)
        with patch.object(model, 'simulate', side_effect=[self.fixture(), failure]):
            study.analysis(self.covers, constants)
            self.assertTrue(constants.last_evaluation['valid'])
            f, g = study.analysis(self.covers, constants)
        self.assertEqual(f, self.config['numerics']['failure_objective'])
        np.testing.assert_array_equal(g[:-1], np.zeros(17))
        self.assertEqual(g[-1], self.config['numerics']['failure_residual'])
        evaluation = constants.last_evaluation
        self.assertIsNone(evaluation['run'])
        self.assertIsNone(evaluation['view'])
        self.assertIsNone(evaluation['physical_peak_m'])
        self.assertFalse(evaluation['valid'])
        self.assertFalse(evaluation['feasible'])
        self.assertEqual(evaluation['failure'], type(failure).__name__ + ': ' + str(failure))
        self.assertEqual(evaluation['diagnostic'], diagnostic)

    def test_analysis_unexpected_errors_propagate_and_clear_previous_context(self):
        constants = study.AnalysisConstants(self.config, self.gates, ('type_a', 'type_b'))
        for failure in (ValueError('invalid configuration'), RuntimeError('programming error')):
            with self.subTest(error=type(failure).__name__):
                with patch.object(model, 'simulate', side_effect=[self.fixture(), failure]):
                    study.analysis(self.covers, constants)
                    self.assertTrue(constants.last_evaluation['valid'])
                    with self.assertRaises(type(failure)) as raised:
                        study.analysis(self.covers, constants)
                self.assertIs(raised.exception, failure)
                self.assertIsNone(constants.last_evaluation)


if __name__ == '__main__':
    unittest.main()
