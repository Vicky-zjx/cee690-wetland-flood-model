"""Analytic stationarity checks and chronological, numerically screened selection."""
import unittest

import numpy as np
import pandas as pd

from C.mimo_penalty import kkt_check, validation_choice


class PenaltyDiagnosticsTests(unittest.TestCase):
    def test_analytic_soft_threshold_solution_and_incorrect_zero(self):
        B = np.eye(3)
        y = np.array([2., -2., .2])
        # Exact minimizer of sum((c-y)^2) + sum(abs(c)).
        correct = kkt_check(B, y, [1.5,-1.5,0.],1.)
        self.assertTrue(correct['kkt_pass'])
        self.assertAlmostEqual(correct['kkt_max_abs_residual'],0.)
        invalid = kkt_check(B,y,[0.,-1.5,0.],1.)
        self.assertFalse(invalid['kkt_pass'])
        self.assertAlmostEqual(invalid['kkt_max_abs_residual'],3.)

    def test_external_scores_and_numerically_failed_candidates_cannot_select(self):
        metrics = pd.DataFrame([
            ['validation_2023','l1_1','P',.5],
            ['validation_2023','l1_0p1','P',.4],
            ['validation_2023','l1_0p01','P',.3],
            ['external_2024','l1_1','P',.01],
            ['external_2024','l1_0p1','P',100.],
        ],columns=['period','model','variable','RMSE'])
        diagnostics = pd.DataFrame([
            ['l1_1','P',True],['l1_0p1','P',True],['l1_0p01','P',False]],
            columns=['model','output','kkt_pass'])
        choice=validation_choice(metrics,diagnostics)
        self.assertEqual(choice['lowest_validation_RMSE_descriptive'],'l1_0p01')
        self.assertEqual(choice['selected_after_numerical_screen'],'l1_0p1')
        diagnostics['kkt_pass']=False
        self.assertIsNone(validation_choice(metrics,diagnostics)['selected_after_numerical_screen'])


if __name__ == '__main__':
    unittest.main()
