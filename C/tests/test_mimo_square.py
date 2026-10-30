"""Focused checks for the single approved physical rainfall-square feature."""
from __future__ import annotations

import hashlib
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

from C import mimo_fit as base
from C import mimo_square as square


SOURCE = Path(__file__).resolve().parents[1] / "results/mimo_rain_history"
RECENT = "P_completed_1h_at_cutoff"
SQUARE = "P_completed_1h_squared"
COUNTS = {"fit_2019_2022": 33412, "validation_2023": 8078, "external_2024": 6951}


class MimoSquareTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.paths = {name: SOURCE/f"{name}_pairs.csv" for name in COUNTS}
        if not all(path.is_file() for path in cls.paths.values()):
            raise unittest.SkipTest("Run the approved rain-history comparison first")
        cls.hashes_before = {name: hashlib.sha256(path.read_bytes()).hexdigest()
                             for name, path in cls.paths.items()}
        cls.periods = square.load_comparison_pairs()

    def test_loader_preserves_exact_cached_rows_values_and_source_files(self):
        self.assertEqual({name: len(frame) for name, frame in self.periods.items()}, COUNTS)
        for name, path in self.paths.items():
            with self.subTest(period=name):
                original = pd.read_csv(path, float_precision="round_trip", low_memory=False)
                loaded = self.periods[name]
                for column in original:
                    if column.endswith("_utc"):
                        original[column] = pd.to_datetime(original[column], utc=True)
                pd.testing.assert_frame_equal(loaded[list(original)], original)
                self.assertEqual(str(loaded.timestamp_utc.dt.tz), "UTC")
                np.testing.assert_array_equal(loaded[SQUARE], loaded[RECENT]**2)
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), self.hashes_before[name])

    def test_square_uses_physical_predictor_only_without_mutating_input_or_targets(self):
        pairs = pd.DataFrame({RECENT: [0., 0.5, 1234.], "P": [100., 5., 0.],
                              "T": [-20., 0., 25.]})
        original = pairs.copy(deep=True)
        enriched = square.add_square(pairs)
        np.testing.assert_array_equal(enriched[SQUARE], [0., 0.25, 1522756.])
        pd.testing.assert_frame_equal(pairs, original)
        pd.testing.assert_frame_equal(enriched[list(original)], original)
        altered_targets = pairs.copy()
        altered_targets["P"] = [0., 1234., 98765.]
        np.testing.assert_array_equal(square.add_square(altered_targets)[SQUARE], enriched[SQUARE])

    def test_invalid_precipitation_is_rejected_without_clipping_or_filling(self):
        for value in [-0.1, np.nan, np.inf, -np.inf]:
            with self.subTest(value=value):
                frame = pd.DataFrame({RECENT: [value], "P": [0.]})
                before = frame.copy(deep=True)
                with self.assertRaises(ValueError):
                    square.add_square(frame)
                pd.testing.assert_frame_equal(frame, before)

    def test_selection_uses_validation_rain_rmse_and_ties_keep_one_hour(self):
        metrics = pd.DataFrame([
            ("validation_2023", "one_hour", "P", 1.0),
            ("validation_2023", "rain_square", "P", 0.9),
            ("validation_2023", "one_hour", "T", 0.01),
            ("validation_2023", "rain_square", "T", 100.0),
            ("external_2024", "one_hour", "P", 0.01),
            ("external_2024", "rain_square", "P", 100.0),
        ], columns=["period", "model", "variable", "RMSE"])
        self.assertEqual(square.choose_on_validation(metrics), "rain_square")
        metrics.loc[metrics.period.eq("external_2024"), "RMSE"] = [1e8, 0.0]
        self.assertEqual(square.choose_on_validation(metrics), "rain_square")
        mask = (metrics.period.eq("validation_2023") & metrics.variable.eq("P") &
                metrics.model.eq("rain_square"))
        for value in [1.0, 1.1]:
            metrics.loc[mask, "RMSE"] = value
            self.assertEqual(square.choose_on_validation(metrics), "one_hour")

    def test_external_changes_do_not_affect_training_values_or_affine_scaling(self):
        features = base.PREDICTORS + [SQUARE]
        train = self.periods["fit_2019_2022"]
        original = train[features].to_numpy().T.copy()
        scale_before = base.fit_affine(original)
        external = self.periods["external_2024"].copy(deep=True)
        external[RECENT] += 10000.
        external[base.OUTPUTS] += 50000.
        external = square.add_square(external)
        self.assertGreater(external[SQUARE].min(), train[SQUARE].max())
        np.testing.assert_array_equal(train[features].to_numpy().T, original)
        scale_after = base.fit_affine(train[features].to_numpy().T)
        for name in scale_before:
            np.testing.assert_array_equal(scale_before[name], scale_after[name])


if __name__ == "__main__":
    unittest.main()
