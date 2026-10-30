"""Focused rain-history checks; modified fixtures are not fitted observations."""
from __future__ import annotations

from pathlib import Path
import unittest

import numpy as np
import pandas as pd

from C import mimo_history as history


WEATHER = Path(__file__).resolve().parents[1] / "results/mimo_srs/causal_weather_hourly.csv"
OLDER = ["P_completed_1h_at_cutoff_minus_1h", "P_completed_1h_at_cutoff_minus_2h"]


class MimoHistoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not WEATHER.is_file():
            raise unittest.SkipTest("Prepare the separate original FIT-only data first")
        cls.weather = pd.read_csv(WEATHER, float_precision="round_trip", low_memory=False)

    def fixture(self):
        # The initial complete cached rows have ordinary contiguous native supports.
        # Alterations below deliberately challenge eligibility, not model skill.
        return self.weather.iloc[1:14].copy().reset_index(drop=True)

    def reason_at(self, audit, timestamp):
        selected = pd.to_datetime(audit.timestamp_utc, utc=True).eq(pd.Timestamp(timestamp))
        self.assertEqual(int(selected.sum()), 1)
        return audit.loc[selected, "reason"].iloc[0]

    def test_cached_counts_and_actual_supports_match_independent_audit(self):
        pairs, audit = history.build_history_pairs(self.weather)
        periods = history.split_periods(pairs)
        self.assertEqual({name: len(frame) for name, frame in periods.items()}, {
            "fit_2019_2022": 33412, "validation_2023": 8078, "external_2024": 6951,
        })
        cutoff = pd.to_datetime(pairs.prediction_cutoff_utc, utc=True)
        current_end = pd.to_datetime(pairs.P_predictor_end_utc, utc=True)
        current_start = pd.to_datetime(pairs.P_predictor_start_utc, utc=True)
        self.assertTrue(current_end.le(cutoff).all())
        self.assertTrue((cutoff-current_end).eq(pd.Timedelta(minutes=8)).all())
        for name, expected_age in zip(OLDER, [68., 128.]):
            start = pd.to_datetime(pairs[f"{name}_start_utc"], utc=True)
            end = pd.to_datetime(pairs[f"{name}_end_utc"], utc=True)
            self.assertTrue(end.le(cutoff).all())
            self.assertTrue((end-start).eq(pd.Timedelta(hours=1)).all())
            self.assertTrue(current_start.eq(end).all())
            np.testing.assert_array_equal(pairs[f"{name}_age_at_cutoff_minutes"], expected_age)
            current_start = start
        self.assertEqual(self.reason_at(audit, "2019-02-11T01:00Z"),
                         "nonconsecutive_native_rain_intervals")

    def test_missing_older_rain_is_excluded_without_zero_filling(self):
        weather = self.fixture()
        target = weather.loc[5, "timestamp_utc"]
        weather.loc[3, "precipitation_mm_per_h"] = np.nan
        pairs, audit = history.build_history_pairs(weather)
        self.assertEqual(self.reason_at(audit, target), "missing_older_rain_history")
        self.assertFalse(pd.to_datetime(pairs.timestamp_utc, utc=True).eq(pd.Timestamp(target)).any())
        self.assertTrue(pd.isna(weather.loc[3, "precipitation_mm_per_h"]))

    def test_missing_grid_hour_cannot_be_replaced_by_an_earlier_row(self):
        weather = self.fixture()
        target = weather.loc[5, "timestamp_utc"]
        # Target and immediate predecessor remain present, but c-1h is absent.
        weather = weather.drop(index=3).reset_index(drop=True)
        pairs, audit = history.build_history_pairs(weather)
        self.assertEqual(self.reason_at(audit, target), "missing_older_rain_history")
        self.assertFalse(pd.to_datetime(pairs.timestamp_utc, utc=True).eq(pd.Timestamp(target)).any())

    def test_native_gap_and_future_ending_source_are_rejected_without_retiming(self):
        for mode, expected_reason in [
            ("gap", "nonconsecutive_native_rain_intervals"),
            ("future", "invalid_native_rain_support"),
        ]:
            with self.subTest(mode=mode):
                weather = self.fixture()
                target = weather.loc[5, "timestamp_utc"]
                # Row2 is the oldest rain source needed by target row5.
                original_end = pd.Timestamp(weather.loc[2, "precipitation_source_end_utc"])
                end = (original_end-pd.Timedelta(minutes=1) if mode == "gap" else
                       pd.Timestamp(weather.loc[2, "timestamp_utc"])+pd.Timedelta(minutes=1))
                weather.loc[2, "precipitation_source_start_utc"] = str(end-pd.Timedelta(hours=1))
                weather.loc[2, "precipitation_source_end_utc"] = str(end)
                weather.loc[2, "precipitation_source_latest_utc"] = str(end)
                weather.loc[2, "precipitation_source_age_minutes"] = (
                    pd.Timestamp(weather.loc[2, "timestamp_utc"])-end).total_seconds()/60
                pairs, audit = history.build_history_pairs(weather)
                self.assertEqual(self.reason_at(audit, target), expected_reason)
                self.assertFalse(pd.to_datetime(pairs.timestamp_utc, utc=True).eq(pd.Timestamp(target)).any())
                self.assertEqual(pd.Timestamp(weather.loc[2, "precipitation_source_end_utc"]), end)

    def test_finite_extreme_and_negative_temperature_are_retained(self):
        weather = self.fixture()
        target = weather.loc[5, "timestamp_utc"]
        weather.loc[5, "precipitation_mm_per_h"] = 1234.0
        weather.loc[3, "precipitation_mm_per_h"] = 987.0
        weather.loc[5, "temperature_C"] = -25.0
        pairs, audit = history.build_history_pairs(weather)
        selected = pairs.loc[pd.to_datetime(pairs.timestamp_utc, utc=True).eq(pd.Timestamp(target))].iloc[0]
        self.assertEqual(self.reason_at(audit, target), "eligible")
        self.assertEqual(selected.P, 1234.0)
        self.assertEqual(selected[OLDER[0]], 987.0)
        self.assertEqual(selected["T"], -25.0)

    def test_selection_uses_only_validation_rain_rmse(self):
        metrics = pd.DataFrame([
            {"period": "validation_2023", "model": "one_hour", "variable": "P", "RMSE": 1.0},
            {"period": "validation_2023", "model": "three_hour", "variable": "P", "RMSE": 0.9},
            {"period": "validation_2023", "model": "one_hour", "variable": "T", "RMSE": 0.01},
            {"period": "validation_2023", "model": "three_hour", "variable": "T", "RMSE": 100.0},
            {"period": "external_2024", "model": "one_hour", "variable": "P", "RMSE": 0.01},
            {"period": "external_2024", "model": "three_hour", "variable": "P", "RMSE": 100.0},
        ])
        self.assertEqual(history.choose_on_validation(metrics), "three_hour")
        metrics.loc[metrics.period.eq("external_2024"), "RMSE"] = [1e8, 0.0]
        self.assertEqual(history.choose_on_validation(metrics), "three_hour")
        metrics.loc[metrics.period.eq("validation_2023") & metrics.variable.eq("P") &
                    metrics.model.eq("three_hour"), "RMSE"] = 1.1
        self.assertEqual(history.choose_on_validation(metrics), "one_hour")


if __name__ == "__main__":
    unittest.main()
