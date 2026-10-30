"""Synthetic fixtures test implementation, never substitute for historical data."""
import json
import unittest

import numpy as np
import pandas as pd

from C.uncertainty.marginals import fit_marginal, randomized_pit, inverse_marginal
from C.uncertainty.weather_generator import COLUMNS, fit_weather, generate_weather, solar_geometry
from C.uncertainty.validation import precipitation_events, validate_weather


def synthetic_training():
    rng = np.random.default_rng(71)
    index = pd.date_range("2020-01-01", periods=24 * 366, freq="h", tz="UTC")
    scale, daylight = solar_geometry(index, 43.0547222, -88.0461111)
    latent = np.zeros((len(index), 4))
    for i in range(1, len(index)):
        latent[i] = .65 * latent[i - 1] + rng.normal(0, .5, 4)
    rain = np.where(latent[:, 0] > .3, rng.gamma(.7, 2, len(index)), 0)
    wind = np.maximum(0, 3 + latent[:, 3])
    wind[rng.random(len(index)) < .02] = 0
    temperature = 8 + 15 * np.sin(2 * np.pi * np.arange(len(index)) / (24 * 366)) + latent[:, 1]
    solar = np.where(daylight, scale * rng.gamma(4, 160, len(index)), 0)
    return pd.DataFrame({"timestamp_utc": index, COLUMNS[0]: rain, COLUMNS[1]: temperature,
                         COLUMNS[2]: solar, COLUMNS[3]: wind})


class WeatherTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = {"weather": {"pit_seed": 123, "block_hours": 48, "burn_in_hours": 336}, "validation": {}}
        cls.training = synthetic_training()
        cls.fit = fit_weather(cls.training, cls.config)
        cls.target = pd.date_range("2022-09-09", periods=144, freq="h", tz="UTC")

    def test_serialization_and_stable_covariance(self):
        json.dumps(self.fit, allow_nan=False)
        self.assertLess(self.fit["var"]["spectral_radius"], 1)
        self.assertGreater(np.min(np.linalg.eigvalsh(self.fit["var"]["Q"])), 0)
        a, q = np.array(self.fit["var"]["A"]), np.array(self.fit["var"]["Q"])
        stationary = np.array(self.fit["var"]["stationary_covariance"])
        np.testing.assert_allclose(stationary, a @ stationary @ a.T + q, atol=1e-10)
        np.testing.assert_allclose(np.diag(self.fit["dependence"]["adjusted"]), 1)

    def test_reproducible_independent_scenario_streams(self):
        def scenario(i):
            rng = np.random.default_rng(np.random.SeedSequence([690, i]))
            return generate_weather(self.fit, self.target, rng=rng)
        first, second = scenario(0), scenario(1)
        pd.testing.assert_frame_equal(first, scenario(0))
        self.assertFalse(first[COLUMNS].equals(second[COLUMNS]))
        # Requesting scenarios in another order cannot change scenario identity.
        pd.testing.assert_frame_equal(second, scenario(1))

    def test_physical_and_night_constraints_both_methods(self):
        for method in ["copula_var", "block_bootstrap"]:
            out = generate_weather(self.fit, self.target, seed=9, method=method, training_df=self.training)
            self.assertTrue(np.all(out[[COLUMNS[0], COLUMNS[2], COLUMNS[3]]].to_numpy() >= 0))
            _, daylight = solar_geometry(self.target, 43.0547222, -88.0461111)
            self.assertTrue(np.all(out.loc[~daylight, COLUMNS[2]].to_numpy() == 0))

    def test_no_per_scenario_normalization(self):
        out = generate_weather(self.fit, self.target, seed=9)
        long = generate_weather(self.fit, pd.date_range(self.target[0], periods=288, freq="h"), seed=9)
        pd.testing.assert_frame_equal(out, long.iloc[:144].reset_index(drop=True))

    def test_missing_hours_not_var_pairs(self):
        gapped = self.training.drop(index=list(range(100, 124)) + list(range(800, 804))).reset_index(drop=True)
        fitted = fit_weather(gapped, self.config)
        expected = len(gapped) - 1 - 2
        self.assertEqual(fitted["var"]["lag_pair_counts"][0][0], expected)
        # Solar has no latent observed value during fully dark hours.
        self.assertLess(fitted["var"]["lag_pair_counts"][2][2], expected)

    def test_missing_values_not_bootstrapped(self):
        broken = self.training.copy()
        broken[COLUMNS[0]] = np.nan
        with self.assertRaisesRegex(ValueError, "No complete"):
            generate_weather(self.fit, self.target, seed=1, method="block_bootstrap", training_df=broken)

    def test_zero_variance_solar_rejected(self):
        constant = self.training.copy()
        constant[COLUMNS[2]] = 0
        with self.assertRaisesRegex(ValueError, "Zero-variance"):
            fit_weather(constant, self.config)

    def test_randomized_mass_and_ties(self):
        marginal = fit_marginal(np.r_[np.zeros(80), np.ones(20)], ["empirical"], mixed_zero=True)
        p = randomized_pit(marginal, np.array([0.] * 30 + [1.] * 10), np.random.default_rng(2))
        self.assertTrue(np.all(p[:30] < .8))
        self.assertTrue(np.all(p[30:] >= .8))
        self.assertGreater(np.std(p[:30]), 0)
        values = inverse_marginal(marginal, np.array([.01, .79, .81]))
        np.testing.assert_equal(values, [0, 0, 1])

    def test_spells_and_storms_break_at_missing_hour(self):
        frame = self.training.iloc[:12].copy()
        frame[COLUMNS[0]] = [1, 2, np.nan, 3, 0, 0, 0, 0, 0, 0, 4, 0]
        events = precipitation_events(frame)
        self.assertEqual(events["wet_spells_hours"], [2, 1, 1])
        self.assertEqual([s["total_mm"] for s in events["storms"]], [3., 3., 4.])
        self.assertEqual([s["duration_hours"] for s in events["storms"]], [2, 1, 1])

    def test_validation_uses_matched_calendar_support(self):
        observed = self.training.iloc[:144].copy()
        generated = observed.copy()
        observed.loc[20:40, COLUMNS[0]] = np.nan
        generated.loc[20:40, COLUMNS[0]] = 1000.
        result = validate_weather(observed, generated, self.config)
        self.assertEqual(result["metrics"]["observed"][COLUMNS[0]], result["metrics"]["generated"][COLUMNS[0]])
        self.assertEqual(result["matched_support"][COLUMNS[0]]["generated_hours_masked_for_comparison"], 21)
        self.assertEqual(result["generated_full_sample_statistics"][COLUMNS[0]]["count"], 144)
        self.assertEqual(result["checks"][-1]["details"]["nonfinite_values"], 0)

    def test_validation_reports_failures_without_tuning(self):
        frame = self.training.iloc[:144].copy()
        altered = frame.copy()
        altered[COLUMNS[0]] = 100
        result = validate_weather(frame, altered, self.config)
        self.assertFalse(result["passed"])
        self.assertIn("wet_hour_frequency", result["failed_checks"])
        json.dumps(result, allow_nan=False)


if __name__ == "__main__":
    unittest.main()
