"""Focused official FIT checks; synthetic fixtures are not weather observations."""
from __future__ import annotations

import hashlib
import inspect
import json
import os
from pathlib import Path
import unittest

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd

from C import mimo_fit as fit
from C.mimo_figures import _residual_acf


def weather_fixture(timestamps):
    """Small deterministic inputs for alignment tests, not fitting evidence."""
    n = len(timestamps)
    return pd.DataFrame({
        "timestamp_utc": pd.to_datetime(timestamps, utc=True),
        "precipitation_mm_per_h": np.arange(n, dtype=float),
        "temperature_C": np.arange(n, dtype=float) - 12,
        "solar_W_per_m2": np.arange(n, dtype=float) * 20,
        "wind_m_per_s": np.arange(n, dtype=float) / 2,
    })


class MimoFitTests(unittest.TestCase):
    def test_pairs_preserve_extremes_negative_temperature_and_real_gaps(self):
        times = pd.to_datetime([
            "2023-01-01T00:00Z", "2023-01-01T01:00Z", "2023-01-01T02:00Z",
            "2023-01-01T03:00Z", "2023-01-01T05:00Z", "2023-01-01T06:00Z",
            "2023-01-01T07:00Z",
        ], utc=True)
        weather = weather_fixture(times)
        weather.loc[2, "precipitation_mm_per_h"] = np.nan
        weather.loc[5, "precipitation_mm_per_h"] = 1234.0
        pairs, audit = fit.make_pairs(weather.iloc[::-1])
        self.assertEqual(pairs["timestamp_utc"].tolist(), list(times[[1, 5, 6]]))
        self.assertTrue((pairs["T"] < 0).all())
        self.assertEqual(pairs.loc[0, "P_completed_1h_at_cutoff"], 0.0)
        self.assertEqual(pairs.loc[1, "P"], 1234.0)
        self.assertEqual(pairs.loc[2, "P_completed_1h_at_cutoff"], 1234.0)
        self.assertEqual(audit.loc[2, "reason"], "missing_or_nonfinite_response")
        self.assertEqual(audit.loc[3, "reason"], "missing_or_nonfinite_lag")
        self.assertEqual(audit.loc[4, "reason"], "no_exact_previous_hour")
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            fit.make_pairs(pd.concat([weather, weather.iloc[[0]]]))

    def test_2024_observations_never_change_development_or_scaling(self):
        weather = weather_fixture(pd.date_range("2023-12-31T20:00Z", periods=10, freq="h"))
        changed = weather.copy()
        external = changed["timestamp_utc"].dt.year.eq(2024)
        changed.loc[external, fit.COLUMNS] += 100000.0
        development, evaluation = fit.split_periods(fit.make_pairs(weather)[0])
        development_changed, evaluation_changed = fit.split_periods(fit.make_pairs(changed)[0])
        pd.testing.assert_frame_equal(development, development_changed)
        self.assertFalse(evaluation[fit.OUTPUTS].equals(evaluation_changed[fit.OUTPUTS]))
        for columns in (fit.PREDICTORS, fit.OUTPUTS):
            before = development[columns].to_numpy().T
            after = development_changed[columns].to_numpy().T
            np.testing.assert_array_equal(before, after)
            scale_before, scale_after = fit.fit_affine(before), fit.fit_affine(after)
            for key in scale_before:
                np.testing.assert_array_equal(scale_before[key], scale_after[key])
        self.assertEqual(evaluation.iloc[0]["timestamp_utc"], pd.Timestamp("2024-01-01T00:00Z"))
        self.assertEqual(evaluation.iloc[0]["predictor_timestamp_utc"], pd.Timestamp("2023-12-31T23:00Z"))

    def test_affine_transform_is_reversible_and_heldout_values_are_not_clipped(self):
        development = np.array([[-20., 0., 20.], [0., 10., 100.], [3., 3., 3.]])
        scaling = fit.fit_affine(development, bound=0.9)
        scaled = fit.transform(development, scaling)
        np.testing.assert_allclose(fit.transform(scaled, scaling, inverse=True), development,
                                   rtol=0, atol=1e-13)
        self.assertLessEqual(np.max(abs(scaled)), 0.9 + 1e-15)
        heldout = np.array([[-100., 100.], [-5., 200.], [0., 6.]])
        transformed = fit.transform(heldout, scaling)
        self.assertGreater(np.max(abs(transformed)), 1.0)
        np.testing.assert_allclose(fit.transform(transformed, scaling, inverse=True), heldout,
                                   rtol=0, atol=1e-13)

    def test_calendar_uses_leap_year_and_civil_time_without_breaking_utc_adjacency(self):
        leap_noon = pd.DatetimeIndex(["2024-02-29T18:00Z"])
        annual_phase = 2*np.pi*59.5/366
        np.testing.assert_allclose(fit.calendar_features(leap_noon)[0],
                                   [np.sin(annual_phase), np.cos(annual_phase), 0., -1.], atol=1e-14)
        spring = pd.date_range("2024-03-10T07:00Z", periods=2, freq="h")
        autumn = pd.date_range("2024-11-03T06:00Z", periods=2, freq="h")
        for times in (spring, autumn):
            pairs, _ = fit.make_pairs(weather_fixture(times))
            self.assertEqual(len(pairs), 1)
            self.assertEqual(pairs.loc[0, "timestamp_utc"] - pairs.loc[0, "predictor_timestamp_utc"],
                             pd.Timedelta(hours=1))
        spring_features = fit.calendar_features(spring)
        np.testing.assert_allclose(spring_features[:, 2], np.sin(2*np.pi*np.array([1., 3.])/24))
        autumn_features = fit.calendar_features(autumn)
        np.testing.assert_array_equal(autumn_features[0], autumn_features[1])

    def test_official_import_has_verified_source_path_and_fingerprint(self):
        module, provenance = fit.official_module()
        expected = fit.REPO / "A/work/multivarious/multivarious/fit/mimo_srs.py"
        self.assertEqual(Path(inspect.getsourcefile(module.mimo_srs)).resolve(), expected.resolve())
        self.assertEqual(Path(provenance["source_path"]).resolve(), expected.resolve())
        self.assertEqual(provenance["source_sha256"], hashlib.sha256(expected.read_bytes()).hexdigest())
        self.assertEqual(provenance["source_sha256"],
                         "36c56bf25006bcc0b9f3af1092c3c0eff07fce150f589341443ab5b1e8ccf38b")
        self.assertEqual(provenance["checkout_revision"], "06d7aed5b541ec7fe26a6dc69be6aa7fb610aa02")

    def test_causal_source_selection_is_backward_and_age_limited(self):
        from C.mimo_causal_data import _asof_values
        targets=pd.DatetimeIndex(['2024-01-01T01:00Z','2024-01-01T03:00Z'])
        times=pd.DatetimeIndex(['2024-01-01T00:55Z','2024-01-01T01:01Z'])
        values,source,age=_asof_values(targets,times,[2.,999.])
        self.assertEqual(values[0],2.)
        self.assertEqual(age[0],5.)
        self.assertTrue(np.isnan(values[1]))
        self.assertTrue((source<=targets).all())

    def test_real_causal_artifact_preserves_solar_and_has_no_future_source(self):
        path=fit.ROOT/'results/mimo_srs/causal_weather_hourly.csv'
        if not path.exists():
            self.skipTest('Run the requested causal-data preparation first')
        frame=pd.read_csv(path,float_precision='round_trip',low_memory=False)
        times=pd.to_datetime(frame.timestamp_utc,utc=True)
        valid=np.isfinite(frame[fit.COLUMNS]).all(axis=1)
        for name in ['precipitation','temperature','wind','solar']:
            sources=pd.to_datetime(frame[f'{name}_source_latest_utc'],utc=True)
            self.assertTrue(sources[valid].le(times[valid]).all())
        for name in ['temperature','wind']:
            self.assertTrue(frame.loc[valid,f'{name}_source_age_minutes'].between(0,60).all())
        self.assertTrue(frame.loc[valid,'precipitation_source_age_minutes'].between(0,60,inclusive='left').all())
        native_start=pd.to_datetime(frame.precipitation_source_start_utc,utc=True)
        native_end=pd.to_datetime(frame.precipitation_source_end_utc,utc=True)
        self.assertTrue((native_end[valid]-native_start[valid]).eq(pd.Timedelta(hours=1)).all())
        self.assertFalse(native_end[valid].duplicated().any())
        pvalid=frame.precipitation_mm_per_h.notna()
        intervals=pd.DataFrame({'start':native_start[pvalid],'end':native_end[pvalid]}).reset_index(drop=True)
        self.assertTrue(intervals.start.iloc[1:].ge(intervals.end.shift(1).iloc[1:]).all())
        original=pd.read_csv(fit.ROOT/'data/processed/weather_hourly.csv',float_precision='round_trip')
        original.index=pd.to_datetime(original.timestamp_utc,utc=True)
        expected=original.solar_W_per_m2.reindex(pd.DatetimeIndex(times)-pd.Timedelta(hours=1)).to_numpy()
        np.testing.assert_array_equal(frame.solar_W_per_m2.to_numpy(),expected)
        pairs,_=fit.make_pairs(frame)
        sources=pd.Series(pd.to_datetime(frame.source_latest_utc,utc=True).to_numpy(),index=times)
        lag_sources=pd.to_datetime(pairs.predictor_timestamp_utc.map(sources),utc=True)
        self.assertTrue(lag_sources.le(pairs.predictor_timestamp_utc).all())

    def test_pair_support_records_actual_intervals_and_rejects_future_predictors(self):
        path=fit.ROOT/'results/mimo_srs/causal_weather_hourly.csv'
        if not path.exists():
            self.skipTest('Prepare the separate FIT-only data first')
        weather=pd.read_csv(path,float_precision='round_trip',low_memory=False).iloc[:6].copy()
        pairs,_=fit.make_pairs(weather)
        aligned,examples,summary=fit.attach_pair_support(pairs,weather)
        self.assertTrue((aligned.prediction_cutoff_utc==aligned.predictor_timestamp_utc).all())
        self.assertTrue((aligned.P_predictor_end_utc<=aligned.prediction_cutoff_utc).all())
        # Actual native targets may begin before the cutoff: record, do not relabel them.
        self.assertTrue((aligned.P_target_start_utc<aligned.prediction_cutoff_utc).all())
        self.assertTrue((aligned.P_target_end_utc>aligned.prediction_cutoff_utc).all())
        self.assertFalse(summary['exact_future_clock_hour_rainfall_available'])
        self.assertTrue(aligned.Rs_support_type.str.contains('assumed').all())
        cutoff=aligned.prediction_cutoff_utc.iloc[0]
        bad=weather.copy()
        index=pd.to_datetime(bad.timestamp_utc,utc=True).eq(cutoff)
        bad.loc[index,'precipitation_source_end_utc']=str(cutoff+pd.Timedelta(minutes=1))
        with self.assertRaisesRegex(ValueError,'ends after its cutoff'):
            fit.attach_pair_support(pairs,bad)

    def test_saved_official_model_reproduces_external_predictions(self):
        results = fit.ROOT / "results/mimo_srs"
        paths = [results / "fit/official_model.npz", results / "scaling.json",
                 results / "external_predictions.csv"]
        if not all(path.is_file() for path in paths):
            self.skipTest("Run python3 -B -m C.run_fitting mimo to create real fitted artifacts first")
        frame = pd.read_csv(paths[2], float_precision="round_trip")
        with np.load(paths[0]) as saved:
            model = {key: saved[key] for key in fit.MODEL_KEYS}
        scaling = json.loads(paths[1].read_text())
        config = json.loads((fit.ROOT / "configs/fitting.json").read_text())
        settings = dict(fit.DEFAULT_SETTINGS, **config.get("mimo", {}))
        predictions = fit.predict(model, frame[fit.PREDICTORS].to_numpy().T, scaling, settings).T
        np.testing.assert_array_equal(predictions, frame[[f"{name}_predicted" for name in fit.OUTPUTS]].to_numpy())
        times = pd.to_datetime(frame["timestamp_utc"], utc=True)
        lagged = pd.to_datetime(frame["predictor_timestamp_utc"], utc=True)
        self.assertTrue(times.dt.year.eq(2024).all())
        self.assertTrue((times-lagged).eq(pd.Timedelta(hours=1)).all())
        execution = json.loads((results / "fit/execution.json").read_text())
        self.assertEqual(execution["status"], "completed")
        self.assertEqual(execution["supplied_samples"], execution["retained_samples"])
        self.assertTrue(execution["exact_internal_split_verified"])

    def test_r_squared_uses_observed_total_sum_of_squares(self):
        result = fit.metrics(np.array([0., 1., 2.]), np.array([0., 2., 4.]))
        self.assertEqual(result["N"], 3)
        self.assertAlmostEqual(result["MAE"], 1.)
        self.assertAlmostEqual(result["mean_prediction_error"], 1.)
        self.assertAlmostEqual(result["RMSE"], np.sqrt(5/3))
        self.assertAlmostEqual(result["R_squared"], -1.5)

    def test_residual_acf_counts_exact_hours_without_compressing_gaps(self):
        frame = pd.DataFrame({"timestamp_utc": pd.to_datetime([
            "2024-01-01T00:00Z", "2024-01-01T01:00Z",
            "2024-01-01T04:00Z", "2024-01-01T05:00Z",
        ], utc=True)})
        for name in fit.OUTPUTS:
            frame[f"{name}_residual"] = [1., 2., 9., -3.]
        acf = _residual_acf(frame, max_lag=3)
        for name in fit.OUTPUTS:
            component = acf[acf["variable"].eq(name)].set_index("lag_hours")
            self.assertEqual(component["pair_count"].tolist(), [2, 0, 1])
            self.assertAlmostEqual(component.loc[1, "autocorrelation"], -1.)
            self.assertTrue(np.isnan(component.loc[2, "autocorrelation"]))
            self.assertTrue(np.isnan(component.loc[3, "autocorrelation"]))


if __name__ == "__main__":
    unittest.main()
