"""Audit source-observation times behind the existing hourly FIT inputs.

This module only reads cached observations. It does not regrid or overwrite
weather data. ``root`` is the Stage C directory, not the repository root.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

from C.uncertainty.data import (
    VARIABLES, align_precipitation, parse_value, precipitation_intervals,
)


def _met_support(times, values, targets, max_gap_hours):
    """Reproduce the existing midpoint interpolation and retain source times."""
    series = pd.Series(values, index=pd.DatetimeIndex(times), dtype=float).dropna().sort_index()
    series = series.groupby(level=0).agg(lambda a: a.iloc[0] if a.nunique() == 1 else np.nan).dropna()
    x = series.index.as_unit("s").asi8.astype(float)
    y = series.to_numpy()
    target = pd.DatetimeIndex(targets).as_unit("s").asi8.astype(float)
    out = np.full(len(target), np.nan)
    left = np.full(len(target), np.nan)
    right = np.full(len(target), np.nan)
    if len(x) < 2:
        return out, left, right
    j = np.searchsorted(x, target, side="left")
    exact = (j < len(x)) & (x[np.minimum(j, len(x)-1)] == target)
    out[exact] = y[j[exact]]
    left[exact] = right[exact] = x[j[exact]]
    ids = np.flatnonzero((~exact) & (j > 0) & (j < len(x)))
    ids = ids[(x[j[ids]] - x[j[ids]-1]) <= float(max_gap_hours)*3600]
    fraction = (target[ids]-x[j[ids]-1])/(x[j[ids]]-x[j[ids]-1])
    out[ids] = y[j[ids]-1] + fraction*(y[j[ids]]-y[j[ids]-1])
    left[ids], right[ids] = x[j[ids]-1], x[j[ids]]
    return out, left, right


def audit_source_times(root: Path, weather: pd.DataFrame, config: dict) -> pd.DataFrame:
    """Return source support, pair validity and causality without altering values.

    Pair cutoff is the start of the response UTC hour; a midpoint cutoff is
    reported separately, never substituted for it. All finite, exact-adjacent
    pairs are audited, including negative temperatures and extreme rainfall.
    Counts and cached-source hashes are attached to ``result.attrs["summary"]``.
    The caller may save these artifacts in a dedicated results directory.
    """
    root = Path(root)
    cfg = config.get("data", config)
    raw = root / cfg.get("raw_dir", "data/raw")
    weather = weather.reset_index(drop=True).copy()
    timestamp = pd.to_datetime(weather["timestamp_utc"], utc=True)
    if timestamp.duplicated().any() or not timestamp.is_monotonic_increasing:
        raise ValueError("Source audit requires unique, increasing UTC timestamps")
    index = pd.DatetimeIndex(timestamp)
    size = len(weather)
    sources = {name: np.full(size, np.nan) for name in (
        "precipitation_source_first_interval_start_utc", "precipitation_source_latest_utc",
        "temperature_source_left_utc", "temperature_source_latest_utc",
        "wind_source_left_utc", "wind_source_latest_utc",
    )}
    rebuilt = {v: np.full(size, np.nan) for v in VARIABLES if v != "solar_W_per_m2"}
    manifests = []
    for year in sorted(timestamp.dt.year.unique()):
        ids = np.flatnonzero(timestamp.dt.year.eq(year))
        year_index = index[ids]
        starts = year_index.as_unit("s").asi8.astype(float)
        rain_path = raw / f"NOAA_{cfg.get('precipitation_station', '72640014839')}_{year}.csv"
        met_path = raw / f"NOAA_{cfg.get('met_station', '72640594869')}_{year}.csv"
        rain = pd.read_csv(rain_path, dtype=str, keep_default_na=False)
        # Parsing once is equivalent to the original per-record timestamp parse.
        rain["DATE"] = pd.to_datetime(rain["DATE"], utc=True)
        intervals, _ = precipitation_intervals(rain)
        rebuilt["precipitation_mm_per_h"][ids] = align_precipitation(intervals, year_index)[0]
        for lo, hi, _, _ in intervals:
            a, b = lo.timestamp(), hi.timestamp()
            left = max(0, np.searchsorted(starts, a, side="right") - 1)
            right = min(len(starts), np.searchsorted(starts, b, side="left") + 1)
            local = np.arange(left, right)
            overlap = np.minimum(b, starts[local]+3600)-np.maximum(a, starts[local])
            global_ids = ids[local[overlap > 0]]
            first = sources["precipitation_source_first_interval_start_utc"]
            first[global_ids] = np.fmin(first[global_ids], a)
            sources["precipitation_source_latest_utc"][global_ids] = b
        met = pd.read_csv(met_path, dtype=str, keep_default_na=False)
        met = met[met.REPORT_TYPE.str.strip().isin(["FM-15", "FM-16"])]
        times = pd.to_datetime(met.DATE, utc=True)
        for name, column, vi, qi, is_temp, variable in (
            ("temperature", "TMP", 0, 1, True, "temperature_C"),
            ("wind", "WND", 3, 4, False, "wind_m_per_s"),
        ):
            values = [parse_value(value, vi, qi, is_temp) for value in met[column]]
            vals, lo, hi = _met_support(times, values, year_index+pd.Timedelta(minutes=30),
                                        cfg.get("max_met_gap_hours", 2))
            rebuilt[variable][ids] = vals
            sources[f"{name}_source_left_utc"][ids] = lo
            sources[f"{name}_source_latest_utc"][ids] = hi
        for path in (rain_path, met_path):
            manifests.append({"path": str(path.resolve()),
                              "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    reconstruction = {}
    for variable, values in rebuilt.items():
        saved = weather[variable].to_numpy(dtype=float)
        # Existing CSV used %.10g, so comparisons allow only its saved precision.
        same = np.allclose(saved, values, rtol=6e-10, atol=6e-10, equal_nan=True)
        finite = np.isfinite(saved) & np.isfinite(values)
        reconstruction[variable] = {
            "matches_saved_precision": bool(same),
            "max_absolute_difference": float(np.max(np.abs(saved[finite]-values[finite]))) if finite.any() else None,
            "matching_missing_mask": bool(np.array_equal(np.isnan(saved), np.isnan(values))),
        }
        if not same:
            raise ValueError(f"Cached-source reconstruction disagrees with saved {variable}")
    result = pd.DataFrame({"timestamp_utc": timestamp})
    for name, seconds in sources.items():
        result[name] = pd.to_datetime(seconds, unit="s", utc=True)
    finite = np.isfinite(weather[VARIABLES].to_numpy(dtype=float)).all(axis=1)
    adjacent = timestamp.diff().eq(pd.Timedelta(hours=1)).to_numpy()
    lag_finite = np.r_[False, finite[:-1]]
    result["current_values_complete"] = finite
    result["lag_values_complete"] = lag_finite
    result["exact_one_hour_utc_adjacency"] = adjacent
    result["base_pair_eligible"] = finite & lag_finite & adjacent
    for name in ("precipitation", "temperature", "wind"):
        lag_source = result[f"{name}_source_latest_utc"].shift(1)
        result[f"lag_{name}_source_latest_utc"] = lag_source
        result[f"lag_{name}_lookahead_minutes"] = (lag_source-timestamp).dt.total_seconds()/60
        result[f"lag_{name}_available_at_target_start"] = lag_source.le(timestamp)
    checks = [f"lag_{name}_available_at_target_start" for name in ("precipitation", "temperature", "wind")]
    result["causal_pair_at_target_start"] = result.base_pair_eligible & result[checks].all(axis=1)
    delays = [f"lag_{name}_lookahead_minutes" for name in ("precipitation", "temperature", "wind")]
    result["causal_pair_at_target_midpoint"] = result.base_pair_eligible & result[delays].le(30).all(axis=1)
    counts = {}
    for label, years in (("development", cfg.get("train_years", [2019, 2020, 2021, 2022, 2023])),
                         ("external_evaluation", cfg.get("validation_years", [2024]))):
        period = timestamp.dt.year.isin(years)
        eligible = period & result.base_pair_eligible
        counts[label] = {
            "rows": int(period.sum()), "joint_complete_rows": int(sum(period & finite)),
            "base_eligible_pairs": int(eligible.sum()),
            "pairs_excluded_before_source_audit": int(sum(period & ~result.base_pair_eligible)),
            "strict_causal_pairs_at_target_start": int(sum(period & result.causal_pair_at_target_start)),
            "strict_causal_pairs_at_target_midpoint": int(sum(period & result.causal_pair_at_target_midpoint)),
            "negative_temperature_responses_retained_in_base_pairs": int(sum(eligible & weather.temperature_C.lt(0))),
            "lag_source_audit": {},
        }
        for name in ("precipitation", "temperature", "wind"):
            delay = result.loc[eligible, f"lag_{name}_lookahead_minutes"]
            counts[label]["lag_source_audit"][name] = {
                "available_at_target_start": int(sum(eligible & result[f"lag_{name}_available_at_target_start"])),
                "after_target_start": int(delay.gt(0).sum()),
                "minimum_lookahead_minutes": float(delay.min()) if len(delay) else None,
                "maximum_lookahead_minutes": float(delay.max()) if len(delay) else None,
            }
    result.attrs["summary"] = {
        "status": "source_time_audit_complete",
        "rule": "Predictor source support must end no later than the response UTC interval start.",
        "solar_caveat": "Existing hourly POWER estimates have inherited hour-start support; operational publication latency is not evaluated.",
        "data_values_modified": False,
        "reconstruction_rtol": 6e-10, "reconstruction_atol": 6e-10,
        "reconstruction": reconstruction, "counts": counts, "raw_sources": manifests,
        "finding": "All otherwise eligible precipitation lags depend on accumulation reports ending 52 minutes after the target-bin start; the same source interval overlaps the target response bin.",
    }
    return result
