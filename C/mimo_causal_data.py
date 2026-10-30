"""Build an isolated FIT-only hourly observation-time-aligned vector from cached data.

No original input is modified. Source support ends at or before each UTC grid
stamp; event-time causality does not establish operational publication latency.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from C.uncertainty.data import VARIABLES, parse_value, precipitation_intervals


# This definition is saved before constructing any new predictor/target pairs.
ALIGNMENT_DEFINITION = {
    "evaluation_type": "Retrospective observation-time-aligned conditional evaluation",
    "prediction_cutoff": "c is an exact UTC hour. Every predictor report or interpreted source interval must end no later than c.",
    "target_reference": "t = c + 1 hour is the next selection/reference time, not a claim that all targets cover (c,t].",
    "nominal_clock_window": "(c, c+1 hour]; retained only as a reference for comparing actual measurement supports.",
    "rain_predictor": "P_completed_1h_at_cutoff: latest accepted native one-hour accumulation ending in (c-1 hour,c]; age is c minus actual accumulation end, strictly less than 60 minutes.",
    "rain_target": "Latest accepted native one-hour accumulation ending in (c,c+1 hour]. Preserve its actual start/end even when its start precedes c. Do not redistribute any accumulation.",
    "temperature_wind": "Predictors are latest QC report values at/before c with age at most 60 minutes. Targets use the same rule at t. DATE is the observation/report timestamp; instrument averaging-window bounds are not in the cache.",
    "solar": "Use the existing hourly grid value labelled c-1 hour as the predictor. Its assumed hour-start support is [c-1 hour,c]; the target is labelled c and assumed to cover [c,c+1 hour]. Cached POWER metadata confirms UTC and Wh/m2 but does not independently specify this interval anchor.",
    "historical_availability_assumption": "For retrospective selection only, assume a station report is usable at its reported observation/accumulation end, and a POWER hourly estimate at the assumed end of its hour. Publication latency is unknown and is not modeled.",
    "limits": "Station report/accumulation timing is checked directly. Solar timing remains conditional on the declared interval-anchor assumption. These mixed-support targets are not an exact one-hour-ahead forecast of four uniform clock-hour quantities.",
    "missing_values": "Leave unavailable or stale values missing; never zero-fill rainfall or carry a source indefinitely.",
}


def _digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _asof_values(target, times, values, max_age_minutes=60, strict=False):
    """Use only the most recent nonconflicting QC value at/before the stamp."""
    series = pd.Series(values, index=pd.DatetimeIndex(times), dtype=float).dropna().sort_index()
    series = series.groupby(level=0).agg(lambda a: a.iloc[0] if a.nunique() == 1 else np.nan).dropna()
    if series.empty:
        return np.full(len(target), np.nan), pd.DatetimeIndex([pd.NaT]*len(target), tz="UTC"), np.full(len(target), np.nan)
    source = pd.DataFrame({"source": series.index, "value": series.to_numpy()})
    merged = pd.merge_asof(pd.DataFrame({"target": target}), source, left_on="target", right_on="source", direction="backward", allow_exact_matches=True)
    age = (merged.target-merged.source).dt.total_seconds().to_numpy()/60
    allowed = (age >= 0) & ((age < max_age_minutes) if strict else (age <= max_age_minutes))
    values = merged.value.to_numpy()
    values[~allowed] = np.nan
    return values, pd.DatetimeIndex(merged.source), age


def build_causal_weather(root: Path, config: dict) -> tuple[pd.DataFrame, dict]:
    """Return the dedicated FIT frame and its provenance; perform no writes.

    ``root`` is the Stage C directory. Each row is indexed by a UTC hourly
    selection/reference stamp. P is an unredistributed native one-hour accumulation,
    T/W are backward-as-of reported observations, and Rs is the completed preceding
    hour's existing QC gridded estimate. Those different supports are explicit.
    """
    root = Path(root)
    cfg = config.get("data", config)
    years = sorted(set(cfg.get("train_years", [2019, 2020, 2021, 2022, 2023])) |
                   set(cfg.get("validation_years", [2024])))
    if years != list(range(min(years), max(years)+1)):
        raise ValueError("Causal FIT builder requires the configured consecutive calendar years")
    start = pd.Timestamp(f"{min(years)}-01-01", tz="UTC")
    stop = pd.Timestamp(f"{max(years)+1}-01-01", tz="UTC")
    target = pd.date_range(start, stop, freq="h", inclusive="left")
    raw = root/cfg.get("raw_dir", "data/raw")
    rain_station = str(cfg.get("precipitation_station", "72640014839"))
    met_station = str(cfg.get("met_station", "72640594869"))
    rain_frames, met_frames, sources = [], [], []
    for year in years:
        for station, frames in ((rain_station, rain_frames), (met_station, met_frames)):
            path = raw/f"NOAA_{station}_{year}.csv"
            frame = pd.read_csv(path, dtype=str, keep_default_na=False)
            if set(frame.STATION) != {station}:
                raise ValueError(f"Unexpected station in {path}")
            frame["DATE"] = pd.to_datetime(frame.DATE, utc=True)
            frames.append(frame)
            sources.append({"path": str(path.resolve()), "sha256": _digest(path)})
    # Global selection avoids inventing a break or independently overlapping
    # intervals at calendar-year boundaries.
    intervals, precipitation_qc = precipitation_intervals(pd.concat(rain_frames, ignore_index=True))
    if not intervals:
        raise ValueError("No QC-approved native one-hour precipitation intervals")
    p_start, p_end, p_depth, p_trace = map(list, zip(*intervals))
    precipitation, p_source, p_age = _asof_values(target, p_end, p_depth, strict=True)
    rain_lookup = pd.DataFrame({"start": p_start, "trace": p_trace}, index=pd.DatetimeIndex(p_end))
    selected_rain = rain_lookup.reindex(p_source).reset_index(drop=True)
    frame = pd.DataFrame({"timestamp_utc": target, "precipitation_mm_per_h": precipitation,
                          "precipitation_source_start_utc": pd.to_datetime(selected_rain.start, utc=True),
                          "precipitation_source_end_utc": p_source,
                          "precipitation_source_latest_utc": p_source,
                          "precipitation_source_age_minutes": p_age,
                          "precipitation_trace_lower_bound": selected_rain.trace.astype("boolean")})
    met = pd.concat(met_frames, ignore_index=True)
    met = met[met.REPORT_TYPE.str.strip().isin(["FM-15", "FM-16"])]
    for name, raw_column, vi, qi, is_temperature, output_column in (
        ("temperature", "TMP", 0, 1, True, "temperature_C"),
        ("wind", "WND", 3, 4, False, "wind_m_per_s"),
    ):
        values = [parse_value(value, vi, qi, is_temperature) for value in met[raw_column]]
        vals, source, age = _asof_values(target, met.DATE, values)
        frame[output_column] = vals
        for suffix in ("start", "end", "latest"):
            frame[f"{name}_source_{suffix}_utc"] = source
        frame[f"{name}_source_age_minutes"] = age
    processed = root/cfg.get("processed_path", "data/processed/weather_hourly.csv")
    old = pd.read_csv(processed, float_precision="round_trip")
    old["timestamp_utc"] = pd.to_datetime(old.timestamp_utc, utc=True)
    if old.timestamp_utc.duplicated().any():
        raise ValueError("Existing solar data has duplicate UTC stamps")
    solar_metadata=[]
    for year in years:
        path=raw/f"POWER_solar_{year}.json"
        document=json.loads(path.read_text())
        header=document['header']
        parameter=document['parameters']['ALLSKY_SFC_SW_DWN']
        if header.get('time_standard') != 'UTC' or parameter.get('units') != 'Wh/m^2':
            raise ValueError(f"Unsupported cached solar time standard or units: {path}")
        solar_metadata.append(dict(path=str(path.resolve()),sha256=_digest(path),
            header=header,parameter=parameter,
            interval_bounds_status='Cached response has hourly labels but no explicit source-interval anchor/bounds.'))
    solar_start = target-pd.Timedelta(hours=1)
    solar = old.set_index("timestamp_utc")["solar_W_per_m2"].reindex(solar_start).to_numpy()
    frame["solar_W_per_m2"] = solar
    solar_available = np.isfinite(solar)
    frame["solar_source_start_utc"] = pd.Series(solar_start).where(solar_available)
    frame["solar_source_end_utc"] = pd.Series(target).where(solar_available)
    frame["solar_source_latest_utc"] = pd.Series(target).where(solar_available)
    frame["solar_source_age_minutes"] = np.where(solar_available, 0.0, np.nan)
    latest_columns = [f"{name}_source_latest_utc" for name in ("precipitation", "temperature", "solar", "wind")]
    frame["source_latest_utc"] = frame[latest_columns].max(axis=1)
    frame["forecast_origin_utc"] = target-pd.Timedelta(hours=1)
    frame["lag_source_latest_utc"] = frame.source_latest_utc.shift(1)
    frame["lag_sources_available_at_forecast_origin"] = frame.lag_source_latest_utc.le(frame.forecast_origin_utc)
    complete = np.isfinite(frame[VARIABLES].to_numpy()).all(axis=1)
    frame["all_values_complete"] = complete
    frame["complete_adjacent_pair"] = complete & np.r_[False, complete[:-1]] & frame.timestamp_utc.diff().eq(pd.Timedelta(hours=1))
    if not frame.loc[complete, "source_latest_utc"].le(frame.loc[complete, "timestamp_utc"]).all():
        raise AssertionError("Future source support entered a causal FIT response")
    if not frame.loc[frame.complete_adjacent_pair, "lag_sources_available_at_forecast_origin"].all():
        raise AssertionError("Future source support entered a causal FIT predictor")
    period_counts = {}
    for label, selected_years in (("development", cfg.get("train_years", [2019, 2020, 2021, 2022, 2023])),
                                  ("external_evaluation", cfg.get("validation_years", [2024]))):
        chosen = frame.timestamp_utc.dt.year.isin(selected_years)
        pairs = chosen & frame.complete_adjacent_pair
        period_counts[label] = {
            "rows": int(chosen.sum()), "joint_complete_rows": int(sum(chosen & complete)),
            "complete_adjacent_pairs": int(pairs.sum()),
            "excluded_pairs": int(sum(chosen & ~frame.complete_adjacent_pair)),
            "negative_temperature_responses_in_pairs": int(sum(pairs & frame.temperature_C.lt(0))),
            "missing_by_variable": {column: int(frame.loc[chosen, column].isna().sum()) for column in VARIABLES},
        }
    manifest = {
        "status": "causal_fit_only_weather_prepared",
        "original_inputs_modified": False,
        "timestamp_definition": "UTC hourly observation-selection reference; not a measured availability or publication time",
        "alignment_definition": ALIGNMENT_DEFINITION,
        "units": {"precipitation_mm_per_h": "mm/h (native one-hour depth divided by one hour)",
                  "temperature_C": "deg C", "solar_W_per_m2": "W/m^2", "wind_m_per_s": "m/s"},
        "precipitation_rule": "Latest QC-approved nonoverlapping native one-hour accumulation whose report end is in (t-1h,t]; no UTC-bin overlap redistribution. Trace remains a flagged zero quantified lower bound.",
        "met_rule": "Last nonconflicting QC-approved FM-15/FM-16 reported observation at or before t, maximum age 60 minutes; no interpolation, forward fill beyond that limit, or centered averaging.",
        "age_rule_status": "Fixed sampling eligibility rule chosen before fitting; not a tuned physical model parameter.",
        "solar_rule": "Existing QC POWER value at original hour-start t-1h is assigned to its completed-hour end t. The original hourly-support anchor assumption is retained explicitly.",
        "causality_scope": "All station predictor reports/native accumulation ends are checked against cutoff c=t-1h. Solar support is conditional on the declared hour-start assumption; full upstream observation support and operational availability are not certified.",
        "support_caveat": "Precipitation typically ends eight minutes before t; at cutoff c=t-1h its target spans (c-8 min,c+52 min], including eight pre-cutoff minutes. T/W are recent reports; solar support is assumed. These are not four uniformly aligned hourly averages.",
        "split_rule": "Split by response timestamp UTC year; 2019-2023 development, 2024 external. Cached years are concatenated before backward-only source selection.",
        "first_boundary": "No 2018 solar observation is fabricated for 2019-01-01 00:00 UTC; this row remains incomplete.",
        "source_table": "Per-variable start/end/latest source timestamps and ages are columns of the returned timestamped weather frame. T/W start=end denotes the reported observation timestamp, not a verified zero-duration instrument average. Solar bounds are interpreted under the documented assumption.",
        "precipitation_qc": precipitation_qc,
        "cached_solar_metadata": solar_metadata,
        "counts": period_counts,
        "raw_sources": sources,
        "existing_processed_solar_source": {"path": str(processed.resolve()), "sha256": _digest(processed)},
    }
    # Keep the canonical timestamp and weather columns first for inspection.
    first = ["timestamp_utc"]+VARIABLES
    return frame[first+[name for name in frame if name not in first]], manifest
