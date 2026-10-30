"""Held-out physical-unit validation; failures remain failures in saved results."""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd

from .weather_generator import COLUMNS, solar_geometry, utc_index, _settings


def _finite(values):
    x = np.asarray(values, dtype=float)
    return x[np.isfinite(x)]


def _stats(values):
    x = _finite(values)
    return {"count": int(len(x)), "mean": float(x.mean()) if len(x) else None,
            "std": float(x.std()) if len(x) else None,
            "quantiles": {str(q): float(np.quantile(x, q)) if len(x) else None for q in [.01, .05, .5, .9, .95, .99]}}


def precipitation_events(frame, wet_threshold=0.0, separation_dry_hours=6):
    """Hourly spells and storms; missing/nonadjacent hours always split records.

    A storm joins wet hours separated by fewer than six dry hours; storm
    duration spans first through last wet hour. No gaps are inferred dry.
    """
    index = utc_index(frame["timestamp_utc"])
    rain = frame[COLUMNS[0]].to_numpy(dtype=float)
    wet_spells, dry_spells, storms = [], [], []
    spell_state, spell_len = None, 0
    storm = None
    dry_gap = 0

    def close_spell():
        if spell_state is not None:
            (wet_spells if spell_state else dry_spells).append(spell_len)

    def close_storm():
        if storm is not None:
            storms.append({"total_mm": storm["total"], "duration_hours": storm["last"] - storm["first"] + 1,
                           "peak_mm_per_h": storm["peak"]})

    for i, amount in enumerate(rain):
        adjacent = i == 0 or index[i].value - index[i - 1].value == pd.Timedelta(hours=1).value
        if not adjacent or not np.isfinite(amount):
            close_spell()
            close_storm()
            spell_state, spell_len, storm, dry_gap = None, 0, None, 0
        if not np.isfinite(amount):
            continue
        wet = amount > wet_threshold
        if spell_state is not None and wet != spell_state:
            close_spell()
            spell_len = 0
        spell_state, spell_len = wet, spell_len + 1
        if wet:
            if storm is None:
                storm = {"first": i, "last": i, "total": float(amount), "peak": float(amount)}
            else:
                storm["last"] = i
                storm["total"] += float(amount)
                storm["peak"] = max(storm["peak"], float(amount))
            dry_gap = 0
        elif storm is not None:
            dry_gap += 1
            if dry_gap >= separation_dry_hours:
                close_storm()
                storm, dry_gap = None, 0
    close_spell()
    close_storm()
    return {"wet_spells_hours": wet_spells, "dry_spells_hours": dry_spells,
            "storms": storms, "wet_threshold_mm_per_h": wet_threshold,
            "storm_separation_dry_hours": separation_dry_hours,
            "censoring": "Segments at year edges/missing gaps are retained but can be censored; no durations bridge missing data"}


def _correlation(frame, method="pearson", anomalies=False, timezone="America/Chicago"):
    data = frame[COLUMNS].copy()
    if anomalies:
        index = utc_index(frame["timestamp_utc"]).tz_convert(timezone)
        # Descriptive diagnostics only: removing each sample's month/hour mean
        # does not refit or influence any generator parameter.
        groups = pd.DataFrame({"month": index.month, "hour": index.hour}, index=data.index)
        means = data.groupby([groups.month, groups.hour]).transform("mean")
        data = data - means
    return data.corr(method=method, min_periods=10).to_numpy()


def _lagcorr(frame, lag):
    index = utc_index(frame["timestamp_utc"])
    x = frame[COLUMNS].to_numpy(dtype=float)
    result = np.full((4, 4), np.nan)
    if len(index) <= lag:
        return result
    adjacent = index[lag:].asi8 - index[:-lag].asi8 == pd.Timedelta(hours=lag).value
    for j in range(4):
        for k in range(4):
            good = adjacent & np.isfinite(x[lag:, j]) & np.isfinite(x[:-lag, k])
            if good.sum() > 10 and np.std(x[lag:, j][good]) > 0 and np.std(x[:-lag, k][good]) > 0:
                result[j, k] = np.corrcoef(x[lag:, j][good], x[:-lag, k][good])[0, 1]
    return result


def _pattern(frame, part, timezone):
    local = (utc_index(frame["timestamp_utc"]) + pd.Timedelta(minutes=30)).tz_convert(timezone)
    keys = local.month if part == "month" else local.hour
    data = frame[COLUMNS].copy()
    data[part] = keys
    return data.groupby(part)[COLUMNS].mean().reindex(range(1, 13) if part == "month" else range(24))


def _clean(value):
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_clean(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    return value


def validate_weather(observed, generated, config, outdir=None):
    """Compare one independent generated calendar against held-out observations.

    Diagnostic criteria are taken from the pre-run configuration. A single
    held-out realization is an adequacy screen, not a probabilistic guarantee.
    """
    settings = _settings(config)
    criteria = {"wet_frequency_absolute": .03, "positive_p99_relative": .5,
                "storm_p95_relative": .5, "wet_spell_p95_relative": .5,
                "lag1_absolute": .25, "physical_violations": 0, **config.get("validation", {})}
    observed_index, generated_index = utc_index(observed.timestamp_utc), utc_index(generated.timestamp_utc)
    if not observed_index.equals(generated_index):
        raise ValueError("Validation observation and generation calendars must match exactly")
    generated_full = generated.copy()
    generated = generated.copy()
    matched_support = {}
    for variable in COLUMNS:
        absent = ~np.isfinite(observed[variable].to_numpy(dtype=float))
        generated.loc[absent, variable] = np.nan
        matched_support[variable] = {"observed_valid_hours": int((~absent).sum()), "generated_hours_masked_for_comparison": int(absent.sum()), "generated_full_valid_hours": int(np.isfinite(generated_full[variable]).sum())}
    metrics, correlations, lags, seasonal, diurnal, event_arrays = {}, {}, {}, {}, {}, {}
    for label, frame in [("observed", observed), ("generated", generated)]:
        rain = _finite(frame[COLUMNS[0]])
        events = precipitation_events(frame)
        event_arrays[label] = events
        metrics[label] = {variable: _stats(frame[variable]) for variable in COLUMNS}
        metrics[label].update({"wet_frequency": float(np.mean(rain > 0)) if len(rain) else None,
                               "zero_precipitation_frequency": float(np.mean(rain == 0)) if len(rain) else None,
                               "zero_frequency_meaning": "Quantified zero includes trace lower bounds; not strictly observed meteorological dryness",
                               "trace_lower_bound_hours": int(frame["precipitation_trace_lower_bound"].fillna(False).astype(bool).sum()) if "precipitation_trace_lower_bound" in frame else None,
                               "positive_precipitation": _stats(rain[rain > 0]),
                               "wet_spell_hours": _stats(events["wet_spells_hours"]),
                               "dry_spell_hours": _stats(events["dry_spells_hours"]),
                               "storm_total_mm": _stats([s["total_mm"] for s in events["storms"]]),
                               "storm_duration_hours": _stats([s["duration_hours"] for s in events["storms"]]),
                               "storm_peak_mm_per_h": _stats([s["peak_mm_per_h"] for s in events["storms"]])})
        correlations[label] = {f"{method}_{'anomaly' if anomaly else 'raw'}": _correlation(frame, method, anomaly, settings["timezone"])
                               for method in ["pearson", "spearman"] for anomaly in [False, True]}
        lags[label] = {str(lag): _lagcorr(frame, lag) for lag in [1, 3, 6, 12, 24, 48, 72]}
        seasonal[label] = _pattern(frame, "month", settings["timezone"])
        diurnal[label] = _pattern(frame, "hour", settings["timezone"])
    checks = []

    def compare(name, a, b, tolerance, relative=False):
        available = a is not None and b is not None and np.isfinite(a) and np.isfinite(b)
        error = abs(b - a) / max(abs(a), 1e-12) if available and relative else abs(b - a) if available else None
        checks.append({"name": name, "observed": a, "generated": b, "error": error,
                       "error_kind": "relative" if relative else "absolute", "tolerance": tolerance,
                       "passed": bool(available and error <= tolerance),
                       "reason": None if available else "Insufficient validation data"})

    obs, gen = metrics["observed"], metrics["generated"]
    compare("wet_hour_frequency", obs["wet_frequency"], gen["wet_frequency"], criteria["wet_frequency_absolute"])
    compare("positive_hourly_precipitation_p99", obs["positive_precipitation"]["quantiles"]["0.99"], gen["positive_precipitation"]["quantiles"]["0.99"], criteria["positive_p99_relative"], True)
    for key in ["storm_total_mm", "storm_duration_hours", "storm_peak_mm_per_h"]:
        compare(key + "_p95", obs[key]["quantiles"]["0.95"], gen[key]["quantiles"]["0.95"], criteria["storm_p95_relative"], True)
    for key in ["wet_spell_hours", "dry_spell_hours"]:
        compare(key + "_p95", obs[key]["quantiles"]["0.95"], gen[key]["quantiles"]["0.95"], criteria["wet_spell_p95_relative"], True)
    for j, variable in enumerate(COLUMNS):
        compare(variable + "_lag1", lags["observed"]["1"][j, j], lags["generated"]["1"][j, j], criteria["lag1_absolute"])
    _, daylight = solar_geometry(generated_index, settings["latitude"], settings["longitude"])
    array = generated_full[COLUMNS].to_numpy(dtype=float)
    violations = {"nonfinite_values": int((~np.isfinite(array)).sum()),
                  "negative_precipitation": int(np.sum(array[:, 0] < 0)),
                  "negative_solar": int(np.sum(array[:, 2] < 0)),
                  "negative_wind": int(np.sum(array[:, 3] < 0)),
                  "nonzero_fully_night_radiation": int(np.sum(np.abs(array[~daylight, 2]) > 1e-10))}
    checks.append({"name": "physical_constraints", "observed": 0, "generated": sum(violations.values()),
                   "error": sum(violations.values()), "tolerance": criteria["physical_violations"],
                   "passed": sum(violations.values()) <= criteria["physical_violations"], "details": violations})
    summary = {"passed": all(item["passed"] for item in checks), "checks": checks,
               "failed_checks": [item["name"] for item in checks if not item["passed"]],
               "criteria": criteria, "metrics": metrics, "correlations": correlations,
               "matched_support": matched_support,
               "matched_support_policy": "Generated diagnostics use the identical per-variable missingness mask as heldout observations. Missing rainfall terminates both event series; full generated array remains unmasked for physical checks.",
               "generated_full_sample_statistics": {variable: _stats(generated_full[variable]) for variable in COLUMNS},
               "lagged_correlations": lags,
               "seasonal": {key: value.to_dict(orient="index") for key, value in seasonal.items()},
               "diurnal": {key: value.to_dict(orient="index") for key, value in diurnal.items()},
               "correlation_order": COLUMNS, "lag_direction": "rows=current variable; columns=variable lag hours earlier",
               "anomaly_correlations": "Each sample's local month/hour means removed for descriptive validation only; generator never refitted on heldout data",
               "storm_definition": "Quantified rainfall >0; storms separated by >=6 quantified-zero hours, including trace lower bounds. Unquantified trace is omitted from storm amount/occurrence; missing hours terminate events.",
               "validation_population": {"start_utc": observed_index.min().isoformat(), "end_utc": observed_index.max().isoformat(),
                                         "hour_count": len(observed), "type": "one full-calendar independent generated realization versus chronologically held-out data"},
               "limitations": ["Failed storm or persistence checks prohibit interpreting this generator as a reliable flood-probability estimator.",
                               "One held-out calendar contains interannual sampling variation and incomplete observations; screen failures are retained, not tuned away.",
                               "Gaussian copula has no asymptotic upper-tail dependence; marginal tail extrapolation does not establish joint-extreme fidelity.",
                               "Empirical marginals/bootstrap cannot extrapolate beyond training component extremes; storm reassembly can nevertheless produce new event totals.",
                               "Monthly/hourly patterns, raw and anomaly correlations are diagnostics; no post hoc pass threshold is invented for them."]}
    summary = _clean(summary)
    if outdir is not None:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        _plots(observed, generated, metrics, lags, seasonal, diurnal, event_arrays, outdir)
        (outdir / "validation_summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
        pd.DataFrame(checks).drop(columns=["details"], errors="ignore").to_csv(outdir / "validation_checks.csv", index=False)
    return summary


def _plots(observed, generated, metrics, lags, seasonal, diurnal, events, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = {"observed": "#0b2545", "generated": "#0e7490"}
    labels = ["Precipitation (mm/h)", "Air temperature (°C)", "Solar radiation (W/m²)", "Wind speed (m/s)"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    for j, ax in enumerate(axes.flat):
        for key, frame in [("observed", observed), ("generated", generated)]:
            x = np.sort(_finite(frame[COLUMNS[j]]))
            if j == 0:
                x = x[x > 0]
            ax.plot(x, np.arange(1, len(x) + 1) / max(len(x), 1), label=key.title(), color=colors[key])
        ax.set(xlabel=labels[j] + ("; positive hours" if j == 0 else ""), ylabel="Empirical cumulative probability")
        ax.grid(alpha=.2)
    axes.flat[0].legend()
    fig.suptitle("Held-out physical marginal distributions")
    fig.savefig(outdir / "marginals.png", dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    for j, ax in enumerate(axes.flat):
        for key in colors:
            lag = list(map(int, lags[key]))
            ax.plot(lag, [lags[key][str(k)][j, j] for k in lag], "o-", label=key.title(), color=colors[key])
        ax.set(xlabel="Lag (hours)", ylabel="Correlation", title=labels[j]); ax.grid(alpha=.2)
    axes.flat[0].legend(); fig.suptitle("Physical-variable autocorrelation")
    fig.savefig(outdir / "autocorrelation.png", dpi=160); plt.close(fig)
    for title, patterns, filename, xlabel in [("Local seasonal means", seasonal, "seasonal.png", "Month (America/Chicago)"),
                                             ("Local diurnal means", diurnal, "diurnal.png", "Hour (America/Chicago)")]:
        fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
        for j, ax in enumerate(axes.flat):
            for key in colors:
                ax.plot(patterns[key].index, patterns[key][COLUMNS[j]], "o-", label=key.title(), color=colors[key])
            ax.set(xlabel=xlabel, ylabel=labels[j]); ax.grid(alpha=.2)
        axes.flat[0].legend(); fig.suptitle(title)
        fig.savefig(outdir / filename, dpi=160); plt.close(fig)
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
    for ax, key, title in zip(axes.flat, ["wet_spells_hours", "dry_spells_hours", "total_mm", "duration_hours", "peak_mm_per_h"],
                             ["Wet spell (h)", "Dry spell (h)", "Storm total (mm)", "Storm duration (h)", "Storm peak (mm/h)"]):
        for label in colors:
            x = np.sort(events[label][key] if key.endswith("spells_hours") else [s[key] for s in events[label]["storms"]])
            ax.plot(x, np.arange(1, len(x) + 1) / max(len(x), 1), color=colors[label], label=label.title())
        ax.set(xlabel=title, ylabel="Empirical cumulative probability"); ax.grid(alpha=.2)
    axes.flat[0].legend(); axes.flat[-1].axis("off")
    fig.suptitle("Storm and spell validation; missing data split events")
    fig.savefig(outdir / "storms_and_spells.png", dpi=160); plt.close(fig)
    # Paired line comparisons avoid hiding signed or lagged correlations in a heatmap.
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    names = ["P–T", "P–solar", "P–wind", "T–solar", "T–wind", "solar–wind"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)
    for label, frame in [("observed", observed), ("generated", generated)]:
        matrix = _correlation(frame)
        axes[0].plot(names, [matrix[a, b] for a, b in pairs], "o-", label=label.title(), color=colors[label])
        lag = list(map(int, lags[label]))
        axes[1].plot(lag, [lags[label][str(k)][0, 2] for k in lag], "o-", label=label.title(), color=colors[label])
    axes[0].set(title="Contemporaneous physical Pearson correlation", ylabel="Correlation")
    axes[1].set(title="Current rainfall vs earlier solar radiation", xlabel="Lag (hours)", ylabel="Correlation")
    for ax in axes: ax.grid(alpha=.2)
    axes[0].legend(); fig.savefig(outdir / "cross_correlations.png", dpi=160); plt.close(fig)
