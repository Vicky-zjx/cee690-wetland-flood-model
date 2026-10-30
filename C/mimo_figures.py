"""Descriptive figures for external, observation-time-aligned official mimo_srs estimates.

This module never fits, rescales, clips, fills, or recursively generates weather.
Residuals are observed minus predicted. Lagged residual pairs must be separated
by the exact requested number of UTC hours; missing timestamps remain missing.
"""
from __future__ import annotations

from datetime import timezone
from pathlib import Path

import numpy as np
import pandas as pd


VARIABLES = (
    ("P", "Precipitation", "mm/h"),
    ("T", "Temperature", "°C"),
    ("Rs", "Solar radiation", "W/m²"),
    ("W", "Wind speed", "m/s"),
)


def _example_window(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Choose the first complete seven-day window, without inspecting fit quality."""
    breaks = frame["timestamp_utc"].diff().ne(pd.Timedelta(hours=1))
    runs = [part for _, part in frame.groupby(breaks.cumsum(), sort=False)]
    qualifying = next((part for part in runs if len(part) >= 168), None)
    fallback = qualifying is None
    if fallback:
        # max() retains the earliest run in a tie because runs are chronological.
        qualifying = max(runs, key=len)
        selected = qualifying
        rule = (
            "No 168-hour eligible contiguous window exists. Select the longest "
            "contiguous eligible UTC run, with the earliest start breaking ties. "
            "Selection uses timestamp coverage only, never fit quality."
        )
    else:
        selected = qualifying.iloc[:168]
        rule = (
            "Select the earliest 168 consecutive eligible UTC hourly timestamps. "
            "Selection uses timestamp coverage only, never fit quality."
        )
    return selected, {
        "rule": rule,
        "start_utc": selected["timestamp_utc"].iloc[0].isoformat(),
        "end_utc": selected["timestamp_utc"].iloc[-1].isoformat(),
        "n_hours": int(len(selected)),
        "fallback": fallback,
        "contiguous_run_length_hours": int(len(qualifying)),
    }


def _residual_acf(frame: pd.DataFrame, max_lag: int = 48) -> pd.DataFrame:
    """Pairwise Pearson residual autocorrelation at exact UTC hour separations."""
    indexed = frame.set_index("timestamp_utc")
    records = []
    for key, _, _ in VARIABLES:
        values = indexed[f"{key}_residual"]
        current = values.to_numpy(dtype=float)
        for lag in range(1, max_lag + 1):
            previous = values.reindex(
                indexed.index - pd.Timedelta(hours=lag)
            ).to_numpy(dtype=float)
            valid = np.isfinite(current) & np.isfinite(previous)
            left, right = current[valid], previous[valid]
            correlation = np.nan
            if len(left) >= 2 and np.std(left) > 0 and np.std(right) > 0:
                correlation = float(np.corrcoef(left, right)[0, 1])
            records.append({
                "variable": key,
                "lag_hours": lag,
                "pair_count": int(valid.sum()),
                "autocorrelation": correlation,
            })
    return pd.DataFrame.from_records(records)


def make_figures(
    predictions: pd.DataFrame, plot_dir: Path, table_dir: Path
) -> dict:
    """Save three four-variable figure groups and an exact-hour ACF table.

    Required columns are timestamp_utc, predictor_timestamp_utc, and observed,
    predicted, persistence, and residual columns prefixed P, T, Rs, and W.
    All supplied rows must already be eligible external validation predictions.
    The caller, rather than this plotting helper, determines model eligibility.
    """
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    required = ["timestamp_utc", "predictor_timestamp_utc"] + [
        f"{key}_{suffix}"
        for key, _, _ in VARIABLES
        for suffix in ("observed", "predicted", "persistence", "residual")
    ]
    missing = sorted(set(required) - set(predictions.columns))
    if missing:
        raise ValueError(f"Missing plotting columns: {missing}")
    frame = predictions[required].copy()
    frame["timestamp_utc"] = pd.to_datetime(frame["timestamp_utc"], utc=True)
    frame["predictor_timestamp_utc"] = pd.to_datetime(
        frame["predictor_timestamp_utc"], utc=True
    )
    frame = frame.sort_values("timestamp_utc").reset_index(drop=True)
    if frame.empty or frame["timestamp_utc"].duplicated().any():
        raise ValueError("Figure inputs require nonempty, unique UTC timestamps")
    if frame[["timestamp_utc", "predictor_timestamp_utc"]].isna().any().any():
        raise ValueError("Figure inputs contain missing timestamps")
    if not (
        frame["timestamp_utc"] - frame["predictor_timestamp_utc"]
    ).eq(pd.Timedelta(hours=1)).all():
        raise ValueError("Every plotted prediction must use the preceding UTC hour")
    if not np.isfinite(frame[required[2:]].to_numpy(dtype=float)).all():
        raise ValueError("Figure inputs contain nonfinite values; no rows were dropped")

    plot_dir, table_dir = Path(plot_dir), Path(table_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)
    window, window_info = _example_window(frame)
    acf = _residual_acf(frame)
    acf_path = table_dir / "residual_acf.csv"
    acf.to_csv(acf_path, index=False, float_format="%.17g")

    paths = {
        "observed_vs_predicted": plot_dir / "observed_vs_predicted.png",
        "continuous_example": plot_dir / "continuous_example.png",
        "residual_diagnostics": plot_dir / "residual_diagnostics.png",
    }
    colors = {"observed": "#20262D", "predicted": "#007A9F", "persistence": "#B56116"}
    style = {
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.18,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
    with plt.rc_context(style):
        fig, axes = plt.subplots(2, 2, figsize=(10.8, 9.3), layout="constrained")
        for ax, (key, label, unit) in zip(axes.flat, VARIABLES):
            observed = frame[f"{key}_observed"].to_numpy()
            predicted = frame[f"{key}_predicted"].to_numpy()
            low = min(observed.min(), predicted.min())
            high = max(observed.max(), predicted.max())
            padding = max((high - low) * 0.045, 0.01)
            limits = (low - padding, high + padding)
            ax.scatter(observed, predicted, s=9, alpha=0.28, color=colors["predicted"],
                       edgecolors="none", rasterized=True)
            ax.plot(limits, limits, "--", linewidth=1.1, color="#20262D", label="1:1")
            ax.set(xlim=limits, ylim=limits, xlabel=f"Observed ({unit})",
                   ylabel=f"Predicted ({unit})", title=f"{label} · N = {len(frame):,}")
            ax.set_aspect("equal", adjustable="box")
            ax.legend(loc="upper left", frameon=False, fontsize=9)
        fig.suptitle("Official mimo_srs · Retrospective 2024 aligned evaluation", fontsize=14)
        fig.savefig(paths["observed_vs_predicted"], dpi=180)
        plt.close(fig)

        fig, axes = plt.subplots(4, 1, figsize=(12.8, 10.8), sharex=True, layout="constrained")
        for ax, (key, label, unit) in zip(axes, VARIABLES):
            for suffix, title, width, linestyle, alpha in (
                ("observed", "Observed", 1.7, "-", 0.9),
                ("predicted", "mimo_srs", 1.5, "-", 0.9),
                ("persistence", "Persistence", 1.1, "--", 0.85),
            ):
                ax.plot(window["timestamp_utc"], window[f"{key}_{suffix}"],
                        color=colors[suffix], label=title, linewidth=width,
                        linestyle=linestyle, alpha=alpha)
            ax.set_ylabel(f"{label}\n({unit})")
        axes[0].legend(ncol=3, loc="upper right", frameon=False)
        locator = mdates.AutoDateLocator(minticks=5, maxticks=9, tz=timezone.utc)
        axes[-1].xaxis.set_major_locator(locator)
        axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%b %d\n%H:%M", tz=timezone.utc))
        axes[-1].set_xlabel("Target selection/reference time (UTC)")
        window_label = "Longest eligible continuous window" if window_info["fallback"] else "Earliest eligible 168-hour window"
        fig.suptitle(f"External 2024 continuous example · {window_label}\n"
                     "Predictors end by cutoff; actual target supports differ by variable", fontsize=13)
        fig.savefig(paths["continuous_example"], dpi=180)
        plt.close(fig)

        fig, axes = plt.subplots(4, 2, figsize=(12.6, 13.0), layout="constrained")
        for row, (key, label, unit) in enumerate(VARIABLES):
            hist_ax, acf_ax = axes[row]
            hist_ax.hist(frame[f"{key}_residual"], bins=50, color=colors["predicted"],
                         alpha=0.82, edgecolor="white", linewidth=0.25)
            hist_ax.axvline(0, color="#20262D", linestyle="--", linewidth=1)
            hist_ax.set(title=f"{label} residuals", xlabel=f"Observed − predicted ({unit})",
                        ylabel="Number of hours")
            component = acf.loc[acf["variable"].eq(key)]
            acf_ax.axhline(0, color="#20262D", linewidth=0.8)
            acf_ax.plot(component["lag_hours"], component["autocorrelation"],
                        color=colors["predicted"], linewidth=1.4, marker="o", markersize=2.7)
            acf_ax.set(title=f"{label} residual autocorrelation", xlabel="Exact UTC lag (hours)",
                       ylabel="Pairwise correlation", xlim=(0, 49), ylim=(-1.05, 1.05))
            acf_ax.set_xticks([1, 12, 24, 36, 48])
        fig.suptitle("External 2024 residual diagnostics\n"
                     "Descriptive diagnostics; missing hours are not compressed or filled", fontsize=14)
        fig.savefig(paths["residual_diagnostics"], dpi=180)
        plt.close(fig)

    per_variable = {}
    for key, _, _ in VARIABLES:
        component = acf.loc[acf["variable"].eq(key)].set_index("lag_hours")
        per_variable[key] = {
            "lag1_pair_count": int(component.loc[1, "pair_count"]),
            "lag48_pair_count": int(component.loc[48, "pair_count"]),
            "undefined_lags": [int(lag) for lag in component.index[component["autocorrelation"].isna()]],
        }
    return {
        "figures": {key: str(path.resolve()) for key, path in paths.items()},
        "example_window": window_info,
        "residual_acf": {
            "definition": "Pearson correlation of observed-minus-predicted residual pairs separated by exactly k UTC hours, using only available pairs for each lag. No gap compression or filling.",
            "max_lag_hours": 48,
            "table_path": str(acf_path.resolve()),
            "by_variable": per_variable,
            "confidence_bands": "Not shown; no residual independence or whiteness test is inferred from these descriptive correlations.",
        },
    }
