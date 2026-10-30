"""Matched figures for one versus three completed rainfall-hour inputs.

The caller fits the models and decides eligibility. This module never fits,
clips, fills observations, or modifies prediction values. All figures use the
same eligible external 2024 timestamps for all three methods.
"""
from __future__ import annotations

from datetime import timezone
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from .mimo_figures import VARIABLES, _example_window
except ImportError:  # Direct execution with C/ on the Python import path.
    from mimo_figures import VARIABLES, _example_window


MODELS = ("one_hour", "three_hour", "persistence")
LABELS = {
    "one_hour": "One rainfall hour",
    "three_hour": "Three rainfall hours",
    "persistence": "Persistence",
}
COLORS = {
    "observed": "#20262D",
    "one_hour": "#007A9F",
    "three_hour": "#C25B27",
    "persistence": "#777777",
}


def _matched_external(predictions: pd.DataFrame) -> dict[str, pd.DataFrame]:
    required = ["timestamp_utc", "prediction_cutoff_utc", "period", "model"] + [
        f"{key}_{kind}"
        for key, _, _ in VARIABLES
        for kind in ("observed", "predicted")
    ]
    missing = sorted(set(required) - set(predictions.columns))
    if missing:
        raise ValueError(f"Missing comparison figure columns: {missing}")
    frame = predictions.loc[predictions["period"].eq("external_2024"), required].copy()
    if frame.empty or set(frame["model"].unique()) != set(MODELS):
        raise ValueError("External figures require all three specified methods")
    for field in ("timestamp_utc", "prediction_cutoff_utc"):
        frame[field] = pd.to_datetime(frame[field], utc=True)
    if frame[["timestamp_utc", "prediction_cutoff_utc"]].isna().any().any():
        raise ValueError("Comparison figure timestamps must be nonmissing")
    if not frame["timestamp_utc"].dt.year.eq(2024).all():
        raise ValueError("external_2024 figure inputs contain another target year")
    if frame.duplicated(["model", "timestamp_utc"]).any():
        raise ValueError("Each method requires unique UTC target timestamps")
    if not np.isfinite(frame[required[4:]].to_numpy(dtype=float)).all():
        raise ValueError("Nonfinite comparison values; no samples were dropped")
    if not (frame["timestamp_utc"] - frame["prediction_cutoff_utc"]).eq(
        pd.Timedelta(hours=1)
    ).all():
        raise ValueError("Target reference timestamp must be cutoff plus one hour")
    models = {
        model: frame.loc[frame["model"].eq(model)].sort_values("timestamp_utc")
        .set_index("timestamp_utc")
        for model in MODELS
    }
    reference = models["one_hour"]
    observed_columns = [f"{key}_observed" for key, _, _ in VARIABLES]
    for model in MODELS[1:]:
        candidate = models[model]
        if not candidate.index.equals(reference.index):
            raise ValueError(f"{model}: unmatched target timestamps")
        if not candidate["prediction_cutoff_utc"].equals(reference["prediction_cutoff_utc"]):
            raise ValueError(f"{model}: unmatched prediction cutoffs")
        if not np.array_equal(candidate[observed_columns].to_numpy(),
                              reference[observed_columns].to_numpy()):
            raise ValueError(f"{model}: observations differ across methods")
    return models


def _rainfall_acf(frame: pd.DataFrame) -> list[dict]:
    """Pearson residual correlations at exact UTC lags, without gap compression."""
    residual = frame["P_observed"] - frame["P_predicted"]
    current = residual.to_numpy(dtype=float)
    records = []
    for lag in range(1, 49):
        previous = residual.reindex(
            residual.index - pd.Timedelta(hours=lag)
        ).to_numpy(dtype=float)
        valid = np.isfinite(current) & np.isfinite(previous)
        left, right = current[valid], previous[valid]
        correlation = None
        if len(left) >= 2 and np.std(left) > 0 and np.std(right) > 0:
            correlation = float(np.corrcoef(left, right)[0, 1])
        records.append({"lag_hours": lag, "pair_count": int(valid.sum()),
                        "autocorrelation": correlation})
    return records


def make_figures(predictions: pd.DataFrame, outdir: Path) -> dict:
    """Save three comparison PNGs and return JSON-serializable figure metadata.

    Input is a long DataFrame: period, model, timestamp_utc,
    prediction_cutoff_utc, and P/T/Rs/W observed/predicted columns. Accepted
    external model labels are one_hour, three_hour, and persistence. The
    one_hour method is the controlled refit, not an overwritten historical fit.
    """
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    models = _matched_external(predictions)
    reference = models["one_hour"]
    window, window_info = _example_window(reference.reset_index())
    window_index = pd.DatetimeIndex(window["timestamp_utc"])
    peak_time = reference["P_observed"].idxmax()  # Earliest timestamp breaks ties.
    peak_index = pd.date_range(peak_time - pd.Timedelta(hours=24),
                               peak_time + pd.Timedelta(hours=24), freq="h")
    acfs = {model: _rainfall_acf(models[model]) for model in MODELS[:2]}
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths = {name: outdir / f"{name}.png" for name in (
        "observed_vs_predicted", "continuous_example", "rainfall_diagnostics"
    )}
    style = {"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
             "axes.spines.top": False, "axes.spines.right": False,
             "axes.grid": True, "grid.alpha": 0.18,
             "figure.facecolor": "white", "savefig.facecolor": "white"}

    def time_axis(ax):
        ax.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=5, maxticks=9,
                                                        tz=timezone.utc))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d\n%H:%M", tz=timezone.utc))
        ax.set_xlabel("Target selection/reference time (UTC)")

    def plot_series(ax, key, index):
        # Reindexing adds NaNs only for display gaps; it never fills a value.
        ax.plot(index, reference[f"{key}_observed"].reindex(index),
                color=COLORS["observed"], label="Observed", linewidth=1.6, zorder=4)
        for model in MODELS:
            ax.plot(index, models[model][f"{key}_predicted"].reindex(index),
                    color=COLORS[model], label=LABELS[model],
                    linewidth=1.3 if model != "persistence" else 1.0,
                    linestyle="--" if model == "persistence" else "-", alpha=0.9)

    with plt.rc_context(style):
        fig, axes = plt.subplots(2, 4, figsize=(15.2, 8.2), sharex="col", sharey="col",
                                 layout="constrained")
        for col, (key, title, unit) in enumerate(VARIABLES):
            observed = reference[f"{key}_observed"].to_numpy()
            all_values = np.concatenate([observed] + [
                models[model][f"{key}_predicted"].to_numpy() for model in MODELS[:2]
            ])
            low, high = float(all_values.min()), float(all_values.max())
            pad = max((high - low) * 0.05, 0.01)
            limits = (low - pad, high + pad)
            denominator = np.sum((observed - observed.mean()) ** 2)
            for row, model in enumerate(MODELS[:2]):
                ax = axes[row, col]
                predicted = models[model][f"{key}_predicted"].to_numpy()
                r2 = 1 - np.sum((observed - predicted) ** 2) / denominator if denominator > 0 else np.nan
                ax.scatter(observed, predicted, color=COLORS[model], s=7,
                           alpha=0.25, edgecolors="none", rasterized=True)
                ax.plot(limits, limits, "--", color=COLORS["observed"], linewidth=1)
                ax.set(xlim=limits, ylim=limits,
                       ylabel=f"{LABELS[model]}\nPredicted ({unit})")
                ax.set_aspect("equal", adjustable="box")
                ax.text(0.04, 0.96, f"$R^2$ = {r2:.3f}", transform=ax.transAxes,
                        ha="left", va="top", bbox={"facecolor": "white", "alpha": 0.8,
                                                   "edgecolor": "none"})
                if row == 0:
                    ax.set_title(title)
                else:
                    ax.set_xlabel(f"Observed ({unit})")
        fig.suptitle("Official first-order mimo_srs · Matched external 2024 comparison\n"
                     f"One versus three completed rainfall-hour inputs · N = {len(reference):,}", fontsize=14)
        fig.savefig(paths["observed_vs_predicted"], dpi=180)
        plt.close(fig)

        fig, axes = plt.subplots(4, 1, figsize=(13.0, 10.4), sharex=True, layout="constrained")
        for ax, (key, title, unit) in zip(axes, VARIABLES):
            plot_series(ax, key, window_index)
            ax.set_ylabel(f"{title}\n({unit})")
        axes[0].legend(ncol=4, frameon=False, loc="upper left")
        time_axis(axes[-1])
        window_label = "Longest eligible contiguous run" if window_info["fallback"] else "Earliest eligible 168 consecutive hours"
        fig.suptitle(f"External 2024 continuous example · {window_label}\n"
                     "Retrospective conditional estimates; actual target supports differ by variable", fontsize=13)
        fig.savefig(paths["continuous_example"], dpi=180)
        plt.close(fig)

        fig = plt.figure(figsize=(13.0, 9.8), layout="constrained")
        grid = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05])
        hist_ax, acf_ax = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])
        event_ax = fig.add_subplot(grid[1, :])
        residuals = {model: models[model]["P_observed"] - models[model]["P_predicted"]
                     for model in MODELS[:2]}
        combined = np.concatenate([value.to_numpy() for value in residuals.values()])
        bins = np.histogram_bin_edges(combined, bins=70)
        for model in MODELS[:2]:
            hist_ax.hist(residuals[model], bins=bins, histtype="step", linewidth=1.6,
                         color=COLORS[model], label=LABELS[model])
            acf_ax.plot([item["lag_hours"] for item in acfs[model]],
                        [np.nan if item["autocorrelation"] is None else item["autocorrelation"]
                         for item in acfs[model]], color=COLORS[model],
                        label=LABELS[model], linewidth=1.4)
        hist_ax.axvline(0, color=COLORS["observed"], linestyle="--", linewidth=0.8)
        hist_ax.set(yscale="log", xlabel="Observed − predicted rainfall (mm/h)",
                    ylabel="Number of samples (log scale)", title="Rainfall residual distribution")
        hist_ax.legend(frameon=False)
        acf_ax.axhline(0, color=COLORS["observed"], linewidth=0.8)
        acf_ax.set(xlim=(0, 49), ylim=(-1.05, 1.05), xlabel="Exact UTC lag (hours)",
                   ylabel="Pairwise residual correlation", title="Rainfall residual autocorrelation")
        acf_ax.set_xticks([1, 12, 24, 36, 48])
        acf_ax.legend(frameon=False)
        plot_series(event_ax, "P", peak_index)
        event_ax.axvline(peak_time, color="#555555", linewidth=0.8, linestyle=":")
        event_ax.set(ylabel="Precipitation (mm/h)",
                     title="Largest observed rainfall in matched 2024 samples · ±24-hour detail")
        event_ax.legend(ncol=4, frameon=False, loc="upper right")
        time_axis(event_ax)
        fig.suptitle("External 2024 rainfall diagnostics\n"
                     "Raw predictions retained; missing hours remain gaps; ACF is descriptive", fontsize=14)
        fig.savefig(paths["rainfall_diagnostics"], dpi=180)
        plt.close(fig)

    return {
        "figures": {key: str(path.resolve()) for key, path in paths.items()},
        "external_common_count": int(len(reference)),
        "comparison_note": "The one-hour comparator is a controlled refit; historical outputs remain separate.",
        "continuous_example": window_info,
        "maximum_rainfall_example": {
            "selection_rule": "Largest observed rainfall among common external 2024 timestamps; earliest timestamp breaks ties. Display ±24 UTC hours without filling gaps.",
            "timestamp_utc": peak_time.isoformat(),
            "observed_mm_per_hour": float(reference.loc[peak_time, "P_observed"]),
            "predicted_mm_per_hour": {model: float(models[model].loc[peak_time, "P_predicted"])
                                      for model in MODELS},
            "available_hours_in_display": int(reference.index.isin(peak_index).sum()),
        },
        "rainfall_residual_acf": {
            "definition": "Pairwise Pearson correlation of observed-minus-predicted rainfall residuals separated by exact UTC hours. No compression or filling; no whiteness test is implied.",
            "by_model": acfs,
        },
    }
