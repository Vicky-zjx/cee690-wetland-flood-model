"""Diagnostic figures for the fixed-feature official FIT initial-L1 comparison.

The caller fits and evaluates all models. This module does not fit, clip,
rescale, impute, discard samples, or choose a model. Missing display hours are
shown as gaps. Initial L1 is distinguished from the official adaptive alpha.
"""
from __future__ import annotations

from datetime import timezone
from pathlib import Path

import numpy as np
import pandas as pd


PENALTIES = {"l1_1": 1.0, "l1_0p1": 0.1, "l1_0p01": 0.01}
MODELS = tuple(PENALTIES) + ("one_hour", "persistence")
PERIODS = ("validation_2023", "external_2024")
COLORS = {"l1_1": "#007A9F", "l1_0p1": "#C25B27", "l1_0p01": "#22835C",
          "one_hour": "#7B60A0", "persistence": "#777777", "observed": "#20262D"}
LABELS = {"l1_1": "Initial L1 = 1", "l1_0p1": "Initial L1 = 0.1",
          "l1_0p01": "Initial L1 = 0.01", "one_hour": "One-hour linear",
          "persistence": "Persistence"}
STYLES = {"l1_1": "-", "l1_0p1": "--", "l1_0p01": "-.",
          "one_hour": ":", "persistence": "--"}


def _require(frame, columns, name):
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"{name}: missing columns {missing}")


def _validated(predictions, metrics, iterations, diagnostics):
    values = [f"{key}_{kind}" for key in ("P", "T", "Rs", "W")
              for kind in ("observed", "predicted")]
    _require(predictions, ["period", "model", "timestamp_utc", "prediction_cutoff_utc"] + values,
             "predictions")
    _require(metrics, ["period", "model", "variable", "N", "MAE", "RMSE",
                       "mean_prediction_error", "R_squared"], "metrics")
    _require(iterations, ["model", "output", "iteration", "alpha_used", "alpha_next",
                          "accepted", "candidate_square", "c_square_before_zeroing", "c_square_returned",
                          "fit_error_candidate"], "iterations")
    _require(diagnostics, ["model", "output", "initial_penalty", "square_coefficient_scaled",
                           "kkt_max_abs_residual", "last_accepted_alpha", "termination"],
             "diagnostics")
    predictions = predictions.copy()
    for field in ("timestamp_utc", "prediction_cutoff_utc"):
        predictions[field] = pd.to_datetime(predictions[field], utc=True)
    periods = {}
    for period in PERIODS:
        frame = predictions.loc[predictions.period.eq(period)].copy()
        if frame.empty or set(frame.model.unique()) != set(MODELS):
            raise ValueError(f"{period}: require all five specified model labels")
        if frame[["timestamp_utc", "prediction_cutoff_utc"]].isna().any().any():
            raise ValueError(f"{period}: missing timestamp or cutoff")
        year = int(period.rsplit("_", 1)[1])
        if not frame.timestamp_utc.dt.year.eq(year).all():
            raise ValueError(f"{period}: target reference timestamp in another year")
        if not (frame.timestamp_utc-frame.prediction_cutoff_utc).eq(pd.Timedelta(hours=1)).all():
            raise ValueError(f"{period}: reference timestamp must equal cutoff plus one hour")
        if frame.duplicated(["model", "timestamp_utc"]).any():
            raise ValueError(f"{period}: duplicate target timestamps")
        if not np.isfinite(frame[values].to_numpy(dtype=float)).all():
            raise ValueError(f"{period}: nonfinite values; no samples were discarded")
        models = {model: frame.loc[frame.model.eq(model)].sort_values("timestamp_utc")
                  .set_index("timestamp_utc") for model in MODELS}
        reference = models["one_hour"]
        observed = [f"{key}_observed" for key in ("P", "T", "Rs", "W")]
        for model, candidate in models.items():
            if not candidate.index.equals(reference.index):
                raise ValueError(f"{period}/{model}: unmatched target rows")
            if not candidate.prediction_cutoff_utc.equals(reference.prediction_cutoff_utc):
                raise ValueError(f"{period}/{model}: unmatched cutoffs")
            if not np.array_equal(candidate[observed].to_numpy(), reference[observed].to_numpy()):
                raise ValueError(f"{period}/{model}: observations differ")
            metric = metrics.loc[metrics.period.eq(period) & metrics.model.eq(model)
                                 & metrics.variable.eq("P")]
            if len(metric) != 1:
                raise ValueError(f"{period}/{model}: require one rainfall metric row")
            row = metric.iloc[0]
            if not np.isfinite(row[["N", "MAE", "RMSE", "mean_prediction_error", "R_squared"]]
                               .to_numpy(dtype=float)).all() or int(row.N) != len(reference):
                raise ValueError(f"{period}/{model}: missing metrics or inconsistent N")
            residual = candidate.P_predicted.to_numpy()-candidate.P_observed.to_numpy()
            if not np.allclose([row.MAE, row.RMSE],
                               [np.abs(residual).mean(), np.sqrt(np.mean(residual**2))],
                               rtol=1e-9, atol=1e-12):
                raise ValueError(f"{period}/{model}: rainfall metrics do not match predictions")
        periods[period] = models
    histories, summaries = {}, {}
    for model, initial in PENALTIES.items():
        history = iterations.loc[iterations.model.eq(model) & iterations.output.eq("P")].copy()
        summary = diagnostics.loc[diagnostics.model.eq(model) & diagnostics.output.eq("P")]
        if history.empty or len(summary) != 1:
            raise ValueError(f"{model}: missing rainfall iterations or diagnostic row")
        history = history.sort_values("iteration")
        if history.iteration.duplicated().any():
            raise ValueError(f"{model}: duplicate rainfall iteration numbers")
        required_numeric = ["iteration", "alpha_used", "candidate_square", "c_square_before_zeroing",
                            "c_square_returned", "fit_error_candidate"]
        if not np.isfinite(history[required_numeric].to_numpy(dtype=float)).all():
            raise ValueError(f"{model}: nonfinite rainfall iteration diagnostics")
        if history.alpha_used.le(0).any():
            raise ValueError(f"{model}: adaptive alpha must be positive for log display")
        accepted = history.accepted
        if accepted.dtype != bool:
            accepted = accepted.map({True: True, False: False, 1: True, 0: False,
                                     "True": True, "False": False, "true": True, "false": False})
        if accepted.isna().any():
            raise ValueError(f"{model}: accepted must contain explicit booleans")
        history["accepted"] = accepted.astype(bool)
        row = summary.iloc[0]
        if not np.isfinite([float(row.initial_penalty), float(row.square_coefficient_scaled),
                            float(row.kkt_max_abs_residual)]).all():
            raise ValueError(f"{model}: missing numerical rainfall diagnostic")
        if not np.isclose(float(row.initial_penalty), initial) or pd.isna(row.termination):
            raise ValueError(f"{model}: invalid initial penalty or missing termination")
        if not np.isclose(history.iloc[-1].c_square_returned, row.square_coefficient_scaled,
                          rtol=1e-8, atol=1e-10):
            raise ValueError(f"{model}: returned coefficient and final diagnostic differ")
        histories[model], summaries[model] = history, row
    return periods, histories, summaries


def make_figures(predictions, metrics, iterations, diagnostics, outdir) -> dict:
    """Save two matched comparison PNGs; return paths and selection metadata.

    Predictions and metrics include validation_2023 and external_2024 for all
    five MODELS. Iterations and diagnostics use output='P' for rainfall and
    include each of the three initial penalties. All inputs remain unchanged.
    """
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    periods, histories, summaries = _validated(predictions, metrics, iterations, diagnostics)
    external = periods["external_2024"]
    reference = external["one_hour"]
    peak_time = reference.P_observed.idxmax()  # Sorted timestamps: earliest tie wins.
    peak_window = pd.date_range(peak_time-pd.Timedelta(hours=24),
                                peak_time+pd.Timedelta(hours=24), freq="h")
    outdir = Path(outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    paths = {name: outdir / f"{name}.png" for name in
             ("penalty_comparison", "rainfall_and_iterations")}
    counts = {model: {"negative": int(external[model].P_predicted.lt(0).sum()),
                       "exact_zero": int(external[model].P_predicted.eq(0).sum())}
              for model in PENALTIES}
    style = {"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
             "axes.spines.top": False, "axes.spines.right": False,
             "axes.grid": True, "grid.alpha": 0.18,
             "figure.facecolor": "white", "savefig.facecolor": "white"}
    with plt.rc_context(style):
        fig, axes = plt.subplots(2, 2, figsize=(12.6, 8.6), layout="constrained")
        width = 0.23
        for ax, metric in zip(axes[0], ("RMSE", "MAE")):
            for offset, model in enumerate(PENALTIES):
                y = [float(metrics.loc[metrics.period.eq(period) & metrics.model.eq(model)
                                       & metrics.variable.eq("P"), metric].iloc[0])
                     for period in PERIODS]
                bars = ax.bar(np.arange(2)+(offset-1)*width, y, width,
                              color=COLORS[model], label=LABELS[model], alpha=0.87)
                ax.bar_label(bars, fmt="%.4f", padding=3, fontsize=9, rotation=0)
            ax.set(xticks=[0, 1], xticklabels=["2023 validation", "2024 evaluation"],
                   ylabel=f"Rainfall {metric} (mm/h)", title=f"Rainfall {metric} · lower is better")
            ax.margins(y=0.22)
            ax.set_ylim(bottom=0)
        axes[0, 0].legend(frameon=False, fontsize=9, loc="upper left")
        initial_labels = ["1", "0.1", "0.01"]
        coefficients = [float(summaries[model].square_coefficient_scaled) for model in PENALTIES]
        coefficient_ax = axes[1, 0]
        bars = coefficient_ax.bar(np.arange(3), coefficients,
                                  color=[COLORS[model] for model in PENALTIES], alpha=0.87)
        coefficient_ax.bar_label(bars, labels=[f"{value:.6g}" for value in coefficients],
                                 padding=4, fontsize=10)
        coefficient_ax.axhline(0, color="#555555", linewidth=0.8)
        coefficient_ax.set(xticks=np.arange(3), xticklabels=initial_labels,
                           xlabel="Supplied initial L1 penalty",
                           ylabel="Final square-feature coefficient (scaled units)",
                           title="Returned rainfall square coefficient")
        coefficient_ax.margins(y=0.3)
        # Keep an exact-zero label inside the panel when other coefficients are negative.
        coefficient_span = max(max(coefficients + [0.0])-min(coefficients + [0.0]), 0.01)
        coefficient_ax.set_ylim(min(coefficients + [0.0])-0.28*coefficient_span,
                                max(coefficients + [0.0])+0.16*coefficient_span)
        count_ax = axes[1, 1]
        negative = count_ax.bar(np.arange(3)-0.18, [counts[m]["negative"] for m in PENALTIES],
                                0.36, color="#A84B4B", label="Negative predictions")
        zero = count_ax.bar(np.arange(3)+0.18, [counts[m]["exact_zero"] for m in PENALTIES],
                            0.36, color="#73838C", label="Exactly zero predictions")
        for bars in (negative, zero):
            count_ax.bar_label(bars, fmt="%d", padding=3, fontsize=9)
        count_ax.set(xticks=np.arange(3), xticklabels=initial_labels,
                     xlabel="Supplied initial L1 penalty", ylabel="2024 prediction count",
                     title=f"Raw rainfall predictions · N = {len(reference):,}")
        count_ax.margins(y=0.3)
        count_ax.legend(frameon=False, fontsize=9, loc="upper left")
        fig.suptitle("Official mimo_srs · Fixed inputs with one squared-rain feature\n"
                     "Initial L1 comparison; the official algorithm adapts its penalty", fontsize=14)
        fig.savefig(paths["penalty_comparison"], dpi=180)
        plt.close(fig)

        fig = plt.figure(figsize=(13.6, 9.3), layout="constrained")
        grid = fig.add_gridspec(2, 2, height_ratios=[1.1, 1.0])
        event_ax = fig.add_subplot(grid[0, :])
        coefficient_ax, alpha_ax = fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])
        event_ax.plot(peak_window, reference.P_observed.reindex(peak_window),
                      color=COLORS["observed"], label="Observed", linewidth=1.8, zorder=5)
        for model in MODELS:
            event_ax.plot(peak_window, external[model].P_predicted.reindex(peak_window),
                          color=COLORS[model], linestyle=STYLES[model], label=LABELS[model],
                          linewidth=1.4 if model in PENALTIES else 1.1, alpha=0.9)
        event_ax.axvline(peak_time, color="#555555", linewidth=0.8, linestyle=":")
        event_ax.set(ylabel="Precipitation (mm/h)", xlabel="Target selection/reference time (UTC)",
                     title="Largest observed rainfall among matched 2024 samples · ±24 hours; gaps retained")
        event_ax.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=5, maxticks=9, tz=timezone.utc))
        event_ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d\n%H:%M", tz=timezone.utc))
        event_ax.legend(ncol=3, frameon=False, fontsize=9, loc="upper right")
        for model, initial in PENALTIES.items():
            frame = histories[model]
            x = frame.iteration.to_numpy()
            coefficient_ax.plot(x, frame.c_square_returned, color=COLORS[model],
                                linestyle=STYLES[model], marker="o", markersize=3,
                                label=LABELS[model], linewidth=1.4)
            # Show rejected candidates separately: they did not replace the returned fit.
            rejected = ~frame.accepted
            coefficient_ax.scatter(x[rejected], frame.loc[rejected, "candidate_square"],
                                   color=COLORS[model], marker="x", s=35, zorder=5)
            alpha_ax.plot(x, frame.alpha_used, color=COLORS[model], linestyle=STYLES[model],
                          marker="o", markersize=3, label=LABELS[model], linewidth=1.4)
            alpha_ax.axhline(initial, color=COLORS[model], linestyle=":", linewidth=0.8, alpha=0.55)
        coefficient_ax.axhline(0, color="#555555", linewidth=0.8)
        coefficient_ax.set(xlabel="Official rainfall-output iteration",
                           ylabel="Square-feature coefficient (scaled units)",
                           title="Returned coefficients; × = rejected proposed coefficient")
        coefficient_ax.legend(frameon=False, fontsize=9, loc="best")
        alpha_ax.set(yscale="log", xlabel="Official rainfall-output iteration",
                     ylabel="Actual alpha_used (log scale)",
                     title="Adaptive alpha used; dotted horizontals mark initial values")
        alpha_ax.legend(frameon=False, fontsize=9, loc="best")
        for ax in (coefficient_ax, alpha_ax):
            ax.xaxis.get_major_locator().set_params(integer=True)
        fig.suptitle("Rainfall response and official iteration diagnostics\n"
                     "Fixed expanded inputs; raw predictions retained; initial penalty is not a fixed alpha", fontsize=14)
        fig.savefig(paths["rainfall_and_iterations"], dpi=180)
        plt.close(fig)

    return {"figures": {key: str(path) for key, path in paths.items()},
            "counts_by_period": {period: len(models["one_hour"]) for period, models in periods.items()},
            "initial_penalties": PENALTIES.copy(), "external_prediction_counts": counts,
            "maximum_rainfall_example": {
                "selection_rule": "Largest observed rainfall at common external_2024 timestamps; earliest tie. Display plus/minus 24 UTC hours without filling gaps.",
                "timestamp_utc": peak_time.isoformat(),
                "observed_mm_per_hour": float(reference.loc[peak_time, "P_observed"]),
                "predicted_mm_per_hour": {model: float(external[model].loc[peak_time, "P_predicted"])
                                          for model in MODELS},
                "available_hours_in_display": int(reference.index.isin(peak_window).sum())},
            "iteration_note": "alpha_used is the actual adaptive penalty. Returned square coefficients are plotted; rejected proposed candidate coefficients have x markers. All raw predictions are retained."}
