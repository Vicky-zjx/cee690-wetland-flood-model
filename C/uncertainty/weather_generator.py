"""Seasonal mixed marginals, Gaussian-copula VAR(1), and joint block bootstrap."""
from __future__ import annotations

import hashlib
import json
import numpy as np
import pandas as pd
from scipy import linalg, stats

from .marginals import fit_marginal, inverse_marginal, randomized_pit

COLUMNS = ["precipitation_mm_per_h", "temperature_C", "solar_W_per_m2", "wind_m_per_s"]
SEASONS = ["DJF", "MAM", "JJA", "SON"]
DEFAULTS = {"latitude": 43.0547222, "longitude": -88.0461111,
            "timezone": "America/Chicago", "pit_seed": 1729, "burn_in_hours": 336,
            "block_hours": 48, "bootstrap_month_radius": 1,
            "minimum_group_samples": 80, "max_var_radius": .98}


def _settings(config):
    return {**DEFAULTS, **(config.get("weather", config) if config else {})}


def utc_index(values):
    index = pd.DatetimeIndex(pd.to_datetime(values))
    if index.tz is None:
        raise ValueError("Timestamps must explicitly be timezone-aware")
    return index.tz_convert("UTC")


def _time_features(index, timezone):
    local = (index + pd.Timedelta(minutes=30)).tz_convert(timezone)
    annual = 2 * np.pi * (local.dayofyear.to_numpy() - 1) / 365.2425
    daily = 2 * np.pi * (local.hour.to_numpy() + local.minute.to_numpy() / 60) / 24
    return np.column_stack([np.ones(len(index)), np.sin(annual), np.cos(annual),
                            np.sin(2 * annual), np.cos(2 * annual), np.sin(daily),
                            np.cos(daily), np.sin(2 * daily), np.cos(2 * daily),
                            np.sin(annual) * np.sin(daily), np.cos(annual) * np.cos(daily)])


def solar_geometry(index, latitude, longitude):
    """NOAA fractional-year approximate solar elevation, integrated hourly.

    The mean positive sine of elevation across 5-minute subinterval midpoints
    scales daylight radiation. A bin is night only if all subintervals are dark.
    This is astronomical conditioning, not a calibrated irradiance model.
    """
    index = utc_index(index)
    elevation = []
    lat = np.deg2rad(latitude)
    for offset in np.arange(2.5, 60, 5):
        t = index + pd.Timedelta(minutes=float(offset))
        hour = t.hour.to_numpy() + t.minute.to_numpy() / 60 + t.second.to_numpy() / 3600
        gamma = 2 * np.pi / 365 * (t.dayofyear.to_numpy() - 1 + (hour - 12) / 24)
        equation = 229.18 * (.000075 + .001868 * np.cos(gamma) - .032077 * np.sin(gamma)
                             - .014615 * np.cos(2 * gamma) - .040849 * np.sin(2 * gamma))
        declination = (.006918 - .399912 * np.cos(gamma) + .070257 * np.sin(gamma)
                       - .006758 * np.cos(2 * gamma) + .000907 * np.sin(2 * gamma)
                       - .002697 * np.cos(3 * gamma) + .00148 * np.sin(3 * gamma))
        angle = np.deg2rad((hour * 60 + equation + 4 * longitude) / 4 - 180)
        elevation.append(np.sin(lat) * np.sin(declination) + np.cos(lat) * np.cos(declination) * np.cos(angle))
    positive = np.maximum(np.asarray(elevation), 0)
    return positive.mean(axis=0), positive.max(axis=0) > 0


def _season(index, timezone):
    month = (index + pd.Timedelta(minutes=30)).tz_convert(timezone).month.to_numpy()
    return ((month % 12) // 3).astype(int)


def _repair(matrix, correlation=False, minimum=1e-6):
    symmetric = (matrix + matrix.T) / 2
    eig, vec = np.linalg.eigh(symmetric)
    adjusted = (vec * np.maximum(eig, minimum)) @ vec.T
    if correlation:
        scale = np.sqrt(np.diag(adjusted))
        adjusted /= np.outer(scale, scale)
        np.fill_diagonal(adjusted, 1.0)
    return adjusted, {"original": matrix.tolist(), "adjusted": adjusted.tolist(),
                      "original_eigenvalues": eig.tolist(),
                      "adjusted_eigenvalues": np.linalg.eigvalsh(adjusted).tolist(),
                      "frobenius_adjustment": float(np.linalg.norm(adjusted - matrix)),
                      "max_absolute_adjustment": float(np.max(np.abs(adjusted - matrix))),
                      "method": "symmetric_eigenvalue_floor_then_diagonal_normalization" if correlation else "symmetric_eigenvalue_floor",
                      "eigenvalue_floor": minimum}


def fit_weather(train_df, config):
    """Use training records only; return a JSON-serializable fitted model."""
    settings = _settings(config)
    index = utc_index(train_df["timestamp_utc"])
    if not index.is_monotonic_increasing or index.has_duplicates:
        raise ValueError("Training timestamps must be unique and increasing")
    x = train_df[COLUMNS].to_numpy(dtype=float)
    if np.any(x[:, [0, 2, 3]][np.isfinite(x[:, [0, 2, 3]])] < 0):
        raise ValueError("Nonnegative weather variable has negative training values")
    design = _time_features(index, settings["timezone"])
    temperature_valid = np.isfinite(x[:, 1])
    if temperature_valid.sum() < design.shape[1] * 3:
        raise ValueError("Insufficient training temperature observations")
    beta = np.linalg.lstsq(design[temperature_valid], x[temperature_valid, 1], rcond=None)[0]
    solar_scale, daylight = solar_geometry(index, settings["latitude"], settings["longitude"])
    transformed = x.copy()
    transformed[:, 1] -= design @ beta
    # Floor only the deterministic divisor in twilight; never alter observations.
    transformed[:, 2] = np.where(daylight, x[:, 2] / np.maximum(solar_scale, .05), np.nan)
    season = _season(index, settings["timezone"])
    for j, variable in enumerate(COLUMNS):
        values = transformed[np.isfinite(transformed[:, j]), j]
        if not len(values) or np.ptp(values) < 1e-10:
            raise ValueError(f"Zero-variance or unavailable training variable: {variable}; do not standardize it")
    marginals, latent = {}, np.full(x.shape, np.nan)
    rng = np.random.default_rng(settings["pit_seed"])
    candidate_sets = [["gamma", "lognorm", "empirical"], ["norm", "empirical"],
                      ["gamma", "lognorm", "empirical"], ["weibull_min", "gamma", "empirical"]]
    for k, label in enumerate(SEASONS):
        marginals[label] = {}
        for j, variable in enumerate(COLUMNS):
            selection = season == k
            fitting = selection & np.isfinite(transformed[:, j])
            pooled = False
            if fitting.sum() < settings["minimum_group_samples"]:
                fitting = np.isfinite(transformed[:, j])
                pooled = True
            if fitting.sum() < 20:
                raise ValueError(f"Insufficient observations for {label}/{variable}")
            marginal = fit_marginal(transformed[fitting, j], candidate_sets[j], mixed_zero=j != 1)
            marginal["pooled_all_seasons"] = pooled
            marginals[label][variable] = marginal
            latent[selection, j] = stats.norm.ppf(randomized_pit(marginal, transformed[selection, j], rng))
    centers, scales = np.nanmean(latent, axis=0), np.nanstd(latent, axis=0)
    if np.any(~np.isfinite(scales)) or np.any(scales < 1e-8):
        raise ValueError("A latent variable has zero or undefined variance")
    standardized = (latent - centers) / scales
    raw_corr = pd.DataFrame(standardized).corr(min_periods=100).to_numpy()
    if not np.all(np.isfinite(raw_corr)):
        raise ValueError("Insufficient pairwise overlap to estimate latent dependence")
    c0, dependence = _repair(raw_corr, correlation=True)
    # A covariance-moment VAR permits observed nighttime P/T/W to contribute
    # without imputing unobserved nighttime cloudiness as data.
    adjacent = np.diff(index.asi8) == pd.Timedelta(hours=1).value
    lag_cov, pair_counts = np.zeros((4, 4)), np.zeros((4, 4), dtype=int)
    for j in range(4):
        for k in range(4):
            good = adjacent & np.isfinite(standardized[1:, j]) & np.isfinite(standardized[:-1, k])
            pair_counts[j, k] = int(good.sum())
            if good.sum() < 100:
                raise ValueError("Insufficient truly adjacent pairs for VAR estimation")
            lag_cov[j, k] = np.mean(standardized[1:, j][good] * standardized[:-1, k][good])
    original_a = lag_cov @ np.linalg.inv(c0)
    original_radius = float(np.max(np.abs(np.linalg.eigvals(original_a))))
    maximum = float(settings["max_var_radius"])
    if not 0 < maximum < 1:
        raise ValueError("max_var_radius must lie strictly between zero and one")
    a = original_a * min(1.0, maximum / max(original_radius, 1e-15))
    q_original = c0 - a @ c0 @ a.T
    q, innovation = _repair(q_original)
    # Cholesky applies ONCE to innovation covariance Q, not again to C0.
    chol = np.linalg.cholesky(q)
    stationary_cov = linalg.solve_discrete_lyapunov(a, q)
    stationary_scale = np.sqrt(np.diag(stationary_cov))
    # Residual diagnostics only where the observed current and all four prior
    # latent variables exist. Nighttime excluded; counts explicitly reported.
    residual_good = adjacent & np.isfinite(standardized[1:]).all(axis=1) & np.isfinite(standardized[:-1]).all(axis=1)
    residual = standardized[1:][residual_good] - standardized[:-1][residual_good] @ a.T
    residual_times = index[1:][residual_good]
    residual_lag1 = []
    pairs = np.diff(residual_times.asi8) == pd.Timedelta(hours=1).value
    for j in range(4):
        residual_lag1.append(float(np.corrcoef(residual[1:, j][pairs], residual[:-1, j][pairs])[0, 1]) if pairs.sum() > 10 else None)
    persistence = {}
    for lag in sorted({1, 6, 12, 24, 48, 72, int(settings["block_hours"])}):
        is_adjacent = (index[lag:].asi8 - index[:-lag].asi8) == pd.Timedelta(hours=lag).value
        correlations = []
        for j in range(4):
            good = is_adjacent & np.isfinite(standardized[lag:, j]) & np.isfinite(standardized[:-lag, j])
            correlations.append(float(np.corrcoef(standardized[lag:, j][good], standardized[:-lag, j][good])[0, 1]) if good.sum() > 10 else None)
        persistence[str(lag)] = correlations
    block_acf = [abs(v) for v in persistence[str(int(settings["block_hours"]))] if v is not None]
    block_adequate = bool(block_acf and max(block_acf) < .2)
    dataframe_hash = hashlib.sha256(pd.util.hash_pandas_object(train_df[["timestamp_utc"] + COLUMNS], index=False).values.tobytes()).hexdigest()
    fit = {"schema_version": 1, "training_start": index.min().isoformat(), "training_end": index.max().isoformat(),
           "training_row_count": len(index), "training_data_hash": dataframe_hash, "settings": settings,
           "marginals": marginals, "temperature_harmonic_coefficients": beta.tolist(),
           "temperature_feature_order": "constant, annual sin/cos, semiannual sin/cos, diurnal sin/cos, semidiurnal sin/cos, annual_sin*diurnal_sin, annual_cos*diurnal_cos",
           "solar_transform": "daylight radiation / max(hourly mean positive sine solar elevation,0.05)",
           "nighttime_solar_training_rows_excluded": int((~daylight).sum()),
           "observed_night_radiation_positive_count": int(np.sum((~daylight) & np.isfinite(x[:, 2]) & (x[:, 2] > 1))),
           "latent_centers": centers.tolist(), "latent_scales": scales.tolist(), "dependence": dependence,
           "var": {"estimator": "pairwise contiguous-hour covariance moments; A=C1 inv(C0)",
                   "A_original": original_a.tolist(), "A": a.tolist(), "intercept": [0.] * 4,
                   "spectral_radius_original": original_radius, "spectral_radius": float(np.max(np.abs(np.linalg.eigvals(a)))),
                   "lag_pair_counts": pair_counts.tolist(), "Q": q.tolist(), "Q_repair": innovation,
                   "innovation_cholesky": chol.tolist(), "stationary_mean": [0.] * 4,
                   "stationary_covariance": stationary_cov.tolist(), "stationary_scale": stationary_scale.tolist(),
                   "residual_complete_daylight_count": int(len(residual)), "residual_lag1": residual_lag1},
           "latent_acf": persistence,
           "bootstrap_block_diagnostics": {"block_hours": int(settings["block_hours"]), "latent_acf_at_block": persistence[str(int(settings["block_hours"]))], "max_absolute_latent_acf_at_block": max(block_acf) if block_acf else None, "acf_below_prespecified_point2_screen": block_adequate},
           "bootstrap_block_basis": "48 h preregistered benchmark covers two diurnal cycles. Training ACF at block length is checked against absolute 0.2; failure explicitly means remaining persistence exceeds this block, not permission to tune after validation. Block boundaries truncate dependence."}
    fit["model_hash"] = hashlib.sha256(json.dumps(fit, sort_keys=True).encode()).hexdigest()
    return fit


def generate_weather(fit, timestamps, rng=None, seed=None, method="copula_var", training_df=None, config=None):
    """Generate a deterministic scenario given its independently assigned RNG."""
    if rng is not None and seed is not None:
        raise ValueError("Pass rng OR seed, not both")
    rng = np.random.default_rng(seed) if rng is None else rng
    settings = {**fit["settings"], **((config or {}).get("weather", {}))}
    index = utc_index(timestamps)
    if len(index) > 1 and not np.all(np.diff(index.asi8) == pd.Timedelta(hours=1).value):
        raise ValueError("Generated weather grid must contain consecutive hourly interval starts")
    if method in {"block_bootstrap", "bootstrap"}:
        if training_df is None:
            raise ValueError("Block bootstrap requires original training observations")
        values, diagnostics = _bootstrap(training_df, index, rng, settings)
    elif method in {"copula_var", "var"}:
        a, chol = np.asarray(fit["var"]["A"]), np.asarray(fit["var"]["innovation_cholesky"])
        burn = int(settings["burn_in_hours"])
        if burn < 1:
            raise ValueError("A positive burn-in is required")
        z, states = np.zeros(4), np.empty((len(index), 4))
        for i in range(burn + len(index)):
            z = a @ z + chol @ rng.standard_normal(4)
            if i >= burn:
                states[i - burn] = z
        # Fixed fitted stationary variance, NEVER a scenario-specific rescaling.
        u = stats.norm.cdf(states / np.asarray(fit["var"]["stationary_scale"]))
        season, values = _season(index, settings["timezone"]), np.empty(states.shape)
        for k, label in enumerate(SEASONS):
            mask = season == k
            for j, variable in enumerate(COLUMNS):
                values[mask, j] = inverse_marginal(fit["marginals"][label][variable], u[mask, j])
        values[:, 1] += _time_features(index, settings["timezone"]) @ np.asarray(fit["temperature_harmonic_coefficients"])
        solar_scale, daylight = solar_geometry(index, settings["latitude"], settings["longitude"])
        values[:, 2] *= np.maximum(solar_scale, .05)
        values[~daylight, 2] = 0.0
        diagnostics = {"burn_in_hours": burn, "stationary_normalization": "fitted_process_not_scenario"}
    else:
        raise ValueError(f"Unknown generation method: {method}")
    out = pd.DataFrame(values, columns=COLUMNS)
    out.insert(0, "timestamp_utc", index)
    out.attrs["generation_diagnostics"] = diagnostics
    return out


def _bootstrap(training_df, target_index, rng, settings):
    train_index = utc_index(training_df["timestamp_utc"])
    x = training_df[COLUMNS].to_numpy(dtype=float)
    block = int(settings["block_hours"])
    if block < 1:
        raise ValueError("block_hours must be positive")
    local = train_index.tz_convert(settings["timezone"])
    target_local = target_index.tz_convert(settings["timezone"])
    valid_row = np.isfinite(x).all(axis=1)
    # Full blocks only: no interpolation or hidden removal of missing values.
    bad = (~valid_row).astype(int)
    cs = np.r_[0, np.cumsum(bad)]
    starts = np.arange(max(0, len(x) - block + 1))
    eligible = (cs[starts + block] - cs[starts]) == 0
    if block > 1:
        eligible &= (train_index[starts + block - 1].asi8 - train_index[starts].asi8) == pd.Timedelta(hours=block - 1).value
    starts = starts[eligible]
    result, sources = np.empty((len(target_index), 4)), []
    for offset in range(0, len(target_index), block):
        target = target_local[offset]
        month_difference = np.abs(local.month.to_numpy()[starts] - target.month)
        month_difference = np.minimum(month_difference, 12 - month_difference)
        matched = starts[(local.hour.to_numpy()[starts] == target.hour) & (month_difference <= settings["bootstrap_month_radius"])]
        if not len(matched):
            raise ValueError(f"No complete {block}-hour training block aligned to month/hour {target.month}/{target.hour}")
        start = int(rng.choice(matched))
        count = min(block, len(target_index) - offset)
        result[offset:offset + count] = x[start:start + count]
        sources.append({"target_offset_hours": offset, "source_start_utc": train_index[start].isoformat(), "hours": count,
                        "eligible_matched_blocks": len(matched)})
    _, daylight = solar_geometry(target_index, settings["latitude"], settings["longitude"])
    changed = (~daylight) & (result[:, 2] != 0)
    removed_energy = float(np.sum(result[changed, 2]))
    result[~daylight, 2] = 0
    return result, {"block_hours": block, "source_blocks": sources,
                    "night_solar_mask_altered_hours": int(changed.sum()), "night_solar_removed_Wh_m2": removed_energy,
                    "night_solar_note": "Joint observed blocks retained except deterministic target-night solar mask; source/target calendar daylight differs"}
