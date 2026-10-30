"""Conditional marginal distributions and reproducible mixed-variable PIT.

Candidate comparison uses an inner chronological validation quantile loss. All
candidates, including the empirical distribution, use the identical score.
"""
from __future__ import annotations

import numpy as np
from scipy import stats

QUANTILES = np.array([.01, .05, .1, .25, .5, .75, .9, .95, .99])


def _fit(values, name):
    x = np.asarray(values, dtype=float)
    if not len(x):
        raise ValueError("No observations available for marginal fitting")
    if name == "empirical":
        return {"name": name, "values": np.sort(x).tolist()}
    if np.ptp(x) < 1e-12:
        raise ValueError("A continuous parametric distribution needs nonzero variance")
    distribution = getattr(stats, name)
    args = distribution.fit(x) if name == "norm" else distribution.fit(x, floc=0)
    if not np.all(np.isfinite(args)) or args[-1] <= 0:
        raise ValueError("Invalid fitted distribution parameters")
    return {"name": name, "parameters": list(map(float, args))}


def ppf(model, probability):
    u = np.clip(np.asarray(probability, dtype=float), 1e-7, 1 - 1e-7)
    if model["name"] == "empirical":
        values = np.asarray(model["values"])
        return values[np.minimum((u * len(values)).astype(int), len(values) - 1)]
    return getattr(stats, model["name"]).ppf(u, *model["parameters"])


def cdf_interval(model, values):
    """Return [F(x-),F(x)]; empirical tied observations have probability mass."""
    x = np.asarray(values, dtype=float)
    if model["name"] == "empirical":
        samples = np.asarray(model["values"])
        return (np.searchsorted(samples, x, side="left") / len(samples),
                np.searchsorted(samples, x, side="right") / len(samples))
    cdf = getattr(stats, model["name"]).cdf(x, *model["parameters"])
    return cdf, cdf


def fit_marginal(values, candidates, mixed_zero=False):
    """Fit on chronological data, score 70%/30% inner split, then refit all.

    The conditional positive-amount comparison excludes zeros identically for
    every candidate. The zero probability is separately estimated from all
    valid training observations. Inner validation is not the held-out year.
    """
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if not len(x):
        raise ValueError("Empty training marginal")
    if mixed_zero and np.any(x < 0):
        raise ValueError("Negative observation in a nonnegative marginal")
    zero_probability = float(np.mean(x == 0)) if mixed_zero else 0.0
    amount = x[x > 0] if mixed_zero else x
    if not len(amount):
        return {"zero_probability": 1.0, "distribution": None, "sample_count": len(x),
                "positive_count": 0, "candidates": [], "selection": "all_zero_training"}
    split = min(len(amount) - 1, max(1, int(.7 * len(amount))))
    inner_train, inner_validation = amount[:split], amount[split:]
    scores = []
    for name in candidates:
        try:
            candidate = _fit(inner_train, name)
            predicted = ppf(candidate, QUANTILES)
            residual = inner_validation[:, None] - predicted[None, :]
            loss = float(np.mean(np.maximum(QUANTILES * residual, (QUANTILES - 1) * residual)))
            if not np.isfinite(loss):
                raise ValueError("Non-finite quantile loss")
            scores.append({"name": name, "mean_quantile_loss_physical_units": loss,
                           "inner_train_count": len(inner_train), "inner_validation_count": len(inner_validation),
                           "inner_fit": candidate if name != "empirical" else {"name": name, "sample_count": len(inner_train)},
                           "status": "ok"})
        except (ValueError, RuntimeError, FloatingPointError) as exc:
            scores.append({"name": name, "status": "not_fitted", "reason": str(exc)})
    successful = [item for item in scores if item["status"] == "ok"]
    if not successful:
        raise ValueError("No marginal candidate could be fitted")
    best = min(successful, key=lambda item: item["mean_quantile_loss_physical_units"])
    return {"zero_probability": zero_probability, "distribution": _fit(amount, best["name"]),
            "sample_count": len(x), "positive_count": len(amount), "candidates": scores,
            "selection": "minimum_inner_chronological_mean_quantile_loss",
            "score_quantiles": QUANTILES.tolist()}


def randomized_pit(model, values, rng):
    x = np.asarray(values, dtype=float)
    out = np.full(x.shape, np.nan)
    valid = np.isfinite(x)
    p0 = model["zero_probability"]
    zero = valid & (x == 0) if p0 > 0 else np.zeros(x.shape, dtype=bool)
    out[zero] = rng.uniform(0, p0, int(zero.sum()))
    remaining = valid & ~zero
    if np.any(remaining):
        if model["distribution"] is None:
            raise ValueError("Positive validation value under all-zero fitted marginal")
        lower, upper = cdf_interval(model["distribution"], x[remaining])
        out[remaining] = p0 + (1 - p0) * (lower + rng.random(lower.shape) * (upper - lower))
    return np.clip(out, 1e-7, 1 - 1e-7)


def inverse_marginal(model, probability):
    u = np.asarray(probability, dtype=float)
    p0 = model["zero_probability"]
    out = np.zeros(u.shape)
    positive = u > p0
    if np.any(positive) and model["distribution"] is not None:
        out[positive] = ppf(model["distribution"], (u[positive] - p0) / (1 - p0))
    return out
