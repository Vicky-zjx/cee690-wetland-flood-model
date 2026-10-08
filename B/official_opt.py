"""Thin, audited adapter to the pinned official multivarious Nelder-Mead solver.

No optimizer or hydrologic equations are reimplemented here. The caller supplies
finite, already conditioned f and a fixed-size dimensionless residual vector.
See docs/08_official_opt_audit.md for the source quirks and feasibility limitations.
"""
from __future__ import annotations

import importlib
import inspect
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

import numpy as np

ROOT = Path(__file__).resolve().parent
OFFICIAL_CHECKOUT = ROOT.parent / "A" / "work" / "multivarious"
VERIFIED_REVISION = "06d7aed5b541ec7fe26a6dc69be6aa7fb610aa02"
DEFAULT_SETTINGS = {
    "max_evals": 200,
    "tol_x": 1e-3,
    "tol_f": 1e-3,
    "feasibility_tol": 1e-6,
    "penalty": 10.0,
    "penalty_exponent": 1.0,
    "display": 0,
    "max_seconds": None,
}


def load_official():
    """Load only the inspected checkout; reject changed tracked vendor files."""
    if not (OFFICIAL_CHECKOUT / "multivarious/opt/nms.py").is_file():
        raise ImportError(
            "B.2 requires the verified official checkout at A/work/multivarious "
            f"revision {VERIFIED_REVISION}; no dependency is installed automatically."
        )
    revision = subprocess.check_output(
        ["git", "-C", str(OFFICIAL_CHECKOUT), "rev-parse", "HEAD"], text=True
    ).strip()
    if revision != VERIFIED_REVISION:
        raise RuntimeError(f"Unreviewed multivarious revision: {revision}")
    check = subprocess.run(
        ["git", "-C", str(OFFICIAL_CHECKOUT), "diff", "--quiet", "HEAD", "--"],
        check=False,
    )
    if check.returncode != 0:
        raise RuntimeError("The official checkout has tracked edits or cannot be checked.")
    checkout_path = str(OFFICIAL_CHECKOUT)
    if checkout_path not in sys.path:
        sys.path.insert(0, checkout_path)
    # Imports must not rewrite vendor bytecode or A.2 plot/cache outputs.
    original_bytecode = sys.dont_write_bytecode
    previous_mpl = os.environ.get("MPLCONFIGDIR")
    with tempfile.TemporaryDirectory(prefix="b2-official-opt-mpl-") as cache:
        sys.dont_write_bytecode = True
        if previous_mpl is None:
            os.environ["MPLCONFIGDIR"] = cache
        try:
            package = importlib.import_module("multivarious")
            routine = importlib.import_module("multivarious.opt.nms").nms
            ode4u = importlib.import_module("multivarious.ode.ode4u").ode4u
            opt_hyp = importlib.import_module("multivarious.utl.opt_hyp").opt_hyp
        finally:
            sys.dont_write_bytecode = original_bytecode
            if previous_mpl is None:
                os.environ.pop("MPLCONFIGDIR", None)
    for obj in (routine, ode4u, opt_hyp):
        if not Path(inspect.getfile(obj)).resolve().is_relative_to(OFFICIAL_CHECKOUT.resolve()):
            raise RuntimeError("A different multivarious installation was already imported.")
    return {
        "nms": routine, "ode4u": ode4u, "opt_hyp": opt_hyp,
        "version": package.__version__, "revision": revision,
        "nms_signature": str(inspect.signature(routine)),
        "ode4u_signature": str(inspect.signature(ode4u)),
    }


def _validated_settings(settings, n):
    supplied = {} if settings is None else dict(settings)
    unknown = set(supplied) - set(DEFAULT_SETTINGS)
    if unknown:
        raise ValueError(f"Unknown OPT settings: {sorted(unknown)}")
    options = {**DEFAULT_SETTINGS, **supplied}
    max_evals = options["max_evals"]
    if not np.isfinite(max_evals) or int(max_evals) != max_evals or max_evals < n + 2:
        raise ValueError(f"max_evals must be an integer >= n+2 ({n+2}).")
    options["max_evals"] = int(max_evals)
    for name in ("tol_x", "tol_f", "penalty"):
        if not np.isfinite(options[name]) or options[name] <= 0:
            raise ValueError(f"{name} must be finite and positive.")
    if not np.isfinite(options["feasibility_tol"]) or options["feasibility_tol"] < 0:
        raise ValueError("feasibility_tol must be finite and nonnegative.")
    if not np.isfinite(options["penalty_exponent"]) or not .001 <= options["penalty_exponent"] <= 5:
        raise ValueError("penalty_exponent must be within the official [.001, 5] range.")
    if options["display"] not in (0, 1):
        raise ValueError("display must be 0 or 1; official surface plots are disabled.")
    if options["max_seconds"] is not None:
        if not np.isfinite(options["max_seconds"]) or options["max_seconds"] <= 0:
            raise ValueError("max_seconds must be finite and positive, or None.")
    return options


class _AdapterBudgetStop(Exception):
    """Internal sentinel only; caller exceptions are deliberately not caught."""


def run_opt(objective_constraints, x0, lower, upper, settings, constants=None):
    """Run the official derivative-free constrained search through an affine map.

    Callback: ``objective_constraints(x, constants) -> (scalar_f, residuals)``.
    Residuals use c<=0; numerical feasibility is max(c)<=feasibility_tol. The
    callback must itself reject invalid integrations through an explicit finite
    response; programming errors and caller budget exceptions escape unchanged.

    ``settings`` keys are DEFAULT_SETTINGS. Evaluation/time limits stop before
    the next callback, never by interrupting an ODE already running. Result is a
    dict; ``best_feasible`` is separate from the penalized solver-returned point.
    ``history`` uses physical design units but its f row is PENALIZED, not raw.
    ``evaluations`` is the genuine callback log and never includes invented data.
    """
    initial = np.array(x0, dtype=float, copy=True)
    lb = np.array(lower, dtype=float, copy=True)
    ub = np.array(upper, dtype=float, copy=True)
    if initial.ndim != 1 or initial.size < 2 or lb.shape != initial.shape or ub.shape != initial.shape:
        raise ValueError("x0/lower/upper must be equal-length vectors with at least two variables.")
    if not all(np.all(np.isfinite(a)) for a in (initial, lb, ub)):
        raise ValueError("Design vectors and finite bounds must be finite.")
    if np.any(ub <= lb) or np.any(initial < lb) or np.any(initial > ub):
        raise ValueError("Require lower < upper and x0 inside the bounds.")
    n = initial.size
    options = _validated_settings(settings, n)
    official = load_official()
    hyp = official["opt_hyp"]()
    hyp[:10] = [options["display"], options["tol_x"], options["tol_f"],
                options["feasibility_tol"], options["max_evals"], options["penalty"],
                options["penalty_exponent"], 1, .1, 0]
    hyp[18] = 0  # All supplied residuals are inequalities; no equality penalties.
    midpoint = (lb + ub) / 2
    halfwidth = (ub - lb) / 2
    scaled_initial = (initial - midpoint) / halfwidth
    log = []
    constraint_count = None
    best_feasible = None
    started = time.monotonic()
    stop_reason = None

    def score(f, g):
        # Exact inspected avg_cov_func expression, for reporting only.
        value = f + options["penalty"] * np.sum(
            g * (g > options["feasibility_tol"])
        ) ** options["penalty_exponent"]
        if not np.isfinite(value):
            raise FloatingPointError("Nonfinite augmented score; revise the finite failure scaling.")
        return float(value)

    def evaluate(scaled_x, forwarded_constants):
        nonlocal constraint_count, best_feasible, stop_reason
        if len(log) >= options["max_evals"]:
            stop_reason = "evaluation_budget"
            raise _AdapterBudgetStop()
        if options["max_seconds"] is not None and time.monotonic()-started >= options["max_seconds"]:
            stop_reason = "time_budget"
            raise _AdapterBudgetStop()
        x = midpoint + halfwidth * np.asarray(scaled_x, dtype=float)
        if np.any(x < lb) or np.any(x > ub):
            raise RuntimeError("Official native bounds produced a point outside the physical box.")
        f, g = objective_constraints(x.copy(), forwarded_constants)
        scalar = np.asarray(f, dtype=float)
        residuals = np.asarray(g, dtype=float)
        if scalar.ndim != 0 or not np.isfinite(scalar):
            raise ValueError("Callback objective must be one finite scalar.")
        if residuals.ndim != 1 or not residuals.size or not np.all(np.isfinite(residuals)):
            raise ValueError("Callback constraints must be a nonempty finite 1-D vector.")
        if constraint_count is None:
            constraint_count = residuals.size
        elif residuals.size != constraint_count:
            raise ValueError("Callback constraint vector size changed during optimization.")
        maximum = float(np.max(residuals))
        record = {
            "evaluation": len(log)+1, "x": x.copy(), "raw_objective": float(scalar),
            "constraints": residuals.copy(), "max_constraint": maximum,
            "max_violation": max(0., maximum), "penalized_objective": score(float(scalar), residuals),
            "feasible": maximum <= options["feasibility_tol"],
        }
        log.append(record)
        if record["feasible"] and (best_feasible is None or record["raw_objective"] < best_feasible["raw_objective"]):
            best_feasible = record
        return float(scalar), residuals.copy()

    try:
        # Nonzero native bounds avoid this revision's modification of zero bounds.
        result = official["nms"](evaluate, scaled_initial, -np.ones(n), np.ones(n), hyp, constants)
    except _AdapterBudgetStop:
        result = None
    elapsed = time.monotonic() - started
    history = np.empty((n+5, 0))
    iterations = None
    official_evaluations = None
    converged = False
    solver_record = None
    if result is not None:
        z, augmented, constraints, raw_history, official_evaluations, iterations = result
        x = midpoint + halfwidth * np.asarray(z)
        history = np.array(raw_history, copy=True)
        history[:n] = midpoint[:, None] + halfwidth[:, None] * history[:n]
        matches = [entry for entry in log if np.array_equal(entry["x"], x)]
        if not matches:
            raise RuntimeError("Official result does not correspond to a recorded callback.")
        solver_record = matches[-1]
        if not np.isclose(augmented, solver_record["penalized_objective"], rtol=1e-12, atol=1e-12):
            raise RuntimeError("Official returned score disagrees with the inspected penalty rule.")
        if not np.array_equal(np.asarray(constraints), solver_record["constraints"]):
            raise RuntimeError("Official returned constraints disagree with the callback record.")
        converged = bool(iterations > n*n and history[n+3, -1] < options["tol_x"]
                         and history[n+4, -1] < options["tol_f"])
        last_update = history[n+2, 0]
        previous_best = history[n, 0]
        for column in history.T[1:]:
            if column[n] < previous_best:
                last_update = column[n+2]
                previous_best = column[n]
        stalled = official_evaluations-last_update > .2*options["max_evals"]
        stop_reason = ("simplex_converged" if converged else "stalled" if stalled
                       else "evaluation_budget" if official_evaluations >= options["max_evals"]
                       else "returned_unclassified")
        if official_evaluations != len(log):
            raise RuntimeError("Official evaluation count disagrees with callback log.")
    selected = solver_record
    if selected is None and log:
        selected = best_feasible or min(log, key=lambda item: item["penalized_objective"])
    return {
        "x": None if selected is None else selected["x"].copy(),
        "raw_objective": None if selected is None else selected["raw_objective"],
        "penalized_objective": None if selected is None else selected["penalized_objective"],
        "constraints": None if selected is None else selected["constraints"].copy(),
        "feasible": False if selected is None else selected["feasible"],
        "best_feasible": best_feasible, "history": history, "evaluations": log,
        "n_evaluations": len(log), "n_iterations": None if iterations is None else int(iterations),
        "official_evaluations": None if official_evaluations is None else int(official_evaluations),
        "official_returned": result is not None, "converged": converged, "status": stop_reason,
        "elapsed_seconds": elapsed, "settings": options, "official_hyp": hyp.copy(),
        "requested_start": initial, "effective_start": midpoint+halfwidth*np.clip(scaled_initial, -.8, .8),
        "source": {k: official[k] for k in ("version", "revision", "nms_signature", "ode4u_signature")},
    }
