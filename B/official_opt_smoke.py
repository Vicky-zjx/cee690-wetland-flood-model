"""Run the actual official Example A and focused adapter checks, without ODEs.

Run from the project root: python3 -B -m B.official_opt_smoke
Only results/official_opt_smoke*.json is written. No figures are generated.
"""
from __future__ import annotations

import contextlib
import hashlib
import importlib.util
import inspect
import io
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

from B.official_opt import OFFICIAL_CHECKOUT, ROOT, load_official, run_opt


def json_value(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def snapshot():
    tracked = subprocess.check_output(
        ["git", "-C", str(OFFICIAL_CHECKOUT), "ls-files", "-z"]
    ).decode().strip("\0").split("\0")
    vendor = {name: hashlib.sha256((OFFICIAL_CHECKOUT/name).read_bytes()).hexdigest()
              for name in tracked if (OFFICIAL_CHECKOUT/name).is_file()}
    a_root = ROOT.parent/"A"
    a2_paths = [a_root/"wetland_flood_model.py", a_root/"algebraic_plots.py",
                a_root/"results/baseline.npz", a_root/"results/baseline_summary.json",
                *sorted((a_root/"plots").glob("[A-J]_*.png"))]
    a2 = {str(p.relative_to(ROOT.parent)): hashlib.sha256(p.read_bytes()).hexdigest()
          for p in a2_paths if p.is_file()}
    return {"vendor_tracked_sha256": vendor, "a2_sha256": a2}


def main():
    sys.dont_write_bytecode = True
    before = snapshot()
    vendor_status_before = subprocess.check_output(
        ["git", "-C", str(OFFICIAL_CHECKOUT), "status", "--short"], text=True
    )
    official = load_official()
    from multivarious.utl import StableNamespace
    import matplotlib
    import scipy
    example_source = OFFICIAL_CHECKOUT/"examples/opt_example.py"
    spec = importlib.util.spec_from_file_location("official_opt_example", example_source)
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)
    constants = StableNamespace(a=[-.4,.2,.5,1.4,1.4],
                                b=[1.,-.5,.5,-1.4,-1.4],c=[0.,.8,.2])
    seed = 20261001
    np.random.seed(seed)
    x0 = np.array([.8,.8])+.10*np.random.rand(2)
    direct_settings = [0,.02,.02,.001,200,.7,.5,1,.05]
    direct_log = []

    def original_analysis(x, consts):
        if len(direct_log) >= 220:
            raise RuntimeError("Official smoke exceeded its emergency 220-call bound.")
        f, g = example.analysis(x, consts)
        direct_log.append({"x":x.copy(), "raw_objective":float(f), "constraints":g.copy()})
        return f, g

    started = time.monotonic()
    stdout = io.StringIO()
    with contextlib.redirect_stdout(stdout):
        direct = official["nms"](original_analysis, x0, np.zeros(2), np.ones(2),
                                  direct_settings, constants)
    x, augmented, g, history, evaluations, iterations = direct
    raw, check_g = example.analysis(x, constants)
    assert np.array_equal(g, check_g)
    assert len(direct_log) == evaluations
    direct_report = {
        "purpose":"Actual official Example A analysis with official NMS; not wetland optimization.",
        "source":{k:official[k] for k in ("version","revision","nms_signature","ode4u_signature")},
        "example_source":"A/work/multivarious/examples/opt_example.py:analysis; Example A constants",
        "example_source_sha256":hashlib.sha256(example_source.read_bytes()).hexdigest(),
        "numpy":np.__version__, "scipy":scipy.__version__, "matplotlib":matplotlib.__version__,
        "seed":seed, "settings":direct_settings, "x0":x0, "x":x,
        "raw_objective":raw, "penalized_objective":augmented, "constraints":g,
        "n_evaluations":evaluations, "n_iterations":iterations,
        "feasible_at_return":bool(np.max(g)<=.001), "constraint_tolerance":.001,
        "elapsed_seconds":time.monotonic()-started, "history":history,
        "evaluation_log":direct_log, "solver_stdout":stdout.getvalue(),
        "note":"Direct official smoke retains its zero-bound adjustment; adapter tests exact bounds separately.",
    }
    options = {"max_evals":200, "tol_x":.02,"tol_f":.02,"feasibility_tol":.001,
               "penalty":.7,"penalty_exponent":.5,"max_seconds":10.}
    with contextlib.redirect_stdout(stdout):
        adapted = run_opt(example.analysis, x0, np.zeros(2), np.ones(2), options, constants)
    assert adapted["n_evaluations"] <= 200
    assert adapted["best_feasible"] is not None
    assert adapted["raw_objective"] < .03
    assert all(np.all((r["x"] >= 0)&(r["x"] <= 1)) for r in adapted["evaluations"])

    # Exact zero remains reachable, and native mutation never reaches caller bounds.
    lower = np.zeros(2)
    upper = np.ones(2)
    def boundary(x, consts):
        return float((x[0]+.2)**2+(x[1]-.4)**2), np.array([-1.])
    with contextlib.redirect_stdout(stdout):
        boundary_result = run_opt(boundary, [.6,.6], lower, upper,
                                  {"max_evals":150,"max_seconds":10.})
    assert np.array_equal(lower, np.zeros(2)) and np.array_equal(upper,np.ones(2))
    assert any(entry["x"][0] == 0. for entry in boundary_result["evaluations"])
    assert all(np.all((r["x"] >= 0)&(r["x"] <= 1)) for r in boundary_result["evaluations"])

    # The official soft evaluation limit must never exceed the adapter hard cap.
    with contextlib.redirect_stdout(stdout):
        budget = run_opt(boundary, [.6,.6], [0.,0.], [1.,1.], {"max_evals":5})
    assert budget["status"] == "evaluation_budget" and budget["n_evaluations"] == 5

    # Only the adapter's sentinel is handled; caller budget/errors escape intact.
    sentinel = RuntimeError("Caller budget sentinel: do not swallow")
    calls = 0
    def throws(x, consts):
        nonlocal calls
        calls += 1
        if calls == 3:
            raise sentinel
        return boundary(x,consts)
    try:
        run_opt(throws,[.6,.6],[0.,0.],[1.,1.],{"max_evals":50})
    except RuntimeError as error:
        assert error is sentinel
    else:
        raise AssertionError("Caller exception was swallowed")
    calls = 0
    def changes_size(x, consts):
        nonlocal calls
        calls += 1
        return 1., np.array([-1.] if calls == 1 else [-1.,-2.])
    try:
        run_opt(changes_size,[.6,.6],[0.,0.],[1.,1.],{"max_evals":50})
    except ValueError as error:
        assert "size changed" in str(error)
    else:
        raise AssertionError("Changing constraint dimension was accepted")
    try:
        run_opt(lambda x,c:(float('nan'),np.array([-1.])), [.6,.6],[0.,0.],[1.,1.],{"max_evals":50})
    except ValueError as error:
        assert "finite scalar" in str(error)
    else:
        raise AssertionError("Nonfinite objective was accepted")

    # A finite failure response is recorded as infeasible, never as a successful solve.
    with contextlib.redirect_stdout(stdout):
        failed = run_opt(lambda x,c:(10.,np.array([1.,1.])),[.6,.6],[0.,0.],[1.,1.],{"max_evals":12})
    assert failed["best_feasible"] is None and not failed["feasible"]
    after = snapshot()
    vendor_status_after = subprocess.check_output(
        ["git", "-C", str(OFFICIAL_CHECKOUT), "status", "--short"], text=True
    )
    assert before == after
    assert vendor_status_before == vendor_status_after
    report = {
        "scope":"Official example and adapter plumbing only; no hydrologic optimization or model execution.",
        "run_command":"python3 -B -m B.official_opt_smoke",
        "checks":{
            "official_example_A":True,"adapter_example_A":True,"exact_zero_bound_reachable":True,
            "caller_bound_arrays_unchanged":True,"all_callbacks_in_bounds":True,
            "hard_evaluation_budget":True,"caller_budget_exception_escapes":True,
            "fixed_constraint_count_enforced":True,"nonfinite_callback_rejected":True,
            "finite_failure_not_feasible":True,"vendor_and_a2_sha256_unchanged":True,
            "vendor_git_status_unchanged":True,
        },
        "adapter_example":adapted,"boundary_test":boundary_result,
        "hard_budget_test":budget,"finite_failure_test":failed,
        "vendor_git_status":vendor_status_after,"protected_before_after":before,
    }
    out = ROOT/"results"
    out.mkdir(parents=True,exist_ok=True)
    for name,value in (("official_opt_smoke_direct.json",direct_report),
                       ("official_opt_smoke_adapter.json",report)):
        (out/name).write_text(json.dumps(value,default=json_value,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"checks":report["checks"],
                      "direct_evaluations":int(evaluations),
                      "adapter_evaluations":adapted["n_evaluations"],
                      "adapter_status":adapted["status"],
                      "adapter_raw_objective":adapted["raw_objective"],
                      "adapter_max_constraint":float(np.max(adapted["constraints"]))},indent=2))


if __name__ == "__main__":
    main()
