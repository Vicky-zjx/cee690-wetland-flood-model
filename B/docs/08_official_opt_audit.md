# Official OPT and ODE audit for HW B.2

Audit date: 2026-10-01. This is an audit of the existing local dependency and a small optimization example, not a completed wetland optimization.

## Source and execution status

The inspected external checkout is `A/work/multivarious`, origin `https://github.com/hpgavin/multivarious.git`, revision **`06d7aed5b541ec7fe26a6dc69be6aa7fb610aa02`**, package version **1.2.1**. No fetch, pull, update, vendor edit, or replacement solver was used. Its tracked files were clean; existing untracked `.DS_Store` and Python bytecode directories were present. The smoke run left both tracked-file hashes and the complete vendor Git-status text unchanged. A.2 scripts, A–J figures, and existing baseline result hashes were also unchanged.

No separate professor OPT Python example was identified among the top-level Downloads filenames. The actual official checkout supplies [`examples/opt_example.py`](https://github.com/hpgavin/multivarious/blob/06d7aed5b541ec7fe26a6dc69be6aa7fb610aa02/examples/opt_example.py). It explicitly presents `ors`, `nms`, and `sqp` as solver choices; its default main block selects `ors`. This project selects the official derivative-free **`nms`** implementation because the peak-stage objective and grouped path extrema are nonsmooth. This is a documented implementation choice, not a claim that the assignment specifically mandates Nelder–Mead. No gradients, smoothing of extrema, or alternate optimizer are introduced. Nelder–Mead remains a local heuristic without a global optimality guarantee.

Executed:

```sh
python3 -B -m B.official_opt_smoke
```

The smoke imports the actual Example A `analysis` function and its published constants. It uses seed `20261001`, the example's random perturbation of its initial guess, `nms`, no plots, a 200-evaluation nominal budget, and an emergency 220-call cap for the direct vendor call. Example noise amplitude remains zero. It does not execute the example's interactive main block or copy its implementation into this project.

| Executed case | Callback evaluations | Iterations | Raw objective | Maximum residual | Constraint tolerance |
|---|---:|---:|---:|---:|---:|
| Direct official Example A with `nms` | 66 | 21 | 0.018410663010427982 | 0.0009287612994149513 | 0.001 |
| Same Example A through the project adapter | 66 | 21 | 0.018406294158279705 | 0.0009541637898663469 | 0.001 |

These points satisfy the stated numerical tolerance, not exact nonpositive residuals. The small numerical difference comes from the adapter preserving the exact physical zero lower bound instead of the vendor's zero-bound adjustment. Both returns meet the inspected simplex convergence criterion for this two-variable test. Genuine histories and callback logs are saved locally in `B/results/official_opt_smoke_direct.json` and `B/results/official_opt_smoke_adapter.json`; they are not wetland results. Runtime versions were NumPy 2.2.6, SciPy 1.16.2, and Matplotlib 3.10.6. The official example additionally imports `rich`; that presentation dependency is not used by the adapter itself.

## Actual interfaces

Inspected source: [`opt/nms.py`](https://github.com/hpgavin/multivarious/blob/06d7aed5b541ec7fe26a6dc69be6aa7fb610aa02/multivarious/opt/nms.py), [`utl/opt_hyp.py`](https://github.com/hpgavin/multivarious/blob/06d7aed5b541ec7fe26a6dc69be6aa7fb610aa02/multivarious/utl/opt_hyp.py), [`utl/avg_cov_func.py`](https://github.com/hpgavin/multivarious/blob/06d7aed5b541ec7fe26a6dc69be6aa7fb610aa02/multivarious/utl/avg_cov_func.py), and [`ode/ode4u.py`](https://github.com/hpgavin/multivarious/blob/06d7aed5b541ec7fe26a6dc69be6aa7fb610aa02/multivarious/ode/ode4u.py).

```python
nms(func, v_init, v_lb=None, v_ub=None, hyp_in=None, consts=1.0)
f, residuals = func(v, consts)
v, penalized_f, residuals, history, function_evals, iteration = nms(...)

ode4u(odefun, time, x0, u=None, c=None)
time, states, derivatives, outputs = ode4u(...)
```

The last two NMS returns are **evaluation count, then iteration count** in executable source; its docstring lists those two in the opposite order. The adapter follows executable source and checks counts against the callback log. The ODE solver retains the original project interface; no hydrologic simulation was needed for this adapter smoke.

`func` supplies one finite scalar and a nonempty fixed-length residual vector. Mathematical inequalities are `residual <= 0`; `avg_cov_func` activates penalties only for `residual > tol_g`. `nms` also has a strict `< tol_g` test in its optional stop-when-feasible mode, which this adapter disables. The reporting criterion is explicitly `max(residuals) <= feasibility_tol`. Safety/applicability margins and all constraint scaling remain the responsibility of the B.2 shooting layer.

## Bounds and starting point

The inspected NMS changes exact-zero lower bounds to `1e-4` and exact-zero upper bounds to `-1e-4`, and can mutate caller-provided bound arrays. That is inappropriate for exact gate closure and zero vegetation cover.

`run_opt` maps each physical variable bijectively to `z = (x-midpoint)/halfwidth`, passes **nonzero native bounds `[-1,1]`**, and maps callback points back to the original units. The original feasible interval and categorical enumeration are unchanged. Native `avg_cov_func` enforces the box by projecting proposed scaled points into `[-1,1]`. The adapter does not clip reconstructed gate histories, hydrologic states, or constraint residuals. Physical bound arrays are copied, and every actual callback point is independently checked against them. A focused smoke verifies that a physical zero bound is actually reached.

NMS moves its initial scaled point into `[-0.8,0.8]` before constructing the simplex. The result records both `requested_start` and `effective_start`; a reference run must be simulated independently at its exact specified decision vector. This does not approximate or change the chosen gate-knot time grid. This adapter supports at least two continuous variables with strictly separated finite bounds; the B.2 full gate/covers vector satisfies that requirement. Fixed parameters are excluded from the optimization vector.

## Settings and penalty semantics

Public adapter:

```python
run_opt(objective_constraints, x0, lower, upper, settings, constants=None)
```

The callback receives `(physical_x, constants)` and returns a conditioned **raw** objective plus already dimensionless residuals. The adapter adds no hydrologic penalty or alternative objective. Any fixed positive objective scaling/additive shift is configured in the shooting layer and must be inverted when reporting peak stage in metres.

| Adapter setting | Default | Official entry | Meaning |
|---|---:|---:|---|
| `display` | 0 | 0 | 0/1 only; surface plotting disabled |
| `tol_x` | 0.001 | 1 | Relative simplex convergence metric, not a physical gate tolerance |
| `tol_f` | 0.001 | 2 | Relative penalized-objective convergence metric |
| `feasibility_tol` | 0.000001 | 3 | Dimensionless residual tolerance |
| `max_evals` | 200 | 4 | Official nominal limit plus adapter hard callback cap |
| `penalty` | 10 | 5 | Official dimensionless penalty coefficient |
| `penalty_exponent` | 1 | 6 | Official exponent, allowed range 0.001–5 |
| — | 1 | 7 | One deterministic callback per evaluation; no sampling/risk adjustment |
| — | 0.1 | 8 | Inactive averaging COV target when count is one |
| — | 0 | 9 | Stop-when-feasible disabled |
| `max_seconds` | null | None | Optional wall budget checked before each callback |

The complete inspected 19-entry `opt_hyp` default vector is:

```text
[1, 1e-3, 1e-3, 0, 1000, 10, 1, 1, .1, 0, 0, 1, 25, 35, 1e-6, 2, 1e-6, .1, 0]
```

Entries 10–13 configure surface plotting; 14 is a finite-difference minimum; 15 is a penalty-type setting; 16–17 are finite-difference change limits; 18 counts equality constraints. These entries are unused by this NMS path except that equality count is explicitly kept zero. In particular, `avg_cov_func` does not branch on entry 15. The actual `official_hyp` passed is returned with each run. `opt_hyp` takes absolute values of input entries and clips the exponent; the adapter validates settings before that normalization so invalid negative settings are not silently accepted.

With one callback per evaluation, the exact official augmented score is:

```text
F = f + penalty * (sum(residual[j] for residual[j] > tol_g)) ** exponent
```

It sums the entire positive residual exceeding tolerance, **not** only the excess over tolerance. This gives a small discontinuity at positive `tol_g`. The penalty is applied once by the official module. The adapter's duplicate expression is used solely to label logs and verify returned scores. All residual scales must be fixed and dimensionally consistent before calling this routine; the penalty coefficient is numerical, not a physical price or a new weighted objective. Low `F`, convergence, or a returned point does not prove feasibility. Independently verify all constraints and preserve the original objective in metres.

## History, termination, and budgets

Official `history` has shape `(n+5, iterations+1)` with rows `[design variables; penalized F; max residual; cumulative callback count; relative variable metric; relative objective metric]`. [`plot_cvg_hst`](https://github.com/hpgavin/multivarious/blob/06d7aed5b541ec7fe26a6dc69be6aa7fb610aa02/multivarious/utl/plot_cvg_hst.py) expects exactly this layout. The adapter converts only design rows back to physical units; the objective row remains explicitly penalized. Its separate `evaluations` log records every actual raw objective and full residual vector, suitable for the requested evaluation-history figure. No interpolation or invented convergence trace is supplied.

The routine checks simplex convergence only after `iteration > n*n`, using both relative metrics. For the full 24 gate values plus two covers, that means more than 676 iterations before its formal convergence test can pass. A modest evaluation budget may therefore find feasible improvements without ever declaring convergence. NMS also stops when no strict improvement in its best penalized score has occurred for more than `0.2*max_evals` evaluations, or at its evaluation limit. It has no maximum wall-time argument and returns no explicit status enum. The adapter reports source-derived `simplex_converged`, `stalled`, `evaluation_budget`, or `returned_unclassified`; none is a physical-validity or global-optimality claim.

Official NMS can exceed its nominal evaluation budget within a simplex operation. The adapter checks its own limit before every callback, catches only its private budget sentinel, and returns `official_returned=False` on that stop. No history is fabricated when the vendor routine did not return; completed callback records remain available. Optional wall time is also checked before callbacks, so a single in-flight ODE may finish after the nominal deadline. There is no unsafe mid-integration interrupt.

`best_feasible` is the best raw-objective point among actual calls satisfying the stated tolerance. On a normal return, top-level `x` describes the solver-returned point, even if it differs from `best_feasible`. On an adapter budget stop, `x` is the recorded best feasible point, or the best augmented-score point if none is feasible; `official_returned` distinguishes these cases. With no successful callback, `x` is null. The driver must re-integrate the chosen candidate afresh and on a finer grid before reporting a solution.

A caller-defined budget exception escapes unchanged. Unexpected programming/integration exceptions are not broadly swallowed. The shooting layer may catch only its explicitly expected invalid-trajectory exceptions and return its documented finite failure objective and fixed-size violation vector. NaN/Inf objectives/residuals or a changed vector length are rejected by the adapter. A finite failure response was tested to remain infeasible, not to masquerade as a solution.

## Checks and remaining limits

Executed checks cover the actual official example, its adapter run, exact-zero bound reachability, unchanged caller bound arrays, bounded callback points, strict evaluation cap, propagation of caller budget exceptions, fixed residual count, rejection of nonfinite objectives, and an always-failing finite response. Vendor and A.2 hashes/status were compared before and after.

This establishes functioning official OPT plumbing only. No supported wetland safety levels, local gate dimensions/coefficient evidence, ecological mapping, local forcing configuration, or final constrained hydrologic result is supplied by this example. Those remain separate preflight requirements. No substantive local-data or illustrative wetland optimization should run until its required parameters and approvals exist.
