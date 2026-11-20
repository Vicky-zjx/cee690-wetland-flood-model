# Professor example and official uncertainty-tool audit

## Scope and verified implementation

The supplied [fair_lite_demo.py](/Users/zhangjiaxuan/Downloads/fair_lite_demo.py) was inspected as source, not imported or executed. Its SHA-256 is `e5f9d52a926bdfe28075d3a3ecd95abcc239a675696540416605f2a2f9ab9ccd`. The example's climate equations, sample count, uncertainty percentages, objective and constraints are not instructions or evidence for this wetland model.

The existing official checkout is [A/work/multivarious](../../../A/work/multivarious), revision **06d7aed5b541ec7fe26a6dc69be6aa7fb610aa02**, package metadata version **1.2.1**. Its tracked files were unchanged at inspection. `multivarious` is not installed as a discoverable package in the checked Miniconda environment; existing workflows use the explicit checkout path. Do not replace this pinned dependency or change A/B/C to match the example.

Nonphysical API smoke checks imported that checkout with bytecode disabled and verified seeded toy draws, singleton shapes, first-order term construction and the single-output scaling behavior described below. They did not fit a model or execute a physical simulation. Available dependencies were NumPy 2.2.6, SciPy 1.16.2 and Matplotlib 3.10.6. These cover the inspected sampling and plotting functions. The package `pyproject.toml` does not declare its runtime dependencies. The professor script additionally imports `rich` and climate-specific local modules; those are not dependencies of a narrowly adapted wetland uncertainty runner.

## What the professor's workflow actually does

```python
analysis(v, c) -> (f, g)
mcs(v_medn, v_covn, N_mcs) -> (dataX, dataY)
```

These two functions are defined in the example, not exported package algorithms. `analysis` mutates `c._v`, constructs emissions, randomly redraws seven climate constants inside every call, runs the full official `ode4u` model, extracts an objective and three constraint residuals, and appends them to an already-open history stream. `mcs` uses global `c` and `pdf_plots`, samples a lognormal design vector for each iteration, calls `analysis`, reloads its CSV and sends transposed arrays to the official scatter routine. Its returned CSV-oriented arrays have observations in rows; official fitting/plotting inputs have observations in columns.

The reusable lesson is **sample declared uncertain inputs, run the complete physical analysis, retain physical outputs and constraints, then summarize their distribution**. For the wetland study, the three selected categories and their numerical specification status are in [parameter_mapping.md](parameter_mapping.md). Neither 10,000 samples nor 30% CV from the climate demonstration is adopted. Those are example settings, not uncertainty evidence for infiltration, routing or outlet resistance.

The example is not directly executable against this checkout without adaptation:

- `load_model_npz` is absent; no equivalent loader appears in the checked fitting module.
- `compute_model` exists inside `multivarious.fit.mimo_srs`, but is not exported by `multivarious.fit`.
- `opt_open_history` is absent from the checked package.
- The example unpacks the last three `mimo_srs` returns in a different order from the actual implementation.

These are version/API differences, not a reason to overwrite the official source or reconstruct missing functions from memory.

## Official random-variable APIs

```python
normal.rnd(meanX=0.0, sdvnX=1.0, N=1, R=None, seed=None)
lognormal.rnd(mednX, covnX, N, R=None, seed=None)
uniform.rnd(a, b, N, R=None, seed=None)
correlated_rvs(R, n, N=1, seed=None) -> (R_used, normal_values, uniform_values)
```

Sources: [normal.py](../../../A/work/multivarious/multivarious/rvs/normal.py), [lognormal.py](../../../A/work/multivarious/multivarious/rvs/lognormal.py), [uniform.py](../../../A/work/multivarious/multivarious/rvs/uniform.py) and [correlated_rvs.py](../../../A/work/multivarious/multivarious/utl/correlated_rvs.py).

| Function | Actual parameter convention | Important limitation |
| --- | --- | --- |
| Normal | Arithmetic mean and standard deviation | Standard deviation must be positive. Positive means do not prevent physically negative draws. |
| Lognormal | Median and coefficient of variation SD/mean | Both must be positive; this is not an arithmetic-mean/SD interface. |
| Uniform | Lower and upper bounds | Each upper bound must exceed its lower bound. A bounded support is still an assumption requiring justification. |

For a lognormal input with median m and CV c, the implementation uses `log(X)` mean `log(m)` and variance `log(1+c*c)`. Its arithmetic mean is `m*sqrt(1+c*c)`. To specify an approved arithmetic mean μ, the input median must be `μ/sqrt(1+c*c)`. Fixed parameters belong outside the sampled vector: CV = 0, SD = 0 and equal uniform bounds are rejected by these APIs. Parameter arrays must have matching lengths; a scalar CV is not automatically broadcast across three medians by validation.

Before squeezing, samples have shape **(n_parameters, N_samples)**. Actual shapes verified were `(3,5)`, `(3,)` for one three-parameter draw, `(5,)` for five scalar draws, and a zero-dimensional array for one scalar draw. Normalize deliberately with `np.asarray(draw).reshape(3, N)` for a known three-variable joint call; never infer the sample axis from a squeezed array.

`R` acts on the latent Gaussian variables. It is the physical Pearson correlation for normal marginals, but generally not for transformed lognormal or uniform marginals. The underlying routine can change the supplied matrix's diagonal in place and shrink a problematic matrix toward identity. Its loose eigenvalue tolerance does not guarantee a safe square root for a remaining slightly negative eigenvalue. Validate dimensions, finiteness, symmetry, diagonal and positive definiteness before a correlated call; pass a copy and record the intended/effective dependence. Do not silently allow a numerical correlation repair to redefine the approved uncertainty law. Independence through `R=None` is an explicit study assumption, not a fact established by absent data.

### Reproducible call patterns

`rnd` constructs `np.random.default_rng(seed)` internally. Therefore, a global `np.random.seed(...)` does not seed these functions. Repeating the same integer seed inside a loop repeats the same draw. The following are interface recipes, not approved distributions or sample counts:

```python
# If all three approved marginals use the same family:
X = np.asarray(lognormal.rnd(approved_medians, approved_CVs,
                            approved_N, R=approved_R, seed=master_seed))
X = X.reshape(3, approved_N)
for sample_id in range(approved_N):
    v = X[:, sample_id].copy()
    # Both physical cases receive this exact v; neither redraws it.

# For independently sampled different marginal families, if approved:
children = np.random.SeedSequence(master_seed).spawn(3)
# Pass a distinct saved child SeedSequence as seed to each marginal call.
# Using the same seed for three separate scalar calls creates dependence.
```

Archive the full sampled parameter table, variable order/units, laws, seed or child-sequence provenance, official revision and source hashes. Sampling occurs once per realization, never inside an ODE derivative or separately for the two paired cases.

For the newly approved log-uniform multipliers, call the official `uniform.rnd` on log(0.5) to log(2) for each coordinate and exponentiate once; `lognormal.rnd` is a different distribution. The log-uniform multiplier has median 1 but arithmetic mean 1.5/log(4), approximately 1.08202. Thus symmetry around the baseline is logarithmic, not arithmetic. Record both multipliers and the resulting physical parameters.

## Official empirical-distribution and scatter APIs

```python
plot_ECDF_ci(data, confidence_level, fig_num,
             x_label='sorted sample values', save_plots=True)
# returns sorted_x, Fx, lower_ci, upper_ci

plot_scatter_dist(dataX, dataY=None, fig_num=100, var_names=None,
                  n_bins=20, font_size=15, ci=0.95)
# returns a matplotlib Figure
```

Sources: [plot_ECDF_ci.py](../../../A/work/multivarious/multivarious/utl/plot_ECDF_ci.py) and [plot_scatter_dist.py](../../../A/work/multivarious/multivarious/utl/plot_scatter_dist.py).

The confidence arguments use **different units**: ECDF `95` means 95%; scatter `0.95` means 95%. The ECDF routine flattens and sorts the sample and uses **rank/(N+1)** plotting positions. It does not implement the ties-aggregated, right-continuous `count(x_i <= x)/N` ECDF used in the earlier weather pilot. Its active limits are `Beta(i, N+1-i)` quantiles for each order statistic. The commented DKW and Dirichlet sections are inactive. These are **pointwise order-statistic intervals**, not simultaneous coverage of the entire curve; their exact continuous-i.i.d. interpretation does not extend automatically to ties or atoms. The official graph must be labeled with these semantics rather than silently relabeled as the earlier ECDF. Bands concern finite-sample uncertainty conditional on the assumed sampling law, not model calibration or empirical certainty about the proposed parameter distributions.

The routine also annotates a median using `x[int(round(N/2))]`, which differs from the conventional sample median for even and some odd N; report an independently computed conventional median in tables. The new wrapper corrects that displayed annotation after the official call, preserves the original index-based value in metadata and does not alter the official source or returned plotting positions/bands. Its CV is unstable/undefined when a signed output's mean is near zero. It has no finite-data or minimum-size guard. Precheck inputs; do not replace a negative paired effect, undefined physical result or missing trajectory with zero to make it plot. Unchanged official plotting positions/bands and any limitations must remain identifiable.

Scatter inputs are `(n_variables, N_samples)`. Optional labels are `{'X': [...], 'Y': [...]}`. The upper triangle displays Pearson correlations and approximate Fisher-z intervals, using `1/sqrt(N-3)`; N must exceed three, and constant or nonfinite variables are unsuitable. It does not compute Sobol indices, causal effects, conditional regression coefficients or prediction intervals. Nonlinear/skewed outcomes and dependent input coordinates limit the interpretation. With `n_bins <= 0`, its diagonal also uses rank/(N+1), not the earlier standard ECDF.

Both routines change Matplotlib state and enable interactive mode. A separate batch plotting process can use `Agg`, an `rc_context` with `text.usetex=False`, explicit output paths and figure closure. The inspected `format_plot` does **not** turn on LaTeX; no monkeypatch is needed. Set `save_plots=False` in `plot_ECDF_ci` to avoid its default PDF in the current directory, then save the actual returned figure handle obtained through its known figure number in the new study folder.

## Official FIT / mimo_srs: review only at this stage

```python
mimo_srs(dataX, dataY, max_order=2, pTrain=70, scaling=1,
         L1_pnlty=1.0, basis_fctn='H', var_names=None)
# Actual return order:
# ordr, coeff, meanX, meanY, invTX, TY, testX, testY, testModelY

# Available inside the module, not the package-level fit namespace:
compute_model(ordr, coeff, meanX, meanY, invTX, TY, dataX,
              scaling, basis_fctn='H') -> (modelY, B)
```

Source: [mimo_srs.py](../../../A/work/multivarious/multivarious/fit/mimo_srs.py). The actual function signature takes precedence over inconsistent docstring defaults.

- Use a two-element mutable scaling specification such as `[1,1]`; the scalar default is immediately indexed and is not usable as written.
- The top-level routine scales **all supplied X and Y before splitting**, then removes columns using a Chauvenet rule applied jointly to X and Y. This can leak holdout information and remove scientifically genuine extremes. It must not be used silently to redefine an uncertainty ensemble.
- The internal split uses `np.random.permutation`, a different RNG interface from `rvs.rnd`; it has no `seed=` argument. Seed and record a separate split, keeping paired sample IDs together. A train-only helper workflow can avoid preprocessing leakage, but is a separate explicitly documented adaptation, not an unchanged call to the top-level function.
- For one variable, scaling option 1 removes its mean but leaves the transform equal to 1; it does not divide by its standard deviation. This was verified with a nonphysical four-value input. Constant variables can also create singular scaling/correlation calculations.
- Three inputs at `max_order=1` produce four terms. With `'H'`, the nonconstant first-order basis is the implementation's attenuated Hermite function, not an unrestricted affine power regression. Preserve the selected basis in any saved model.
- `L1_pnlty` is a numerical fitting choice; 100 from the demonstration is not automatically appropriate. The positive-penalty path uses the official L1 routine; the zero-penalty helper branch passes `Zy` to least squares without the transpose required by its documented shape, so it must not be assumed to work as a drop-in unregularized fit.
- A surrogate cannot establish physical feasibility or replace full paired trajectories. Evaluate a frozen surrogate on independent held-out physical results, report errors near limits, and keep rejected/unsafe cases in the full analysis accounting. No new FIT run is part of this parameter-uncertainty task. Sampling, paired full-ODE evaluation, empirical distributions and scatter diagnostics are the minimum relevant adaptation of the demonstration. A new surrogate is unnecessary for those primary results; a separately requested later fit would need completed physical training and held-out results.

## Narrow adaptation to the proposed three-parameter paired study

For an approved parameter specification, `analysis(v,c)` should be deterministic for a given parameter vector, saved storm and model variant. Apply infiltration capacity, the approved routing coordinate and the approved drainage coordinate to a fresh per-case parameter object. Keep geometry, initial storages, boundaries, wetland design, safety limits and feedback policy unchanged. Run both full models with the exact same realization and preserve their actual gate histories. If outlet resistance varies, compute initial-discharge-relative metrics from that realization's own `Q_out(0)`; do not reuse the fixed baseline reference.

Store one row per parameter realization containing all three inputs, both cases' physical metrics, signed paired differences, hydraulic/numerical completion flags and safety residuals. A valid unsafe trajectory or negative wetland effect is a retained result. An unsupported or numerically invalid trajectory has a recorded cause and unavailable physical outputs, not an invented penalized peak. Never silently resample failures, clip parameter draws or filter extreme physical outcomes to improve the distribution. Any conditional valid-output plots must expose their denominator and exclusions.

The user subsequently approved the existing `K_u`, a common multiplier of the three runoff-routing rates, and `R_r` as the coordinates, with exploratory multipliers 0.5, 1 and 2. Their 27-combination factorial grid is a deterministic sensitivity design, not an iid Monte Carlo sample. Descriptive fractions across those 27 prescribed points can be labeled design fractions; they do not justify the official Beta or Fisher sampling intervals or a parameter-event probability. A subsequent user approval also permits 64 paired iid exploratory log-uniform multiplier realizations on [0.5, 2]; this is a separate random-sampling study, not a reinterpretation of the 27-point grid. Its distributions are explicit exploratory assumptions, not inferred physical uncertainty. Use the recorded approved joint sampling configuration and retain the distinction between these two designs.

Completed application: the immutable exploratory 64-pair parameter ensemble has now been simulated with the full physical models at 15 s. All 64 paired rows are valid and classified safety-feasible under the existing upstream limits; no draw was filtered or resampled. The runner reports a passed preservation check for 2,659 pre-existing files. The official plotting wrapper verified the configuration hash, planned sample IDs, parameter mappings and signed differences, then created six physical-outcome distribution figures and the five-variable parameter/effect scatter figure under `plots/official_mc/`. All six output columns have 64 distinct values in this realization. The only plotting warnings were seven expected noninteractive Agg display notices; saved files were generated successfully and the paired-stage distribution and scatter layout were visually inspected. Mean paired peak-stage reduction is 0.0915734832 m (range 0.0342844646–0.1592176541 m); mean paired peak-discharge reduction is 0.1765063613 m³/s. All 64 reductions are positive, which is a result for these exploratory draws, not proof of universal benefit or calibrated flood-risk reduction. No surrogate fit was performed. The earlier 16 weather pairs were not used as parameter draws. The wrapper continues to require complete configured-draw provenance for any repeat analysis.

The adapter was API-smoke-tested on eight explicitly labeled nonphysical toy rows saved only under `/private/tmp`; negative and constant-zero outputs were retained, incorrect study provenance was rejected, and official plotting source files were unchanged. This checks the plotting interface, not the wetland physical model or any parameter-uncertainty result.

Within each realization, the linear outlet law uses the same resistance and boundary for both cases, so the paired peaks satisfy `delta_H_peak = R_r * delta_Q_peak` (the saved 64-row check has roundoff-level residual). Consequently, the high correlation between `R_r` and stage reduction is partly inherent in the stage–discharge mapping. It must not be interpreted as an independent causal sensitivity estimate; the discharge reduction and deterministic controlled comparisons provide different context.
