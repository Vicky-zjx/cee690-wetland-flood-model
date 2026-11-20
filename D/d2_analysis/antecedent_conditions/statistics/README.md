# Independent importance-weight audit

This is a saved-result statistical supplement to the antecedent-condition study. It audits the **previous five-event, 64-vector ensemble**, not the new antecedent trajectories. No ODE, new physical parameter vector, random draw, fit or model change is performed. All original results and parameter values remain unchanged.

The proposal density is q(m)=∏[1/(m_j ln 4)] on [0.5,2]³. The exploratory alternative is independent uniform, p(m)=1/(1.5³), so p/q=∏[m_j ln(4)/1.5]. Normalizing these density ratios gives the same per-draw weights for every event. The three categories remain infiltration, common routing, and river resistance; no routing rate is split into a new independent parameter.

The independently recomputed weight ESS is **43.99349928/64** (68.74%). The largest normalized weight is **4.9785%**, for `mc_00018`. The largest 1, 5 and 10 weights contain respectively **4.9785%, 20.8863%, and 35.6483%** of the weight. Equal weighting would give 1.5625%, 7.8125%, and 15.625%. The smallest weight is 0.2302%. The mean unnormalized p/q is 0.964298, rather than exactly one in this finite sample; self-normalization divides by the realized sum.

Two independent ESS formulas agree. The largest difference from the saved normalized weights is 1.39e-17; all previously saved reweighted means and zero masses agree within the tolerances recorded in `prior_result_verification.csv`. All 320 original pairs remain valid and upstream-safety-feasible; no zero or unfavorable metric is excluded.

## Weight concentration versus outcome concentration

For an outcome y_i, its signed contribution to the self-normalized mean is c_i=w_i y_i. The tables distinguish samples ranked by **weight** from samples ranked by **absolute outcome contribution**. The latter share is |c_i|/Σ|c_j|, so it remains interpretable for duration and volume differences that can be negative. Signed contributions are retained separately; ratios to a signed mean can exceed one or be negative when effects cancel, and are not probabilities.

| Fixed event | Uniform-weighted mean ΔH (mm) | Largest outcome contributor | Top 1 / 5 / 10 share of ΔH contribution | Outcome-concentration diagnostic | Largest leave-one-out mean change (mm) |
|---|---:|---|---:|---:|---:|
| September 2022 benchmark | 107.881 | `mc_00018` | 7.2% / 28.8% / 46.1% | 32.25 | 2.557 |
| August 2023 | 117.158 | `mc_00018` | 7.7% / 30.2% / 49.4% | 28.72 | 3.353 |
| July 2022 | 97.922 | `mc_00018` | 8.0% / 31.1% / 50.8% | 27.48 | 3.066 |
| August 2019 | 38.840 | `mc_00018` | 13.4% / 49.7% / 73.9% | 14.42 | 3.432 |
| August 2024 | 59.616 | `mc_00018` | 11.4% / 43.5% / 65.9% | 17.72 | 4.002 |

No single observation carries a majority of the total input weight or of a peak-stage mean. This rules out literal single-point dominance in this saved set; it is not a general reliability certificate. Outcome concentration is stronger for storms with many exact-zero peak differences, since only the remaining draws contribute to the mean. In 2019, the five largest outcome contributors supply about half the weighted mean, and ten supply nearly three-quarters; in 2024, ten supply about two-thirds. A relatively small subset therefore materially drives the weaker-event means even though the weight-only ESS is about 44. Weight ESS alone would understate this caution. Leave-one-out calculations merely remove one existing contribution and renormalize the remaining weights; they are influence diagnostics, not new model results or confidence intervals.

`outcome_contribution_summary.csv` includes ΔH, ΔQ, maximum-rise differences, duration above the initial reference, excess outlet volume above the initial discharge, total outlet-volume differences, and the exact-zero indicator. The zero indicator has identically zero contributions in the benchmark; its outcome-concentration diagnostic is undefined, not zero or 64. `top_weight_outcome_contributions.csv` reports signed and absolute contributions for the same top-weight samples in every event. `per_sample_outcome_contributions.csv` keeps every underlying contribution.

## Interpretation limits

Weight ESS is 1/Σw_i². An outcome-specific concentration diagnostic is 1/Σv_i², with v_i=|w_i y_i|/Σ|w_j y_j| when the denominator is nonzero. These are concentration diagnostics, not literal numbers of new iid target-distribution draws. An apparently adequate ESS does not establish that unsampled influential regions or rare tails were covered. Their formulation and limitations follow [Owen, importance-sampling chapter, Section 9.3](https://artowen.su.domains/mc/Ch-var-is.pdf). No ESS-based binomial intervals, weighted confidence bands, significance claims or sampling-convergence claims are made.

All alternative-law estimates remain exploratory because the uniform law is not calibrated. Positive weights cannot reveal a negative peak difference when all stored peak differences are nonnegative; that arithmetic does not validate the law or rule out a negative effect at an unvisited parameter vector. The five storms are fixed cases and reuse the same vectors, so they do not provide 320 independent draws of historical weather. No pooled weather-risk statement is made. Initial-reference metrics are not flood-threshold metrics, and the retained additional 5,500 m³ of wetland initial water remains part of the structural comparison.

The narrower new antecedent study answers a different question about initial conditions. This weight audit does not turn its designed 12 pairs into Monte Carlo samples, nor attach the old 64 weights to those pairs.

## Reproduce

From the repository root: `python D/d2_analysis/antecedent_conditions/statistics/verify_importance_weights.py`. The script reads existing CSV/JSON files and writes only this statistics directory. `audit.json` records input hashes and package versions. No applicable AGENTS.md was found in project/ancestor paths when this task began.
