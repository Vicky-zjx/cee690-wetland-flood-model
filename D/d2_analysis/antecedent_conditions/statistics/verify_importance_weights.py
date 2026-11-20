#!/usr/bin/env python3
"""Independently audit existing 64-draw importance weights and contributions.

Uses only saved five-event results. No model import, ODE, random draw or fit.
Outputs belong exclusively to antecedent_conditions/statistics.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
sys.dont_write_bytecode = True
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
D2 = HERE.parents[1]
PARAM = D2 / "parameter_uncertainty"
HIST = D2 / "historical_storm_robustness"
PRIOR = D2 / "robustness_checks/tables"
FACTORS = ["m_K_u", "m_routing", "m_R_r"]
METRICS = ["delta_peak_stage_m", "delta_peak_discharge_m3_s", "delta_maximum_rise_m",
           "delta_initial_level_exceedance_h", "delta_initial_discharge_excess_m3",
           "delta_total_outlet_volume_m3", "exact_zero_peak_effect_indicator"]
N = 64
TOP_K = [1, 5, 10]
SOURCE_HASHES = {}
LABELS = {"benchmark_september_2022": "September 2022 benchmark",
          "historical_20230813T0800Z": "August 2023",
          "historical_20220704T2100Z": "July 2022",
          "historical_20190825T0200Z": "August 2019",
          "historical_20240804T1500Z": "August 2024"}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def record(p):
    SOURCE_HASHES[str(p.relative_to(D2))] = sha(p)
    return p


def read_csv(p):
    return pd.read_csv(record(p), float_precision="round_trip")


def save_csv(rows, name):
    df = rows if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows)
    df.to_csv(HERE / name, index=False)
    return df


def save_json(obj, name):
    (HERE / name).write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def main():
    plan = json.loads(record(PARAM / "results/study_plan.json").read_text())
    draws = pd.DataFrame(plan["mc"]).sort_values("draw_index").reset_index(drop=True)
    assert len(draws) == N and draws.sample_id.is_unique
    m = draws[FACTORS].to_numpy()
    assert ((m >= .5) & (m <= 2)).all()
    # Source q: independent log-uniform; alternative p: independent uniform.
    # Keep the normalizing constants in p/q, even though they cancel in w.
    log_q = -np.log(m).sum(axis=1) - 3*np.log(np.log(4))
    log_p = -3*np.log(1.5)
    ratio = np.exp(log_p - log_q)
    direct_ratio = ((np.log(4)/1.5)*m).prod(axis=1)
    np.testing.assert_allclose(ratio, direct_ratio, rtol=1e-15, atol=1e-15)
    w = ratio / ratio.sum()
    ess = float(1/np.dot(w,w))
    ess_raw = float(ratio.sum()**2/np.dot(ratio,ratio))
    assert abs(ess - ess_raw) < 1e-12
    by_weight = np.argsort(-w, kind="stable")
    rank_weight = np.empty(N,dtype=int); rank_weight[by_weight] = np.arange(1,N+1)
    prior_weights = read_csv(PRIOR / "statistical_law_weights.csv")
    prior_weights = prior_weights[prior_weights.law == "uniform_multipliers"].set_index("sample_id").loc[draws.sample_id]
    weight_error = float(np.max(abs(w - prior_weights.normalized_weight.to_numpy())))
    assert weight_error < 1e-14
    prior_diag = read_csv(PRIOR / "statistical_law_diagnostics.csv").set_index("law").loc["uniform_multipliers"]
    assert abs(ess-prior_diag.weight_ess_diagnostic) < 1e-12
    weight_rows=[]
    for i in by_weight:
        weight_rows.append({"rank_by_weight":int(rank_weight[i]),"sample_id":draws.sample_id.iloc[i],
            "draw_index":int(draws.draw_index.iloc[i]),**{f:float(m[i,j]) for j,f in enumerate(FACTORS)},
            "density_ratio_p_over_q":float(ratio[i]),"normalized_weight":float(w[i]),
            "cumulative_weight_ranked_largest_first":float(w[by_weight[:rank_weight[i]]].sum())})
    save_csv(weight_rows,"weights_ranked.csv")
    mass_rows=[{"top_k":k,"normalized_weight_mass":float(w[by_weight[:k]].sum()),
                "equal_weight_mass_for_k_of_64":k/N,
                "sample_ids":";".join(draws.sample_id.iloc[by_weight[:k]])} for k in TOP_K]
    save_csv(mass_rows,"largest_weight_mass.csv")
    summary={"original_parameter_vectors":N,"original_seed":plan["seed"],
        "proposal":"Independent log-uniform multipliers [0.5,2] in the three approved categories",
        "target":"Exploratory independent uniform multipliers [0.5,2], not calibrated",
        "ratio":"p/q = product_j[m_j ln(4)/1.5]",
        "normalized_weight":"(p/q)_i / sum_j (p/q)_j",
        "ess_diagnostic":ess,"ess_fraction_of_original_64":ess/N,
        "maximum_normalized_weight":float(w.max()),"maximum_weight_sample_id":draws.sample_id.iloc[int(np.argmax(w))],
        "minimum_normalized_weight":float(w.min()),"mean_raw_ratio":float(ratio.mean()),
        "top_weight_mass":mass_rows,"prior_normalized_weights_max_absolute_difference":weight_error,
        "prior_ess_absolute_difference":float(abs(ess-prior_diag.weight_ess_diagnostic)),
        "ess_is_confidence_interval_sample_size":False,
        "new_parameter_vectors":0,"new_ode_runs":0}
    save_json(summary,"weight_summary.json")

    manifest=json.loads(record(HIST / "events/manifest.json").read_text())
    prior_out=read_csv(PRIOR / "statistical_law_outcomes.csv")
    prior_out=prior_out[prior_out.law=="uniform_multipliers"].set_index(["event_id","metric"])
    outcome_rows=[];top_rows=[];sample_rows=[];accounting=[];comparison=[]
    for event in manifest["events"]:
        eid=event["event_id"]
        f=read_csv(HIST / "results/ensemble" / eid / "paired_results.csv")
        assert f.sample_id.is_unique and len(f)==N
        f=f.set_index("sample_id").loc[draws.sample_id].reset_index()
        np.testing.assert_allclose(f[FACTORS],m,rtol=0,atol=2e-15)
        valid=f.paired_valid.astype(bool)
        accounting.append({"event_id":eid,"planned_pairs":N,"saved_pairs":len(f),"valid_pairs":int(valid.sum()),
            "invalid_pairs":int((~valid).sum()),"wetland_unsafe":int(f.wetland_unsafe.sum()),
            "no_wetland_unsafe":int(f.no_wetland_unsafe.sum()),"excluded_for_zero_or_negative_outcome":0})
        save_csv(accounting,"saved_pair_accounting.csv")
        if not valid.all(): raise RuntimeError(f"{eid}: incomplete original matrix; stopped without filtering failed rows")
        f["exact_zero_peak_effect_indicator"]=(f.delta_peak_stage_m==0).astype(float)
        for metric in METRICS:
            y=f[metric].to_numpy(dtype=float)
            assert np.isfinite(y).all()
            c=w*y; mean=float(c.sum()); abs_c=abs(c); abs_total=float(abs_c.sum())
            by_contrib=np.argsort(-abs_c,kind="stable")
            ranks=np.empty(N,dtype=int);ranks[by_contrib]=np.arange(1,N+1)
            positive=float(c[c>0].sum());negative=float(c[c<0].sum())
            if abs_total>0:
                f_weights=abs_c/abs_total
                f_ess=float(1/np.dot(f_weights,f_weights))
                cancellation=float(1-abs(mean)/abs_total)
                max_share=float(f_weights.max())
            else:
                f_ess=cancellation=max_share=None
            loo=(mean-c)/(1-w)
            max_loo=int(np.argmax(abs(loo-mean)))
            row={"event_id":eid,"metric":metric,"n_saved_pairs":N,"weighted_mean":mean,
                "positive_weighted_contribution_sum":positive,"negative_weighted_contribution_sum":negative,
                "absolute_contribution_sum":abs_total,"contribution_concentration_diagnostic":f_ess,
                "maximum_absolute_outcome_contribution_share":max_share,
                "largest_outcome_contributor_sample_id":draws.sample_id.iloc[int(by_contrib[0])] if abs_total>0 else "",
                "nonzero_outcome_count":int((y!=0).sum()),"cancellation_fraction":cancellation,
                "maximum_leave_one_out_absolute_mean_change":float(abs(loo[max_loo]-mean)),
                "maximum_leave_one_out_change_sample_id":draws.sample_id.iloc[max_loo],
                "maximum_leave_one_out_relative_change":float(abs(loo[max_loo]-mean)/abs(mean)) if mean!=0 else None,
                "loo_diagnostic_not_a_new_simulation":True,"all_zero_contributions":abs_total==0}
            for k in TOP_K:
                row[f"top_{k}_absolute_outcome_contribution_share"]=float(abs_c[by_contrib[:k]].sum()/abs_total) if abs_total>0 else None
                ix=by_weight[:k]
                top_rows.append({"event_id":eid,"metric":metric,"top_k_ranked_by_weight":k,
                    "weight_mass":float(w[ix].sum()),"signed_contribution_to_weighted_mean":float(c[ix].sum()),
                    "absolute_contribution_share":float(abs_c[ix].sum()/abs_total) if abs_total>0 else None,
                    "signed_contribution_divided_by_mean":float(c[ix].sum()/mean) if mean!=0 else None,
                    "signed_ratio_is_not_a_probability":True,"sample_ids":";".join(draws.sample_id.iloc[ix])})
            outcome_rows.append(row)
            for i in range(N):
                sample_rows.append({"event_id":eid,"metric":metric,"sample_id":draws.sample_id.iloc[i],
                    "normalized_weight":float(w[i]),"rank_by_weight":int(rank_weight[i]),
                    "saved_outcome":float(y[i]),"signed_weighted_contribution":float(c[i]),
                    "rank_by_absolute_outcome_contribution":int(ranks[i]),
                    "absolute_outcome_contribution_share":float(abs_c[i]/abs_total) if abs_total>0 else None})
            old_metric="delta_maximum_rise_above_initial_m" if metric=="delta_maximum_rise_m" else metric
            old_value=float(prior_out.loc[(eid,old_metric),"self_normalized_weighted_mean"]) if metric!="exact_zero_peak_effect_indicator" else float(prior_out.loc[(eid,"delta_peak_stage_m"),"weighted_zero_mass"])
            difference=abs(mean-old_value)
            assert difference<1e-8
            comparison.append({"event_id":eid,"metric":metric,"recomputed":mean,"previous":old_value,"absolute_difference":difference})
    outcomes=save_csv(outcome_rows,"outcome_contribution_summary.csv")
    save_csv(top_rows,"top_weight_outcome_contributions.csv")
    save_csv(sample_rows,"per_sample_outcome_contributions.csv")
    save_csv(comparison,"prior_result_verification.csv")

    lines=["# Independent importance-weight audit", "",
        "This is a saved-result statistical supplement to the antecedent-condition study. It audits the **previous five-event, 64-vector ensemble**, not the new antecedent trajectories. No ODE, new physical parameter vector, random draw, fit or model change is performed. All original results and parameter values remain unchanged.", "",
        "The proposal density is q(m)=∏[1/(m_j ln 4)] on [0.5,2]³. The exploratory alternative is independent uniform, p(m)=1/(1.5³), so p/q=∏[m_j ln(4)/1.5]. Normalizing these density ratios gives the same per-draw weights for every event. The three categories remain infiltration, common routing, and river resistance; no routing rate is split into a new independent parameter.", "",
        f"The independently recomputed weight ESS is **{ess:.8f}/64** ({ess/N:.2%}). The largest normalized weight is **{w.max():.4%}**, for `{summary['maximum_weight_sample_id']}`. The largest 1, 5 and 10 weights contain respectively **{mass_rows[0]['normalized_weight_mass']:.4%}, {mass_rows[1]['normalized_weight_mass']:.4%}, and {mass_rows[2]['normalized_weight_mass']:.4%}** of the weight. Equal weighting would give 1.5625%, 7.8125%, and 15.625%. The smallest weight is {w.min():.4%}. The mean unnormalized p/q is {ratio.mean():.6f}, rather than exactly one in this finite sample; self-normalization divides by the realized sum.", "",
        f"Two independent ESS formulas agree. The largest difference from the saved normalized weights is {weight_error:.3g}; all previously saved reweighted means and zero masses agree within the tolerances recorded in `prior_result_verification.csv`. All 320 original pairs remain valid and upstream-safety-feasible; no zero or unfavorable metric is excluded.", "",
        "## Weight concentration versus outcome concentration", "",
        "For an outcome y_i, its signed contribution to the self-normalized mean is c_i=w_i y_i. The tables distinguish samples ranked by **weight** from samples ranked by **absolute outcome contribution**. The latter share is |c_i|/Σ|c_j|, so it remains interpretable for duration and volume differences that can be negative. Signed contributions are retained separately; ratios to a signed mean can exceed one or be negative when effects cancel, and are not probabilities.", "",
        "| Fixed event | Uniform-weighted mean ΔH (mm) | Largest outcome contributor | Top 1 / 5 / 10 share of ΔH contribution | Outcome-concentration diagnostic | Largest leave-one-out mean change (mm) |",
        "|---|---:|---|---:|---:|---:|",
    ]
    peaks=outcomes[outcomes.metric=="delta_peak_stage_m"]
    for r in peaks.itertuples():
        shares=" / ".join(f"{getattr(r,f'top_{k}_absolute_outcome_contribution_share'):.1%}" for k in TOP_K)
        lines.append(f"| {LABELS[r.event_id]} | {r.weighted_mean*1000:.3f} | `{r.largest_outcome_contributor_sample_id}` | {shares} | {r.contribution_concentration_diagnostic:.2f} | {r.maximum_leave_one_out_absolute_mean_change*1000:.3f} |")
    lines += ["", "No single observation carries a majority of the total input weight or of a peak-stage mean. This rules out literal single-point dominance in this saved set; it is not a general reliability certificate. Outcome concentration is stronger for storms with many exact-zero peak differences, since only the remaining draws contribute to the mean. In 2019, the five largest outcome contributors supply about half the weighted mean, and ten supply nearly three-quarters; in 2024, ten supply about two-thirds. A relatively small subset therefore materially drives the weaker-event means even though the weight-only ESS is about 44. Weight ESS alone would understate this caution. Leave-one-out calculations merely remove one existing contribution and renormalize the remaining weights; they are influence diagnostics, not new model results or confidence intervals.", "",
        "`outcome_contribution_summary.csv` includes ΔH, ΔQ, maximum-rise differences, duration above the initial reference, excess outlet volume above the initial discharge, total outlet-volume differences, and the exact-zero indicator. The zero indicator has identically zero contributions in the benchmark; its outcome-concentration diagnostic is undefined, not zero or 64. `top_weight_outcome_contributions.csv` reports signed and absolute contributions for the same top-weight samples in every event. `per_sample_outcome_contributions.csv` keeps every underlying contribution.", "",
        "## Interpretation limits", "",
        "Weight ESS is 1/Σw_i². An outcome-specific concentration diagnostic is 1/Σv_i², with v_i=|w_i y_i|/Σ|w_j y_j| when the denominator is nonzero. These are concentration diagnostics, not literal numbers of new iid target-distribution draws. An apparently adequate ESS does not establish that unsampled influential regions or rare tails were covered. Their formulation and limitations follow [Owen, importance-sampling chapter, Section 9.3](https://artowen.su.domains/mc/Ch-var-is.pdf). No ESS-based binomial intervals, weighted confidence bands, significance claims or sampling-convergence claims are made.", "",
        "All alternative-law estimates remain exploratory because the uniform law is not calibrated. Positive weights cannot reveal a negative peak difference when all stored peak differences are nonnegative; that arithmetic does not validate the law or rule out a negative effect at an unvisited parameter vector. The five storms are fixed cases and reuse the same vectors, so they do not provide 320 independent draws of historical weather. No pooled weather-risk statement is made. Initial-reference metrics are not flood-threshold metrics, and the retained additional 5,500 m³ of wetland initial water remains part of the structural comparison.", "",
        "The narrower new antecedent study answers a different question about initial conditions. This weight audit does not turn its designed 12 pairs into Monte Carlo samples, nor attach the old 64 weights to those pairs.", "",
        "## Reproduce", "",
        "From the repository root: `python D/d2_analysis/antecedent_conditions/statistics/verify_importance_weights.py`. The script reads existing CSV/JSON files and writes only this statistics directory. `audit.json` records input hashes and package versions. No applicable AGENTS.md was found in project/ancestor paths when this task began.", ""]
    (HERE/"README.md").write_text("\n".join(lines))
    for path,digest in SOURCE_HASHES.items(): assert sha(D2/path)==digest,path
    save_json({"source_hashes_unchanged":True,"source_sha256":SOURCE_HASHES,"script_sha256":sha(__file__),
        "numpy_version":np.__version__,"pandas_version":pd.__version__,"new_physical_vectors":0,"new_ode_runs":0,
        "saved_events":5,"saved_pairs":320,"weight_summary":summary,
        "maximum_prior_outcome_error":max(x["absolute_difference"] for x in comparison)},"audit.json")
    print(json.dumps({"ESS":ess,"top_weight_mass":{str(k):float(w[by_weight[:k]].sum()) for k in TOP_K},
        "new_ode_runs":0,"source_hashes_unchanged":True,"output_directory":str(HERE)},indent=2))


if __name__ == "__main__":
    main()
