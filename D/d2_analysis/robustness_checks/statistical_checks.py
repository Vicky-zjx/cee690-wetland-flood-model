#!/usr/bin/env python3
"""Read-only audit of saved five-event outputs; write only robustness_checks.

No physical model, solver, parameter sampler, or surrogate is imported or run.
Random permutations only select subsets of the 64 already evaluated vectors.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

sys.dont_write_bytecode = True
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/cee690-robustness-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import beta, binomtest, multivariate_normal, norm, spearmanr

HERE = Path(__file__).resolve().parent
D2 = HERE.parent
HIST = D2 / "historical_storm_robustness"
PARAM = D2 / "parameter_uncertainty"
TABLES, PLOTS = HERE / "tables", HERE / "plots"
FACTORS = ["m_K_u", "m_routing", "m_R_r"]
N = 64
SUBSET_SEED = 69020261119
SUBSET_REPETITIONS = 1000
SIZES = [16, 32, 48, 64]
LABELS = {
    "benchmark_september_2022": "Sep 2022\nbenchmark",
    "historical_20230813T0800Z": "Aug 2023",
    "historical_20220704T2100Z": "Jul 2022",
    "historical_20190825T0200Z": "Aug 2019",
    "historical_20240804T1500Z": "Aug 2024",
}
METRICS = [f"{case}_{metric}" for case in ["wetland", "no_wetland"] for metric in [
    "peak_stage_m", "peak_discharge_m3_s", "maximum_rise_above_initial_m",
    "initial_level_exceedance_h", "initial_discharge_excess_m3", "total_outlet_volume_m3"]]
METRICS += ["delta_peak_stage_m", "delta_peak_discharge_m3_s", "delta_maximum_rise_above_initial_m",
            "delta_initial_level_exceedance_h", "delta_initial_discharge_excess_m3", "delta_total_outlet_volume_m3"]
LAW_LABELS = {
    "baseline_loguniform": "Independent log-uniform",
    "uniform_multipliers": "Independent uniform",
    "copula_negative": "Log-uniform; latent rho = -0.25",
    "copula_positive": "Log-uniform; latent rho = +0.25",
}
LAW_COLORS = ["#333333", "#0072B2", "#D55E00", "#009E73"]
SOURCE_HASHES = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def record(path):
    SOURCE_HASHES[str(Path(path).relative_to(D2))] = sha(path)
    return path


def read_csv(path):
    return pd.read_csv(record(path), float_precision="round_trip")


def write_csv(rows, name):
    df = rows if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows)
    df.to_csv(TABLES / f"statistical_{name}.csv", index=False)
    return df


def write_json(data, name):
    (TABLES / f"statistical_{name}.json").write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def cp(k, n):
    k = np.asarray(k)
    lo = np.zeros(k.shape, dtype=float)
    hi = np.ones(k.shape, dtype=float)
    lo[k > 0] = beta.ppf(.025, k[k > 0], n - k[k > 0] + 1)
    hi[k < n] = beta.ppf(.975, k[k < n] + 1, n - k[k < n])
    return lo, hi


def load():
    manifest_path = HIST / "events/manifest.json"
    manifest = json.loads(record(manifest_path).read_text())
    plan = json.loads(record(PARAM / "results/study_plan.json").read_text())
    draws = pd.DataFrame(plan["mc"]).sort_values("draw_index").reset_index(drop=True)
    assert len(draws) == N and not draws.sample_id.duplicated().any()
    data, accounting = {}, []
    for event in manifest["events"]:
        eid = event["event_id"]
        f = read_csv(HIST / "results/ensemble" / eid / "paired_results.csv")
        f = f.set_index("sample_id").loc[draws.sample_id].reset_index()
        assert len(f) == N and not f.sample_id.duplicated().any()
        np.testing.assert_allclose(f[FACTORS], draws[FACTORS], rtol=0, atol=2e-15)
        # Source CSV's shorter alias is mapped without changing source data.
        f["delta_maximum_rise_above_initial_m"] = f["delta_maximum_rise_m"]
        valid = f.paired_valid.astype(bool)
        accounting.append({"event_id": eid, "planned_pairs": N, "saved_pairs": len(f),
            "valid_pairs": int(valid.sum()), "invalid_pairs": int((~valid).sum()),
            "wetland_unsafe": int(f.wetland_unsafe.sum()), "no_wetland_unsafe": int(f.no_wetland_unsafe.sum()),
            "different_actual_gate_histories": int((~f.paired_gate_history_agreement.astype(bool)).sum()),
            "weather_matches": int(f.paired_weather_inputs_exact.sum()), "excluded_for_risk_sign": 0})
        # This specific finite-population audit requires the published complete matrix.
        # Stop rather than silently drop any failure and change the 64-vector population.
        if not valid.all() or not f.paired_weather_inputs_exact.all():
            write_csv(accounting, "accounting")
            raise RuntimeError(f"{eid}: incomplete/invalid matrix; audit stopped without filtering rows")
        assert np.isfinite(f[METRICS].to_numpy()).all()
        data[eid] = f
    return data, draws, write_csv(accounting, "accounting")


def check_ecdf(data):
    check_rows = []
    k = np.arange(N + 1)
    lo, hi = cp(k, N)
    for i in k:
        official = binomtest(int(i), N).proportion_ci(.95, method="exact")
        check_rows.append({"k": int(i), "n": N, "beta_lower": lo[i], "beta_upper": hi[i],
                           "scipy_exact_lower": official.low, "scipy_exact_upper": official.high,
                           "maximum_error": max(abs(lo[i]-official.low), abs(hi[i]-official.high))})
    checks = write_csv(check_rows, "cp_scipy_checks")
    assert checks.maximum_error.max() < 1e-10
    saved = read_csv(HIST / "tables/fixed_grid_pointwise_cdf_intervals.csv.gz")
    grid_rows = []
    for (eid, metric), g in saved.groupby(["event_id", "metric"], sort=False):
        vals = np.sort(data[eid][metric].to_numpy())
        counts = np.searchsorted(vals, g.threshold.to_numpy(), side="right")
        low, high = cp(counts, N)
        count_error = int(np.max(abs(counts - g.count_leq.to_numpy())))
        cdf_error = float(np.max(abs(counts / N - g.empirical_cdf.to_numpy())))
        band_error = float(max(np.max(abs(low-g.pointwise_lower)), np.max(abs(high-g.pointwise_upper))))
        assert count_error == 0 and cdf_error < 1e-14 and band_error < 1e-12
        grid_rows.append({"event_id": eid, "metric": metric, "checked_fixed_thresholds": len(g),
                          "maximum_count_error": count_error, "maximum_cdf_error": cdf_error,
                          "maximum_cp_endpoint_error": band_error})
    write_csv(grid_rows, "saved_grid_checks")
    original = read_csv(HIST / "tables/official_pointwise_ecdf_arrays.csv")
    rank_rows = []
    for (eid, metric), g in original.groupby(["event_id", "metric"], sort=False):
        g = g.sort_values("rank")
        vals = np.sort(data[eid][metric].to_numpy())
        assert len(vals) == len(g)
        value_error = float(np.max(abs(vals-g.sorted_value.to_numpy())))
        rank_error = float(np.max(abs(np.arange(1, N+1)/(N+1)-g.official_rank_position.to_numpy())))
        assert value_error < 1e-12 and rank_error < 1e-14
        rank_rows.append({"event_id": eid, "metric": metric, "maximum_sorted_value_error": value_error,
            "maximum_rank_position_error": rank_error, "unique_values": len(np.unique(vals)),
            "is_count_over_n_ecdf": False, "ties_require_count_over_n": len(np.unique(vals)) < N})
    write_csv(rank_rows, "official_rank_checks")
    atom_rows = []
    for eid, f in data.items():
        y = f.delta_peak_stage_m.to_numpy()
        z = y == 0.0
        low, high = cp(np.array([z.sum()]), N)
        both_t0 = ((f.wetland_peak_stage_time_s == 0) & (f.no_wetland_peak_stage_time_s == 0))
        both_h0 = ((f.wetland_peak_stage_m == f.wetland_initial_stage_m) &
                   (f.no_wetland_peak_stage_m == f.no_wetland_initial_stage_m))
        assert np.array_equal(z, both_t0 & both_h0)
        atom_rows.append({"event_id": eid, "n": N, "positive_count": int((y > 0).sum()),
            "exact_zero_count": int(z.sum()), "negative_count": int((y < 0).sum()),
            "both_peaks_at_initial_condition": int((both_t0 & both_h0).sum()),
            "zero_mass": float(z.mean()), "zero_mass_cp95_lower": low[0], "zero_mass_cp95_upper": high[0],
            "count_at_or_below_zero": int((y<=0).sum()), "empirical_cdf_at_zero": float(np.mean(y<=0)),
            "cdf_at_zero_cp95_lower": float(cp(np.array([(y<=0).sum()]),N)[0][0]),
            "cdf_at_zero_cp95_upper": float(cp(np.array([(y<=0).sum()]),N)[1][0]),
            "interval_contains_half": bool(low[0] <= .5 <= high[0]), "empirical_median_m": float(np.median(y))})
    atoms = write_csv(atom_rows, "zero_atoms")
    assert atoms.set_index("event_id").loc["historical_20190825T0200Z", "exact_zero_count"] == 37
    assert atoms.set_index("event_id").loc["historical_20240804T1500Z", "exact_zero_count"] == 30
    # Recompute every saved outcome summary, not only the peak-benefit metric.
    summaries = read_csv(HIST / "tables/outcome_summaries.csv")
    errors = []
    for r in summaries.itertuples():
        y = data[r.event_id][r.metric].to_numpy()
        computed = {"mean": np.mean(y), "median": np.median(y), "minimum": np.min(y), "maximum": np.max(y),
                    "sample_sd": np.std(y, ddof=1), "q05": np.quantile(y,.05), "q95": np.quantile(y,.95)}
        for stat, value in computed.items():
            delta = abs(value - getattr(r, stat))
            assert delta < 1e-8
            errors.append({"event_id": r.event_id, "metric": r.metric, "statistic": stat,
                           "recomputed": value, "saved": getattr(r,stat), "absolute_error": delta})
        if r.metric.startswith("delta_"):
            assert (int((y>0).sum()),int((y==0).sum()),int((y<0).sum())) == (r.positive_count,r.zero_count,r.negative_count)
    write_csv(errors, "full_summary_checks")
    uncertainty=[]
    for eid,f in data.items():
        for metric in METRICS:
            y=f[metric].to_numpy(); sd=float(np.std(y,ddof=1))
            uncertainty.append({"event_id":eid,"metric":metric,"n":N,
                "mean":float(y.mean()),"sample_sd_of_physical_outcomes":sd,
                "estimated_monte_carlo_standard_error_of_mean":sd/np.sqrt(N),
                "outcome_q05":float(np.quantile(y,.05)),"outcome_q95":float(np.quantile(y,.95)),
                "scope":"Baseline iid assumed parameter law, one fixed event; across-event estimates are dependent",
                "outcome_quantiles_are_confidence_limits":False,"mean_confidence_interval_computed":False})
    write_csv(uncertainty,"sampling_uncertainty")
    return atoms, saved, checks.maximum_error.max()


def subset_stats(f, ix):
    y = f.delta_peak_stage_m.to_numpy()[ix]
    return {"mean_delta_peak_stage_m": float(y.mean()), "median_delta_peak_stage_m": float(np.median(y)),
            "zero_fraction": float(np.mean(y == 0)), "q95_delta_peak_stage_m": float(np.quantile(y,.95)),
            "mean_delta_peak_discharge_m3_s": float(f.delta_peak_discharge_m3_s.to_numpy()[ix].mean()),
            "mean_delta_initial_level_exceedance_h": float(f.delta_initial_level_exceedance_h.to_numpy()[ix].mean()),
            "mean_delta_total_outlet_volume_m3": float(f.delta_total_outlet_volume_m3.to_numpy()[ix].mean())}


def stability(data):
    rng = np.random.default_rng(SUBSET_SEED)
    permutations = np.array([rng.permutation(N) for _ in range(SUBSET_REPETITIONS)])
    write_json({"seed": SUBSET_SEED, "rng": "numpy.default_rng / PCG64", "numpy_version": np.__version__,
        "repetitions": SUBSET_REPETITIONS, "selection": "Prefixes of each random permutation, without replacement",
        "same_indices_across_events": True, "sizes": SIZES,
        "permutations_of_original_draw_index": permutations.tolist(),
        "interpretation": "Subset diagnostics of existing finite 64-draw population; no new parameter values or model runs"}, "subset_design")
    prefix_rows, random_rows = [], []
    for eid, f in data.items():
        for n in SIZES:
            prefix_rows.append({"event_id": eid, "n": n, **subset_stats(f, np.arange(n))})
            count = SUBSET_REPETITIONS if n < N else 1
            for rep in range(count):
                ix = permutations[rep,:n] if n < N else np.arange(N)
                random_rows.append({"event_id": eid, "n": n, "replicate": rep, **subset_stats(f, ix)})
    prefixes = write_csv(prefix_rows, "nested_prefixes")
    randoms = write_csv(random_rows, "subset_replicates")
    summaries = []
    for (eid,n), g in randoms.groupby(["event_id","n"], sort=False):
        for metric in g.columns[3:]:
            y = g[metric].to_numpy()
            summaries.append({"event_id": eid, "n": n, "metric": metric, "subset_repetitions": len(g),
                "mean_over_subsets": y.mean(), "q025": np.quantile(y,.025), "median": np.median(y),
                "q975": np.quantile(y,.975), "minimum": y.min(), "maximum": y.max(),
                "not_a_population_confidence_interval": True})
    return prefixes, write_csv(summaries, "subset_summaries")


def weighted_quantile(y,w,p):
    ix=np.argsort(y,kind="stable")
    cumulative=np.cumsum(w[ix]); cumulative[-1]=1.0
    return float(y[ix[min(np.searchsorted(cumulative,p,side="left"),len(ix)-1)]])


def laws(data, draws):
    m = draws[FACTORS].to_numpy()
    u = np.log(m/.5) / np.log(4)
    assert ((u>0)&(u<1)).all()
    z = norm.ppf(u)
    raw = {"baseline_loguniform": np.ones(N),
           "uniform_multipliers": ((np.log(4)/1.5)*m).prod(axis=1)}
    # Same original marginals. Only latent normal correlation K_u/routing changes.
    for label,rho in [("copula_negative",-.25),("copula_positive",.25)]:
        R = np.eye(3); R[0,1] = R[1,0] = rho
        logw = multivariate_normal.logpdf(z, mean=np.zeros(3), cov=R) - norm.logpdf(z).sum(axis=1)
        closed = -.5*np.log(np.linalg.det(R)) - .5*np.einsum("ni,ij,nj->n",z,np.linalg.inv(R)-np.eye(3),z)
        np.testing.assert_allclose(logw,closed,rtol=0,atol=4e-15)
        raw[label] = np.exp(logw)
    diagnostics, weight_rows, rows, marginal_rows = [], [], [], []
    normalized = {label: w/w.sum() for label,w in raw.items()}
    for label,w in normalized.items():
        ess = 1/(w*w).sum()
        diagnostics.append({"law":label,"n_original_draws":N,"weight_ess_diagnostic":ess,
            "ess_fraction":ess/N,"maximum_normalized_weight":w.max(),"minimum_normalized_weight":w.min(),
            "mean_raw_density_ratio":raw[label].mean(),"sum_normalized_weights":w.sum(),
            "target_rho_latent_Ku_routing": -.25 if label=="copula_negative" else .25 if label=="copula_positive" else 0,
            "calibrated_distribution":False,"weighted_confidence_bands_computed":False})
        for i in range(N):
            weight_rows.append({"law":label,"sample_id":draws.sample_id.iloc[i],"draw_index":int(draws.draw_index.iloc[i]),
                               "raw_density_ratio":raw[label][i],"normalized_weight":w[i]})
        for j,factor in enumerate(FACTORS):
            exact_mean = 1.25 if label=="uniform_multipliers" else 1.5/np.log(4)
            marginal_rows.append({"law":label,"factor":factor,"reweighted_sample_mean":np.dot(w,m[:,j]),
                "target_marginal_mean":exact_mean,"difference":np.dot(w,m[:,j])-exact_mean})
        for eid,f in data.items():
            for metric in METRICS:
                y=f[metric].to_numpy()
                rows.append({"law":label,"event_id":eid,"metric":metric,"n_original_draws":N,
                    "self_normalized_weighted_mean":np.dot(w,y),
                    "weighted_inverse_cdf_median":weighted_quantile(y,w,.5),
                    "weighted_inverse_cdf_q05":weighted_quantile(y,w,.05),
                    "weighted_inverse_cdf_q95":weighted_quantile(y,w,.95),
                    "weighted_positive_mass":w[y>0].sum(),"weighted_zero_mass":w[y==0].sum(),
                    "weighted_negative_mass":w[y<0].sum(),"weight_ess_diagnostic":ess,
                    "not_new_iid_target_draws":True})
    write_csv(weight_rows,"law_weights")
    write_csv(marginal_rows,"law_marginal_diagnostics")
    definitions={
        "baseline_loguniform": "q(m)=product[1/(m_j*ln(4))], 0.5<=m_j<=2; independent; all 64 weights equal",
        "uniform_multipliers": "p(m)=1/(1.5^3) on the same cube; p/q=product[m_j*ln(4)/1.5]",
        "copula_negative": "Same log-uniform marginals; u=ln(m/0.5)/ln(4), z=Phi^-1(u); R12=R21=-0.25, other off-diagonals zero",
        "copula_positive": "Same log-uniform marginals; u=ln(m/0.5)/ln(4), z=Phi^-1(u); R12=R21=+0.25, other off-diagonals zero",
        "copula_ratio": "p/q=det(R)^(-1/2)*exp[-z^T(R^-1-I)z/2]",
        "estimate": "sum_i normalized_weight_i * saved_output_i; normalized_weight=(p/q)/sum(p/q)",
        "why_dependence": "Infiltration and routing share the upland compartment. Neither sign nor magnitude is established; symmetric modest latent correlations are descriptive assumptions, not measured physical dependence.",
        "scope": "No samples created, no ODE run, no distribution fitted. Common routing multiplier remains one category; no extra independent rates.",
        "ess": "1/sum(normalized_weight^2), a weight-concentration diagnostic only. Does not replace n for binomial intervals or certify tail coverage.",
        "quantile": "Inverse weighted empirical CDF; differs from numpy's interpolated sample median at an even sample size.",
    }
    write_json(definitions,"law_definitions")
    return write_csv(rows,"law_outcomes"), write_csv(diagnostics,"law_diagnostics")


def parameter_patterns(data):
    saved = read_csv(HIST / "tables/parameter_rank_associations.csv")
    rows=[]
    for r in saved.itertuples():
        rho=float(spearmanr(data[r.event_id][r.factor],data[r.event_id][r.outcome]).statistic)
        assert abs(rho-r.spearman_rho)<1e-13
        rows.append({"event_id":r.event_id,"factor":r.factor,"outcome":r.outcome,"spearman_rho":rho,
            "absolute_difference_from_saved":abs(rho-r.spearman_rho),"p_value_reported":False,
            "interpretation":"Descriptive marginal rank association; correlated output conversion and interactions remain"})
    ranks=write_csv(rows,"rank_checks")
    grid=read_csv(PARAM / "results/sensitivity/paired_results.csv")
    assert len(grid)==27 and grid.paired_valid.all()
    endpoints=[]
    for factor in FACTORS:
        other=[c for c in FACTORS if c!=factor]
        for context,g in grid.groupby(other,sort=True):
            g=g.sort_values(factor)
            np.testing.assert_array_equal(g[factor].to_numpy(),[.5,1,2])
            low,mid,high=[g.iloc[i] for i in range(3)]
            for metric in ["delta_peak_stage_m","delta_peak_discharge_m3_s"]:
                endpoints.append({"event_id":"benchmark_september_2022","varied_factor":factor,
                    "fixed_factor_1":other[0],"fixed_value_1":context[0],"fixed_factor_2":other[1],"fixed_value_2":context[1],
                    "metric":metric,"value_at_half":low[metric],"value_at_one":mid[metric],"value_at_two":high[metric],
                    "endpoint_change_two_minus_half":high[metric]-low[metric],
                    "monotone_at_three_grid_points":bool((low[metric]<=mid[metric]<=high[metric]) or (low[metric]>=mid[metric]>=high[metric])),
                    "scope":"Existing benchmark grid only; nine controlled contexts per factor; no new runs"})
    end=write_csv(endpoints,"grid_endpoint_contrasts")
    summary=[]
    for (factor,metric),g in end.groupby(["varied_factor","metric"],sort=False):
        y=g.endpoint_change_two_minus_half
        summary.append({"factor":factor,"metric":metric,"contexts":len(g),"positive":int((y>0).sum()),
            "zero":int((y==0).sum()),"negative":int((y<0).sum()),"minimum_endpoint_change":y.min(),
            "maximum_endpoint_change":y.max(),"monotone_three_point_contexts":int(g.monotone_at_three_grid_points.sum())})
    return ranks,write_csv(summary,"grid_endpoint_summary")


def savefig(fig,name):
    fig.savefig(PLOTS/name,dpi=180,facecolor="white")
    plt.close(fig)


def figures(data,atoms,saved,prefixes,subsets,outcomes,diagnostics):
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,"text.usetex":False})
    fig,axs=plt.subplots(1,3,figsize=(14,4.5),gridspec_kw={"width_ratios":[1,1,1.08]})
    for ax,eid in zip(axs[:2],["historical_20190825T0200Z","historical_20240804T1500Z"]):
        y=data[eid].delta_peak_stage_m.to_numpy()*1000
        x,counts=np.unique(y,return_counts=True)
        F=np.cumsum(counts)/N
        ax.step(np.r_[min(-3,x[0]-1),x,x[-1]+5],np.r_[0,F,1],where="post",color="#0072B2",lw=2,label="Exact count / 64")
        g=saved[(saved.event_id==eid)&(saved.metric=="delta_peak_stage_m")]
        g=g[(g.threshold>=-.005)&(g.threshold<=y.max()/1000+.005)]
        ax.fill_between(g.threshold*1000,g.pointwise_lower,g.pointwise_upper,color="#0072B2",alpha=.16,
                        label="95% pointwise CP")
        k=int((y==0).sum()); ax.plot([0],[k/N],"o",color="#D55E00")
        ax.annotate(f"Exact zero: {k}/64\nF(0) = {k/N:.3f}",xy=(0,k/N),xytext=(.42,.20),textcoords="axes fraction",
                    arrowprops={"arrowstyle":"->","color":"#555"},fontsize=10)
        ax.set(title=LABELS[eid],xlabel="Peak-stage difference, no wetland − wetland (mm)",ylabel="Empirical cumulative fraction",ylim=(0,1.02),xlim=(-5,y.max()+5))
        ax.axvline(0,color="#999",lw=.7);ax.grid(alpha=.2);ax.legend(loc="lower right",fontsize=8)
    a=atoms.set_index("event_id").loc[list(data)]
    xx=np.arange(len(a)); yy=a.zero_mass.to_numpy()
    axs[2].errorbar(xx,yy,yerr=[yy-a.zero_mass_cp95_lower,a.zero_mass_cp95_upper-yy],fmt="o",capsize=4,color="#333")
    axs[2].axhline(.5,color="#D55E00",ls="--",lw=1,label="One-half")
    axs[2].set(xticks=xx,xticklabels=[LABELS[e] for e in data],ylim=(-.02,.82),ylabel="Mass at exactly zero peak benefit",title="Pointwise uncertainty at zero")
    for i,r in enumerate(a.itertuples()): axs[2].annotate(f"{r.exact_zero_count}/64",(i,r.zero_mass),xytext=(0,9),textcoords="offset points",ha="center",fontsize=9)
    axs[2].tick_params(axis="x",labelsize=8);axs[2].grid(axis="y",alpha=.2)
    fig.suptitle("Zero effects are retained: both whole-window river maxima occur at the initial state",fontsize=13,y=.98)
    fig.text(.5,.02,"Intervals describe finite-sample uncertainty under the assumed parameter law, at fixed thresholds; not simultaneous bands or weather-risk probabilities.",ha="center",fontsize=9)
    fig.tight_layout(rect=(0,.07,1,.91));savefig(fig,"04_ecdf_zero_atoms.png")

    fig,axs=plt.subplots(2,5,figsize=(15,6.5),sharex=True)
    for col,eid in enumerate(data):
        p=prefixes[prefixes.event_id==eid].sort_values("n")
        for row,(metric,factor,ylabel) in enumerate([("mean_delta_peak_stage_m",1000,"Mean peak-stage difference (mm)"),("zero_fraction",1,"Exact-zero fraction")]):
            s=subsets[(subsets.event_id==eid)&(subsets.metric==metric)].sort_values("n")
            ax=axs[row,col]
            ax.fill_between(s.n.to_numpy(),s.q025.to_numpy()*factor,s.q975.to_numpy()*factor,color="#56B4E9",alpha=.3,label="Middle 95% of subset results")
            ax.plot(s.n,s["median"]*factor,color="#0072B2",ls="--",lw=1.4,label="Median across subsets")
            ax.plot(p.n,p[metric]*factor,"o-",color="#D55E00",lw=1.6,label="Original-order prefix")
            ax.axhline(float(p.loc[p.n==64,metric].iloc[0])*factor,color="#555",lw=.8,ls=":")
            ax.set_xticks(SIZES);ax.grid(alpha=.18)
            if row==0: ax.set_title(LABELS[eid])
            else: ax.set_xlabel("Number of existing draws")
            if col==0: ax.set_ylabel(ylabel)
            if row==1: ax.set_ylim(-.025,.95)
    handles,labels=axs[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc="lower center",bbox_to_anchor=(.5,.045),ncol=3,frameon=False,fontsize=9)
    fig.suptitle("Sample-size diagnostics reuse the same 64 evaluated vectors",fontsize=14,y=.99)
    fig.text(.5,.013,"1,000 matched random permutations; prefixes sampled without replacement. Ribbons are conditional subset ranges, not confidence intervals or convergence proof.\nAt n = 64 every subset is the full existing set, so the ribbon collapses by construction.",ha="center",fontsize=9)
    fig.tight_layout(rect=(0,.13,1,.94));savefig(fig,"05_sample_size_stability.png")

    fig,axs=plt.subplots(1,3,figsize=(14,5),gridspec_kw={"width_ratios":[1.25,1.25,.85]})
    xx=np.arange(len(data))
    for label,color in zip(LAW_LABELS,LAW_COLORS):
        p=outcomes[(outcomes.law==label)&(outcomes.metric=="delta_peak_stage_m")].set_index("event_id").loc[list(data)]
        axs[0].plot(xx,p.self_normalized_weighted_mean*1000,"o-",color=color,label=LAW_LABELS[label],lw=1.5)
        axs[1].plot(xx,p.weighted_zero_mass,"o-",color=color,lw=1.5)
    for ax in axs[:2]:
        ax.set_xticks(xx,[LABELS[e] for e in data],fontsize=9);ax.grid(alpha=.2)
    axs[0].set(title="Mean peak-stage difference",ylabel="No wetland − wetland (mm)")
    axs[1].set(title="Weight assigned to exact-zero effects",ylabel="Normalized weighted mass",ylim=(-.02,.85))
    d=diagnostics.set_index("law").loc[list(LAW_LABELS)]
    short=["Baseline","Uniform","rho −0.25","rho +0.25"]
    axs[2].barh(np.arange(4),d.weight_ess_diagnostic,color=LAW_COLORS,height=.65)
    axs[2].set_yticks(np.arange(4),short);axs[2].invert_yaxis();axs[2].set_xlim(0,75)
    axs[2].set(title="Weight concentration",xlabel="ESS diagnostic (out of 64)")
    for i,r in enumerate(d.itertuples()): axs[2].text(r.weight_ess_diagnostic+1,i,f"{r.weight_ess_diagnostic:.1f}",va="center")
    axs[2].grid(axis="x",alpha=.2)
    handles,labels=axs[0].get_legend_handles_labels();fig.legend(handles,labels,ncol=2,loc="lower center",bbox_to_anchor=(.5,.08),frameon=False,fontsize=9)
    fig.suptitle("Exploratory parameter-law stress tests: reweight saved outputs on the same [0.5, 2] support",fontsize=13,y=.98)
    fig.text(.5,.02,"No new physical draws or simulations. Latent correlation concerns infiltration and common routing only; all laws are assumptions.\nESS is a diagnostic, not an interval sample size; no weighted confidence bands are asserted.",ha="center",fontsize=9)
    fig.tight_layout(rect=(0,.22,1,.91));savefig(fig,"06_parameter_law_checks.png")


def findings(data,atoms,checks_error,prefixes,subsets,outcomes,diagnostics,ranks,grid_summary,accounting):
    lines=["# Statistical robustness checks using saved results", "",
        "This analysis reads the existing five fixed-event, 64-pair result matrices and the existing September 2022 27-point parameter grid. It does not run an ODE, generate a new physical parameter vector, fit a surrogate, revise a parameter or alter an earlier result. All 320 saved pairs remain in the analysis: 320 valid, zero invalid and zero upstream safety violations. Zero or unfavorable response metrics are not excluded.", "",
        "## ECDF audit and zero effects", "",
        f"All 42,030 saved fixed-grid CDF rows match fresh count/N and Clopper–Pearson calculations. All 65 possible success counts for N=64 agree with SciPy's exact-binomial interface to a maximum absolute endpoint difference of {checks_error:.3g}. The original official rank/(N+1) arrays remain correct as plotting positions, but they are not the ties-aware empirical CDF. The statistical audit checks all saved outcome-summary metrics as well as the peak-stage differences.", "",
        "| Fixed event | Positive / zero / negative peak differences | Empirical median (mm) | Exact-zero fraction; pointwise 95% CP interval |",
        "|---|---:|---:|---|",
    ]
    for r in atoms.itertuples():
        lines.append(f"| {LABELS[r.event_id].replace(chr(10),' ')} | {r.positive_count}/{r.exact_zero_count}/{r.negative_count} | {r.empirical_median_m*1000:.4f} | {r.zero_mass:.4f}; [{r.zero_mass_cp95_lower:.4f}, {r.zero_mass_cp95_upper:.4f}] |")
    lines += ["", "Every exact-zero peak difference coincides with both models attaining the same full-window maximum at their initial river level. This does not mean the hydrographs or all other response metrics coincide. No near-zero tolerance was used to manufacture zero counts. Because no negative differences were observed, the empirical zero mass and empirical F(0) have identical counts and numerical intervals here; their population meanings remain different. Absence of an observed negative value does not establish zero population probability below zero. The pointwise F(0) intervals for 2019 [0.4482, 0.7006] and 2024 [0.3428, 0.5977] both include one-half. Their empirical medians are valid descriptions of these 64 draws, but the data do not precisely distinguish a nonpositive population median from a positive one under the assumed law.", "",
        "The CDF intervals describe Monte Carlo sampling uncertainty at predeclared fixed thresholds under an assumed parameter distribution. They are not the spread of physical outcomes, not simultaneous CDF bands, and not probabilities of historical weather occurrence. The five selected weather windows are fixed comparisons, not an iid weather sample. No 320-pair pooled risk bound is computed because events reuse the same 64 vectors and are not a random weather sample. See the [SciPy exact interval API](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats._result_classes.BinomTestResult.proportion_ci.html).", "",
        "`statistical_sampling_uncertainty.csv` separates the sample standard deviation and 5th–95th outcome quantiles from the estimated Monte Carlo standard error of each mean, s/√64. The latter is a precision diagnostic for iid draws under the assumed law, not a physical variability interval or a verified confidence interval for a tail quantile. Estimates across events are dependent because their vectors are matched.", "",
        "## Sample-size diagnostic", "",
        f"Original-order prefixes use 16, 32, 48 and 64 saved vectors. A separate diagnostic uses {SUBSET_REPETITIONS:,} fixed-seed random permutations (seed {SUBSET_SEED}), taking the same nested prefixes for every event. No replacement is used within a subset; the identical subset indices are shared across events. All permutations and replicate statistics are saved. The middle 95% ranges describe possible subsets of this specific finite set, not confidence intervals for an unknown population. At N=64 the range collapses mechanically; it does not establish convergence.", "",
        "| Fixed event | Mean ΔH, original prefixes 16 / 32 / 48 / 64 (mm) | 48-draw subset range for mean ΔH (mm) |",
        "|---|---:|---:|",
    ]
    for eid in data:
        p=prefixes[prefixes.event_id==eid].sort_values("n")
        s=subsets[(subsets.event_id==eid)&(subsets.n==48)&(subsets.metric=="mean_delta_peak_stage_m")].iloc[0]
        vals=" / ".join(f"{v*1000:.2f}" for v in p.mean_delta_peak_stage_m)
        lines.append(f"| {LABELS[eid].replace(chr(10),' ')} | {vals} | [{s.q025*1000:.2f}, {s.q975*1000:.2f}] |")
    lines += ["", "The subset tables also retain median and 95th-percentile peak differences, zero fractions, discharge differences, duration differences and total outlet-volume differences. These diagnostics can reveal dependence on which saved vectors are included, but cannot discover unvisited parameter-space behavior or verify a rare-event tail.", "",
        "## Parameter-law sensitivity without new simulations", "",
        "The original law is independent log-uniform on [0.5, 2] for three categories: upland infiltration, a single common routing-rate multiplier, and river resistance. Independent uniform multipliers on the same support are compared by self-normalized importance weights proportional to the product of the three multipliers. This changes the assumed input measure, not any physical parameter vector or response.", "",
        "Two modest descriptive dependence stresses retain the original log-uniform marginals and use latent-normal correlations −0.25 and +0.25 between infiltration and common routing; river resistance stays independent. The upland connection motivates examining dependence, but neither its sign nor magnitude is measured. Symmetric stress values were specified as hypothetical settings before reweighted outcome summaries, not estimated from the response. They do not represent new calibrated priors. No dependence involving river resistance is tested, so small changes in these two stresses cannot establish that independence assumptions are generally unimportant.", "",
        "For u_j=ln(m_j/0.5)/ln(4), z_j=Φ⁻¹(u_j), the copula density ratio is det(R)⁻¹/² exp[−zᵀ(R⁻¹−I)z/2]. The implementation independently checks this formula against the [SciPy multivariate-normal log-density](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.multivariate_normal.html). Full formulas, per-draw weights and marginal checks are saved in the statistical law tables.", "",
        "| Law | Weight ESS / 64 | Largest normalized weight | Mean raw density ratio |",
        "|---|---:|---:|---:|",
    ]
    for r in diagnostics.itertuples(): lines.append(f"| {LAW_LABELS[r.law]} | {r.weight_ess_diagnostic:.2f} | {r.maximum_normalized_weight:.4f} | {r.mean_raw_density_ratio:.4f} |")
    lines += ["", "These estimates reuse only 64 proposal draws. Self-normalization is finite-sample biased; a high ESS cannot establish adequate coverage of every output tail. ESS=1/Σw² measures weight concentration only; it is not substituted for N in binomial intervals. No weighted confidence bands are computed. The reweighting and its limitations follow [Owen, Monte Carlo theory, importance-sampling chapter, Sections 9.2–9.3](https://artowen.su.domains/mc/Ch-var-is.pdf).", "",
        "| Fixed event | Baseline mean ΔH (mm) | Uniform mean ΔH (mm) | Latent rho −0.25 / +0.25 mean ΔH (mm) | Uniform zero mass |",
        "|---|---:|---:|---:|---:|",
    ]
    for eid in data:
        p=outcomes[(outcomes.event_id==eid)&(outcomes.metric=="delta_peak_stage_m")].set_index("law")
        lines.append(f"| {LABELS[eid].replace(chr(10),' ')} | {p.loc['baseline_loguniform','self_normalized_weighted_mean']*1000:.3f} | {p.loc['uniform_multipliers','self_normalized_weighted_mean']*1000:.3f} | {p.loc['copula_negative','self_normalized_weighted_mean']*1000:.3f} / {p.loc['copula_positive','self_normalized_weighted_mean']*1000:.3f} | {p.loc['uniform_multipliers','weighted_zero_mass']:.4f} |")
    uniform_2019=outcomes[(outcomes.event_id=="historical_20190825T0200Z")&(outcomes.law=="uniform_multipliers")&(outcomes.metric=="delta_peak_stage_m")].iloc[0]
    lines += ["", f"The 2019 exact-zero mass changes from 0.5781 under equal weights to {uniform_2019.weighted_zero_mass:.4f} under the uniform-multiplier stress; its weighted inverse-CDF median becomes {uniform_2019.weighted_inverse_cdf_median*1000:.3f} mm. Thus the zero/positive median distinction is sensitive both to finite-sample uncertainty and to the uncalibrated parameter-law assumption. This is not evidence that one of the laws is more realistic.", "",
        "All reweighted peak differences stay nonnegative because every original saved peak difference is nonnegative and every weight is positive. That sign cannot be an independent validation of another parameter law or evidence about parameter combinations absent from the saved set. Output quantile ranges describe variability across hypothetical parameter values. Their estimation uncertainty is a different quantity. A weighted inverse-CDF median uses the first value reaching cumulative weight 0.5, which can differ from the interpolated even-N empirical median; both definitions are recorded explicitly.", "",
        "## Parameter patterns and limits of attribution", "",
        "Saved Spearman rank associations are independently reproduced, without significance testing or a new fit. They are marginal descriptions within each event, not causal effects or Sobol indices. The saved 27-point benchmark grid additionally supplies nine matched endpoint contrasts per factor: change the factor from 0.5 to 2 while holding the other two exactly fixed. No extra grid was evaluated for another event.", "",
        "| Factor | Outcome | Positive / zero / negative endpoint changes over nine contexts | Endpoint range |",
        "|---|---|---:|---:|",
    ]
    for r in grid_summary.itertuples():
        scale=1000 if r.metric=="delta_peak_stage_m" else 1
        unit="mm" if scale==1000 else "m³/s"
        lines.append(f"| {r.factor} | {'ΔH' if scale==1000 else 'ΔQ'} | {r.positive}/{r.zero}/{r.negative} | {r.minimum_endpoint_change*scale:.5f} to {r.maximum_endpoint_change*scale:.5f} {unit} |")
    lines += ["", "The grid shows why one rank correlation is insufficient: increasing infiltration raises the stage/discharge peak difference in eight contexts but lowers it in one; increasing common routing raises it in four and lowers it in five. River resistance raises ΔH in all nine contexts, while its ΔQ effect changes sign. These are controlled finite differences within the saved hypothetical benchmark model, not measured field sensitivities or general monotonicity results.", "",
        "The outlet rating is uncapped at these peaks, so ΔH_peak=R_r ΔQ_peak for a paired realization. A strong association of resistance with ΔH is partly the conversion from discharge to stage; it must not be read as an independent increase in hydraulic mitigation. Common routing remains a single category rather than three separately perturbed rates.", "",
        "All duration and volume signs remain available in the full outcome tables. Negative differences are normal modeled outcomes, not invalidity. In particular, all four additional storms have higher total outlet volume in the wetland arm. The wetland initial state contains 5,500 m³ more water, so the 144-hour volume contrast is not pure avoided floodwater. H0 is an initial-state reference, not a flood threshold. Initial states, event windows and the benchmark/non-benchmark weather-processing differences are unchanged.", "",
        "## What remains unresolved", "",
        "The saved data support a conditional robustness description, not calibrated distributions, a proven sample-size convergence result, or an environmental risk probability. A useful next approved experiment would append an independent batch of parameter draws under the original law, using exactly the same five fixed storms and paired model settings, and prospectively compare the mean, zero mass and tail quantiles. This would test whether the current estimates persist beyond reuse of the same finite 64-point set. It is a proposal only: no new draws or ODE evaluations were made here.", "",
        "## Figures and reproducibility", "",
        "- `plots/04_ecdf_zero_atoms.png`: exact jumps at zero and pointwise binomial intervals; whole-window maxima retained.",
        "- `plots/05_sample_size_stability.png`: matched finite-set subset diagnostics; ranges are not population confidence intervals.",
        "- `plots/06_parameter_law_checks.png`: reweighted means, zero mass and weight concentration; no weighted intervals.", "",
        "Run `python D/d2_analysis/robustness_checks/statistical_checks.py` from the repository root. Inputs are existing CSV/JSON files; outputs are confined to this directory. `tables/statistical_audit.json` records source hashes, library versions and statistical checks. The original results are never rewritten.", ""]
    (HERE/"statistical_findings.md").write_text("\n".join(lines))


def main():
    TABLES.mkdir(exist_ok=True);PLOTS.mkdir(exist_ok=True)
    data,draws,accounting=load()
    atoms,saved,checks_error=check_ecdf(data)
    prefixes,subsets=stability(data)
    outcomes,diagnostics=laws(data,draws)
    ranks,grid_summary=parameter_patterns(data)
    figures(data,atoms,saved,prefixes,subsets,outcomes,diagnostics)
    findings(data,atoms,checks_error,prefixes,subsets,outcomes,diagnostics,ranks,grid_summary,accounting)
    for name,digest in SOURCE_HASHES.items(): assert sha(D2/name)==digest, name
    audit={"all_source_hashes_unchanged":True,"source_sha256":SOURCE_HASHES,"script_sha256":sha(__file__),
        "libraries":{"numpy":np.__version__,"pandas":pd.__version__,"scipy":scipy.__version__,"matplotlib":matplotlib.__version__},
        "events":len(data),"existing_parameter_vectors":N,"saved_pairs":sum(len(x) for x in data.values()),
        "new_physical_parameter_vectors":0,"new_ode_evaluations":0,"statistical_subset_seed":SUBSET_SEED,
        "checked_fixed_cdf_rows":len(saved),"maximum_scipy_cp_endpoint_error":float(checks_error),
        "zero_atoms":atoms.to_dict("records"),"law_weight_diagnostics":diagnostics.to_dict("records"),
        "weighted_ci_computed":False,"sample_size_convergence_claimed":False}
    write_json(audit,"audit")
    print(json.dumps({"events":len(data),"pairs":320,"new_ode_evaluations":0,"cp_check_max_error":checks_error,
                      "source_hashes_unchanged":True,"output_directory":str(HERE)},indent=2))


if __name__=="__main__":main()
