# Rainfall coverage and the next conditional high-rainfall experiment

## Decision and scope

**The existing 16 non-freezing pairs do not provide a representative heavy-rainfall uncertainty experiment for wetland flood mitigation.** They are a successful hydraulic/numerical pilot under ordinary seasonal weather, with one substantially wetter member. This is a statement about the delivered sample and its sampling design, not a claim that every dry scenario is defective or that the weather generator has been validated for rare storms.

This note analyzes existing cached weather and saved results only. No weather was generated, no new observations were downloaded, and no ODE was run for this assessment. The proposed experiment below has **not** been implemented or approved. No downstream flood threshold is introduced.

## What the current weather sample represents

The [recovered configuration](../g5_diagnostics/recovered/original_config.json) specifies ordinary September six-day weather windows, without storm/flood/outcome selection. It retains ten Copula/VAR and ten block-bootstrap IDs. Four bootstrap IDs contain freezing temperatures and remain explicit exclusions; the conditional hydrologic pilot has ten Copula/VAR and six bootstrap members. These generator groups and the cold exclusion must not be hidden by a pooled flood-probability claim. Weather-generator burn-in is not hydrologic spin-up.

Depths below come from the unchanged hourly scenario files. The maximum 1-hour depth is the largest supplied interval-mean rate times one hour; 6/24-hour values are sums over consecutive complete bins with integer-hour starts. They are not subhourly intensity measurements.

| Existing non-freezing sample | Six-day depth, mm | Maximum 1-hour depth, mm | Maximum 6-hour depth, mm | Maximum 24-hour depth, mm |
|---|---:|---:|---:|---:|
| Median over 16 IDs | 13.222599 | 4.098124 | 5.410917 | 6.764320 |
| Largest value over 16 IDs | 80.760000 | 15.700000 | 46.320000 | 77.960000 |
| Largest value among the other 15 IDs, excluding `block_bootstrap_00000` | 22.075058 | 6.818240 | 12.800000 | 12.800000 |

`block_bootstrap_00000` supplies all four maxima in the second row. One warm member has no quantified rain, and `block_bootstrap_00008` has only 0.04 mm. The accepted deterministic September 2022 forcing totals **145.3 mm**. Thus even the wettest pilot member has about 55.6% of that event's total; the other 15 have at most about 15.2%. Rain totals alone do not characterize a storm, which is why duration-specific depths are reported separately. The full [weather summary](results/weather_summary.csv) also retains the four cold IDs and the deterministic event, with its distinct 120 s forcing reconstruction identified.

In the saved 15 s paired simulations, 15/16 full-window river maxima occur at the prescribed initial 97.35 m. The outlet is a monotone function of river stage under the fixed boundary settings, so substituting peak discharge does not remove this initial-condition dominance. These zero peak differences do not demonstrate zero wetland routing effects. See [metric definitions](metric_definitions.md) for separately labeled response diagnostics that preserve the full-window maxima.

## Historical reference and availability

The local source is `C/data/processed/weather_hourly.csv`, SHA-256 `7edf24675b891d680863079323efb5da46e3d0c89a0640ed5852378a4be58fc4`. It combines KMKE rain, KMWC temperature/wind and POWER solar estimates; it is a regional proxy, not measured areal weather over the hypothetical installation. Its [source/QC report](../../C/docs/data_sources.md) documents overlap-weighted hourly rain, trace lower bounds, interpolation, and missingness. This is the original standardized weather dataset, not the separately aligned FIT dataset. Redistribution can alter clock-hour maxima; it does not reveal unobserved subhour timing.

The [historical reference table](results/historical_weather_reference.csv) reports complete rolling windows wholly inside September UTC. Missing rain is never replaced by zero. A joint non-freezing window additionally requires all four variables and T > 0 C throughout. Rolling windows overlap; their counts are **not independent storm counts**.

| Year | Largest complete 24-hour rain depth, mm | Rain-complete 144-hour windows | Joint-complete non-freezing 144-hour windows |
|---|---:|---:|---:|
| 2019 | 38.100000 | 364 | 364 |
| 2020 | 5.526667 | 241 | 241 |
| 2021 | 29.000000 | 432 | 432 |
| 2022 | 141.400000 | 196 | 196 |
| 2023 | 26.000000 | 239 | 0 |
| 2024, held-out context | 25.700000 | 427 | 326 |

The 2022 standardized September maximum 6-hour depth is 58.493333 mm, and maximum 144-hour depth is 145.3 mm. These clock-hour statistics characterize the same regional storm period, but they are not replacements for the separately prepared B forcing. The absence of a complete joint 144-hour September 2023 window is a data-availability limitation, not absence of rain.

For 2019–2023 September alone, 3,027 of 3,485 possible rolling 24-hour windows have complete rain; 778 have positive quantified depth. Their empirical 95th percentile is **37.806000003 mm**. For the existing bootstrap's adjacent-month scope, August–October 2019–2023, the corresponding counts are **9,229/10,925**, with **2,535 positive** windows and empirical 95th percentile **35.1379999994 mm**. Quantiles use linear interpolation over the positive complete rolling-depth values. These are descriptive short-record percentiles, not annual exceedance probabilities or an official definition of extreme precipitation. The 2024 observations do not determine either threshold.

Reproduction of these reference calculations requires only the cached CSV: group by UTC year and selected month set; compute `P.rolling(h, min_periods=h).sum()` for h = 6, 24 and 144; count a joint window only if the rolling sum of the complete-four-variable/non-freezing indicator equals h. Calculate the stated pooled threshold with `np.quantile(positive_complete_24h_depths, 0.95, method='linear')`. No windows cross the selected season/year boundaries.

## One recommended next design: an intact joint historical high-rainfall catalogue

Use a **rainfall-conditioned historical-event stress test**, rather than continuing to sample ordinary six-day windows and hoping that enough severe rain appears. The following is one concrete proposed design, requiring approval before extraction or simulation:

1. **Selection population and criterion.** Use cached August–October 2019–2023 records, matching the existing adjacent-month seasonal pool. Select 24-hour rain windows at or above the fixed empirical upper-tail criterion **35.1379999994 mm** above. This is a proposed sample-selection rule, not a flood threshold, a NOAA return-period depth, or an annual-frequency estimate. Record 1/6/24-hour depths and six-day totals to show which storm shapes are represented. Do not change the criterion after seeing hydraulic outcomes.
2. **One representative per rain episode.** Before hydrologic runs, group qualifying window starts separated by less than 72 hours into consecutive clusters; retain the greatest 24-hour depth in each cluster, using the earliest start for ties. The 72-hour rule is a declared declustering assumption, not proof of meteorological independence. Report cluster membership and any overlapping context between retained windows. Do not count repeated/overlapping windows as independent extreme events.
3. **Preserve the joint sequence.** For each representative 24-hour interval starting at s, propose the intact 144-hour source window `[s - 24 h, s + 120 h)`. Carry P, T, wind and solar together, retaining source UTC dates, original ordering, source flags, native hourly supports and existing spatial-transfer assumptions. No rain-only multiplier, shuffled covariates, stitched storm blocks, shifted solar clock or fabricated tails. This catalogued variation samples observed storm timing, duration and associated weather within the available record; it cannot generate unprecedented extremes.
4. **Eligibility and denominators.** Select episodes using rain first, then report unavailable four-variable context, freezing context and record-edge exclusions separately. Do not silently replace an excluded episode with a less severe one. Retain source gaps and the original liquid-water domain. If too few complete distinct episodes remain, report inadequate support; do not lower the cutoff, shorten the horizon or claim a large uncertainty sample by duplicating events. Freeze the design before inspecting 2024, which can serve only as a separately reported chronological check.
5. **Matched physical cases and declared initialization.** Each eligible sequence drives both physical models with the same source timestamps, accepted shared initial states, fixed design, boundary assumptions, controller policy and numerical criteria. Retain the additional 5,500 m³ initial wetland water and verify realized gate histories. The 24 hours of preceding context is part of a hypothetical event simulation beginning from the accepted assumed state; it is **not** an observed initialization or demonstrated spin-up. Retain the full 144-hour metrics, and identify any storm-response subwindow in advance without replacing those metrics.
6. **Report conditional outcomes, including unfavorable ones.** Report each event's physical hydrographs, peaks, timing, volume/storage accounting, applicability, mass-balance acceptance and safety margins. Keep negative wetland effects and complete unsafe trajectories. Incomplete/invalid cases have no invented full-event peaks. This is a finite conditional stress-test catalogue, not a joint flood-return-period estimate or an unconditional flood-benefit probability. No downstream flood duration or exceedance probability is available without a separately justified threshold.

The resulting study would address the observed high-rainfall tail more directly while preserving observed multivariate sequences. It would remain conditional on a short, incomplete regional record and hypothetical initial conditions, geometry and operation. It should not be advertised as a complete rare-event risk analysis.

## Why rainfall frequency must not become a flood-return-period label

NOAA distinguishes annual exceedance probability from event recurrence and also distinguishes clock-constrained observations from unconstrained maxima. Consequently, the overlapping-window percentile used above is not a NOAA precipitation-frequency estimate. [NOAA HDSC glossary](https://www.weather.gov/owp/hdsc_glossary).

If a later study needs frequency-rated storms, the applicable NOAA documentation lists Wisconsin in Atlas 14 Volume 8. No local Atlas 14 depth was retrieved or inserted here. A point depth would require the chosen location/duration, confidence limits and an explicit spatial transfer or point-to-area treatment for the 0.66 km² installation; it must not be replaced by the roughly 46.879 km² USGS catchment. HEC also requires temporal-pattern and storm-area choices, which cannot be inferred from a single 24-hour depth. [NOAA precipitation-frequency documents](https://www.weather.gov/owp/hdsc_currentpf), [USACE hypothetical-storm technical reference](https://www.hec.usace.army.mil/confluence/hmsdocs/hmstrm/meteorology/precipitation/hypothetical-storm).

USACE explicitly warns that assigning a precipitation frequency does not establish the same frequency for simulated flood runoff; antecedent conditions, losses, temporal/spatial patterns and regulation also matter. Therefore neither a future Atlas 14 point 24-hour AEP nor this catalogue's selection frequency would be the probability of an entire joint weather trajectory, the wetland benefit probability, or a river flood return period. [USACE precipitation-frequency and flow-frequency guidance](https://www.hec.usace.army.mil/confluence/sspdocs/ssptutorialsguides/combining-precipitation-frequency-and-flow-frequency-results/bulletin-17c-flow-frequency-augmented-by-precipitation-frequency-with-hec-hms).

Primary guidance was consulted on 2026-11-18. Only reference pages were read; no new meteorological dataset was acquired.
