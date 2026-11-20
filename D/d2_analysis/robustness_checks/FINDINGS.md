# Robustness checks: what the saved results establish

This is an analysis note, **not the final HW D.2 report**. No ODE, new physical parameter draw, weather generation, alternative initial condition, optimization, commit or push was performed. All five events and 320 paired outcomes remain included. Their model classifications remain complete, valid and safe under the existing checks; this task audits the saved evidence, not a new safety certification.

## Zero peaks are an initial-condition floor, not identical hydrographs

All **87 exact zero full-window peak differences** correspond to both models attaining their first maximum at the common initial river level, 97.35 m: nine in August 2023, 11 in July 2022, **37 in August 2019 and 30 in August 2024**. There are none in the September 2022 benchmark. Every zero pair has full saved arrays, and every pair's later river-stage waveform differs; the largest within-pair absolute differences range from **49.722 to 102.834 mm**. Peak-discharge differences are zero for the same pairs because of the shared monotone outlet rating.

These are genuine **zero values of the accepted full-horizon metric**. They do not show that wetland routing has no dynamic effect, and a secondary peak after rain begins must not be substituted to manufacture a positive primary benefit. The separate, weather-defined response window shows lower and delayed wetland peaks, followed by intervals when wetland-case river stage is higher. Its timing is explicitly descriptive.

All 320 peak/timing summaries were checked against 640 model diagnostics. Waveform checks cover 258 pairs: all 256 added-event pairs and two retained benchmark pairs. The remaining benchmark waveforms were not reconstructed. Actual paired gates agree to floating-point precision and have no latch disagreements, so these contrasts do not result from a changed policy. See [trajectory findings](trajectory_findings.md), [event/model peak magnitudes and times](tables/trajectory_peak_summary.csv), and [rainfall/event summary](tables/trajectory_event_summary.csv).

## Initial water: exact accounting, incomplete causal attribution

The wetland model starts with **518,100 m³** across all compartments, versus **512,600 m³** in the no-wetland model. All 320 pairs have recorded final states and complete block-integrated balance diagnostics, including the benchmark. Earlier summary-table blanks did not imply absent endpoint data; this check reads the original JSON endpoints without changing the old tables.

With Δ meaning no-wetland minus wetland, S total system storage, V river outlet volume, N net non-outlet external volume and r the recorded numerical residual:

`ΔV = −5,500 − ΔS_final + ΔN + Δr`.

N includes rainfall, external inflow, groundwater-boundary exchange and evapotranspiration; internal transfers are not counted as additional external water. Mean accounting terms across 64 draws per event are:

| Event | Final system stock, wet−no (m³) | Net non-outlet term, no−wet (m³) | Outlet volume, no−wet (m³) |
|---|---:|---:|---:|
| September 2022 | 6,218.46 | 772.21 | 1,490.67 |
| August 2023 | 4,218.18 | 758.11 | −523.71 |
| July 2022 | 4,047.40 | 784.48 | −668.12 |
| August 2019 | 3,889.92 | 701.13 | −908.95 |
| August 2024 | 3,980.68 | 777.93 | −741.40 |

The fixed −5,500 m³ initial term plus these terms and a mean residual of approximately +0.0005 m³ recover the outlet contrast. The largest paired arithmetic discrepancy is **1.24×10⁻¹⁰ m³**; this is consistency of recorded accounting, not an independent new numerical-convergence test.

The added events' wetlands themselves are almost empty at 144 h, while their whole systems still contain about 3,890–4,218 m³ more water. Most of that endpoint difference is in groundwater storage (mean differences about 3,866–4,192 m³), not water remaining in wetland compartments. Thus lower peaks, delayed recession, almost-empty wetland units and higher cumulative outlet volume can coexist. Neither final-storage differences nor net stock decline identify the fraction of original water retained: initial water and rainfall are not tagged by origin. The saved contrast cannot isolate the causal effect of exactly 5,500 m³ or justify subtracting 5,500 m³ from an output as a correction. See [storage findings and balance tables](storage_findings.md).

## ECDFs are correct; precision and prior dependence remain material

All **42,030** saved fixed-grid CDF rows and **630** summary values match independently recalculated results. Repeated zeros form one exact count/N jump. Continuous order-statistic plotting positions are not substituted for that jump. The existing pointwise Clopper–Pearson CDF intervals handle ties; they match SciPy's exact-binomial interface within **4.35×10⁻¹³**.

At zero, the 2019 interval is **[0.4482, 0.7006]** and the 2024 interval **[0.3428, 0.5977]**. Both span one-half, so the assumed-law median being zero versus positive is not precisely determined by 64 draws. These intervals quantify finite parameter-sampling uncertainty at fixed thresholds. They do not quantify all physical-model uncertainty, form simultaneous bands, or describe storm occurrence probabilities. No negative sampled effect does not prove an absence of negative effects elsewhere in parameter space.

Original-order mean peak reductions from 16 to 64 draws change **42.16→26.18 mm** for 2019 and **61.01→41.00 mm** for 2024. The 1,000 matched random-subset checks describe variation within the same finite set; their ranges necessarily collapse at 64 and cannot establish Monte Carlo convergence or reliable tails.

Exploratory reweighting to independent **uniform multipliers on the same [0.5, 2] bounds** changes the estimated means as follows. No physical response is recomputed.

| Event | Original log-uniform mean (mm) | Uniform reweighted mean (mm) |
|---|---:|---:|
| September 2022 | 91.573 | 107.881 |
| August 2023 | 91.045 | 117.158 |
| July 2022 | 73.814 | 97.922 |
| August 2019 | 26.178 | 38.840 |
| August 2024 | 41.003 | 59.616 |

For 2019, reweighted zero mass falls from 0.5781 to 0.4127 and the inverse-CDF median becomes 20.651 mm. Consequently, magnitude and median conclusions depend materially on the uncalibrated input law. Weight ESS is about 44/64, a concentration diagnostic, not the sample size for a new confidence interval. Self-normalized estimates remain approximate. Modest ±0.25 latent infiltration/routing correlation stresses change means less, but do not test dependence involving river resistance. Nonnegative reweighted effects follow arithmetically from the nonnegative saved responses; they are not independent validation of the new laws.

The existing 27-point benchmark grid reinforces the rating-equation limitation: increasing river resistance from 0.5 to 2 times baseline increases ΔH in all nine controlled contexts, but ΔQ increases in four and decreases in five. On the uncapped branch, `ΔH = R_r ΔQ`. Infiltration and routing effects are also context dependent; marginal rank associations are not causal sensitivity indices. See [statistical findings](statistical_findings.md) for definitions, uncertainty types, source references and all retained signs.

## Limited next experiments: proposals only

| Question | Proposed design before any execution | Scope and decision needed |
|---|---|---|
| Does antecedent state explain the initial-peak floor and persistence of extra initial water? | For August 2019 and August 2024, compare 24 h and 72 h of the same existing complete nonfreezing antecedent forcing. Use the three existing pilot vectors `mc_00016`, `mc_00024`, `mc_00018`; carry all resulting water states, gate positions and latch memory continuously into each unchanged 144 h event. | **12 new pairs / 24 model trajectories**, only after approval and an explicit opt-in state/memory handoff implementation. Start conditioning with the accepted states, including extra wetland water; delete or redistribute none. This tests antecedent conditioning, not a pure isolated 5,500 m³ intervention or proven equilibrium. Keep all failures/unsafe outcomes. |
| Do means, zero mass and upper quantiles persist beyond the same 64 points? | Append one separately recorded batch of **32** iid vectors using the existing official sampler, original log-uniform bounds and shared routing multiplier; proposed new seed **690321**. Apply each vector to both models under all five unchanged events. Keep the old 64 records and analyses intact. | **160 new pairs**, only after approval. Report the independent batch and combined 96-draw conditional estimates, MCSEs and pointwise CDF intervals separately. No automatic further batch, distribution tuning or convergence claim; new sampling cannot calibrate the assumed law. |

The [antecedent-data audit](tables/antecedent_continuous_coverage.csv) found 351, 93, 0, 254 and 109 consecutive supported hours before the benchmark, August 2023, July 2022, August 2019 and August 2024 windows, respectively. Both 24 h and 72 h proposals are supported for the two zero-heavy events. July 2022 lacks the immediately preceding rainfall hour. No missing input is filled, and complete station-derived/interpolated weather plus gridded solar is not described as direct observation of every variable. Conditioning durations are numerical experiment choices, not established equilibration times.

The first proposal addresses the main interpretation gap more directly; the second addresses sampling precision only. Neither has been executed. Within the saved tests, nonnegative primary peak differences and delayed routing are supported. Their magnitude, median and outlet-volume interpretation depend on initial conditions and parameter assumptions. Universal mitigation, causal attribution to initial water, calibrated priors, population-tail precision and flood-frequency claims remain unsupported.
