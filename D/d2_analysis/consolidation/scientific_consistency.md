# Cross-study scientific consistency

This review reconciles the saved Stage D studies without running the physical models, fitting a model, drawing parameters, generating weather, or changing previous results. The portable [audit script](cross_study_audit.py) independently recomputes primary CSV contrasts, published numerical summaries, tied empirical CDF counts, pointwise confidence limits, storage identities and antecedent-window summaries. It reads saved JSON records for study plans and prior independent audits. Its **3,355 checks pass** across 32 input files; the input collection is unchanged during the check. Detailed full-array, controller and solver-stage checks remain documented in the original independent audits, rather than being claimed as newly repeated here.

## Study populations and reuse

[study_ledger.csv](study_ledger.csv) distinguishes primary physical configurations from copied results, numerical repeats and saved-result reanalysis.

| Study | Correct analysis population | Reuse and interpretation |
|---|---:|---|
| Original weather pilot | 20 generated IDs; 16 complete warm-weather pairs | Four cold scenarios fail the liquid-water applicability preflight. The 7.5 s repeat is a numerical check of the same 16 configurations. The deterministic controls are separate. |
| Fixed-storm factorial | 27 parameter combinations | A deterministic sensitivity design, with no assigned probabilities. Two pilot corners and the existing deterministic center were reused. |
| Fixed-storm Monte Carlo | 64 parameter draws | One exploratory iid sample under the declared priors. The first two pilot results were reused. |
| Five selected historical storms | 320 paired evaluations of the same 64 vectors | The September 2022 benchmark is an exact copy of the preceding 64-case Monte Carlo results. The four additional storms contribute 256 additional parameter/event configurations; their 12 pilot configurations are included, not counted again. |
| Robustness and importance-weight checks | The existing 320 paired results | No new model run or parameter draw. Endpoint accounting covers all 320 pairs; retained full waveforms cover 258. |
| Antecedent conditioning | 12 new physical configurations | Two events × two prefix lengths × three existing parameter vectors. Six unconditioned comparisons are reused. One of the 12 pairs is also repeated at 7.5 s: 24 primary model runs plus two refinement runs. |

The reviewed union contains **375 primary physical configurations**: 16 + 27 + 64 + 256 + 12. This is neither an iid sample size nor an inventory of every historical ODE execution. The factorial and Monte Carlo plans contain 91 distinct physical parameter vectors, but only 64 are iid draws. No proposed additional 32-vector batch has been executed in these results. Reweighting the existing 64 draws does not create another ensemble.

## Consistent signs and downstream metrics

Paired downstream-response contrasts use **no-wetland minus wetland**. Extra-storage figures explicitly use the reverse, wetland minus no-wetland, so positive values denote additional stored water. A positive river-peak contrast means that the wetland case has the lower maximum. It does not imply that wetland discharge is lower at every time, that discharge duration is shorter, or that its finite-window outlet volume is smaller.

The initial river stage is 97.35 m. Fifteen of the 16 warm-weather pilot pairs reach both maxima at the initial state. In the five historical-event ensembles, exact zero peak contrasts occur in 0, 9, 11, 37 and 30 of the 64 pairs, respectively: **87 total**. Every zero coincides with both maxima occurring at initialization. The prior retained-waveform audit found later nonidentical responses in all 87 cases. A zero maximum contrast is not evidence that the two hydrographs coincide.

The unchanged outlet rating is

\[
Q_{out}=\min\!\left(3,\max\!\left[\frac{H_r-97.15}{R_r},0\right]\right).
\]

All compared maxima are below the outlet cap. Within a matched pair, the same resistance therefore gives \(\Delta H_{peak}=R_r\Delta Q_{peak}\). The stage and discharge peak comparisons are linked by this rating, rather than being independent evidence of mitigation.

There is **no accepted downstream flood threshold**. Duration above the initial stage, discharge above the initial discharge, and maximum rise are explicitly response benchmarks. They cannot be labeled flood exceedance probabilities, flood durations or return periods.

## Results that must remain distinct

| Event | Mean peak reduction (mm) | Median peak reduction (mm) | Zero contrasts / 64 | Mean outlet-volume contrast (m³) |
|---|---:|---:|---:|---:|
| September 2022 benchmark | 91.573 | 94.055 | 0 | +1490.669 |
| August 2023 | 91.045 | 86.859 | 9 | −523.707 |
| July 2022 | 73.814 | 69.598 | 11 | −668.121 |
| August 2019 | 26.178 | 0.000 | 37 | −908.954 |
| August 2024 | 41.003 | 7.341 | 30 | −741.396 |

All 320 unconditioned peak contrasts are nonnegative. All 256 additional-event outlet-volume contrasts are negative. The two statements are compatible: maxima and accumulated volume measure different aspects of a finite response. Wetland duration above the common original initial stage is longer in 27, 8, 4, 1 and 1 cases, respectively. These adverse or mixed outcomes remain in the results.

## Initial water and exact accounting

The original wetland model begins with 3,000 + 2,500 = **5,500 m³** in the wetland compartments. These compartments are absent from the no-wetland model; shared physical initial storages are unchanged. The paired total initial stock contrast is therefore −5,500 m³.

For each completed simulation, let \(r=S_T-S_0-E\), where \(E\) is the saved integral of all net external water fluxes, including outlet loss. Let \(N=E+V_{out}\), the net external contribution excluding the outlet. Then

\[
\Delta V_{out}=\Delta S_0-\Delta S_T+\Delta N+\Delta r.
\]

Internal gate, wetland, soil and groundwater transfers are not added as independent external inputs. Recomputing the saved decomposition closes to less than \(3\times10^{-10}\) m³. All 64 benchmark endpoint pairs are directly recorded, despite the absence of 62 benchmark full waveforms.

Across the five events, the mean extra final whole-system wetland-case storage is approximately 6,218, 4,218, 4,047, 3,890 and 3,981 m³. Most is in groundwater; mean remaining water in the wetland compartments is only about 13.19 m³ for the benchmark and below 0.002 m³ for each additional event. Whole-system retention must not be described as water still stored in the wetland units. Neither stock differences nor their ratio to 5,500 identify how much original wetland water remains: no tracer separates initial water, rainfall and other exchanges. The saved net external term does not identify a separate ET effect.

## Antecedent conditioning and window definitions

The conditioning experiment carries every physical storage and the actual gate/controller memory continuously across the event boundary. It does not empty the wetlands, redistribute water or restart the controller. The 144 h event window and the full 168 h or 216 h conditioned horizon answer different questions. Their water volumes add by phase; their maxima combine by taking the larger phase maximum, not by addition.

For the three selected vectors, the 24 h and 72 h prefix cases give mean **event-window** peak reductions of 71.970 and 69.827 mm in August 2019, and 92.763 and 90.285 mm in August 2024. All twelve conditioned event-window contrasts are positive. The corresponding **full-horizon** means are only 38.869 and 29.407 mm, and 63.662 and 56.920 mm; each three-vector group still contains an exact zero. Both results are retained. Small negative peak contrasts during some antecedent phases are also retained.

The event-window outlet contrasts become positive while all twelve full-horizon outlet contrasts remain negative. This is a temporal accounting change with a different event-start condition, not proof that conditioning removes an initial-water bias or estimates a pure wetland-storage effect. Carried groundwater and river storage differ between the two model topologies. Controller memory was allowed to evolve continuously; in these particular pairs the actual gate histories agree to roundoff and there are no latch mismatches. Duration measured above each model's own window-initial stage is not a common flood-threshold comparison. Three selected parameter vectors do not estimate a probability distribution, and the two prefix lengths do not establish physical equilibration.

## Statistical and numerical limits

The three multiplier groups use explicitly assumed, independent log-uniform distributions on [0.5, 2]: infiltration capacity, one common routing-rate factor and river drainage resistance. The three individual routing rates have perfect relative dependence through that common factor. These are sensitivity/uncertainty assumptions for a hypothetical installation, not measured Milwaukee parameter bounds or an empirically fitted joint distribution.

The five historical events reuse the same parameter vectors and were selected as fixed storm tests. They are not independent weather draws. The benchmark and additional events retain their documented, different weather-preparation conventions. No annual event probability or joint flood return period follows from these experiments.

The audit recomputes all 42,030 saved fixed-grid empirical CDF rows, counting ties with \(\#\{Y_i\le y\}/64\). Their Clopper–Pearson limits are pointwise intervals under the assumed iid parameter sampling, not simultaneous bands. The official rank-based plotting positions remain separately identified. The independent-uniform alternative is an importance reweighting of the same 64 draws: its effective sample size is about 43.9935 and largest normalized weight about 0.0497851. ESS measures weight concentration; it is neither a new sample count nor a replacement binomial denominator. Resampling subsets of the same finite sample does not demonstrate convergence of an enlarged Monte Carlo run.

All primary outcomes in these completed cohorts are numerically valid and satisfy their saved upstream safety checks. This does not justify silently discarding a future unsafe or failed sample. Safety evidence covers saved output nodes and successful solver-stage evaluations; it is not continuous-time certification. The unchanged relative mass tolerance is \(10^{-4}\). Specific 15 s versus 7.5 s comparisons passed with the controller period held at 120 s. Earlier coarse-resolution failed balance checks remain preserved. These checks do not prove convergence over all possible parameters, weather or hydraulic regimes.

The accepted feedback trajectories exercise closed and energy-gate flow. Their success does not validate the unused low-head transition law, freezing/snow processes, real facility dimensions, empirical calibration, or local engineering safety. The studies support conditional comparisons of the scoped hypothetical model.
