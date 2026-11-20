# Saved-result storage and initial-water accounting

This check reads existing results only. No ODE, new initial-condition experiment, parameter sampling, forcing adjustment or optimization was performed. All primary comparisons retain the accepted extra **5,500 m³** in the wetland compartments and unchanged shared initial states. The analysis does not identify the trajectory of those initial water molecules.

## Coverage and source distinction

All **320 paired endpoint balances** are available. In particular, all 64 benchmark pairs retain `initial_storage_m3`, `last_completed_state_m3` and complete `completed_blocks_mass_diagnostics` in their diagnostic JSONs. The earlier historical summary CSV left final storage blank; that omission does not mean the underlying endpoint data are absent. This new analysis reads those recorded endpoints directly and does not rewrite the previous tables.

Full saved waveforms are available for **258 pairs**: all 256 additional-event pairs and benchmark `mc_00000`/`mc_00001` from the original parameter pilot. Benchmark timing summaries based on two predetermined retained draws are descriptive examples, not estimates representative of all 64 draws. No other benchmark waveform is reconstructed.

## Paired balance identity

For each case, let S0 and ST denote the sum of all water-storage states at the beginning and end; V the integrated river outlet discharge; N the net non-outlet external flux integral; and r the saved numerical balance residual. The unchanged model identity is

```
ST - S0 = N - V + r
N = integral(P*A_rain + Q_ext + Q_gx - E_total) dt
  = saved_external_net + V

Delta V = Delta S0 - Delta ST + Delta N + Delta r
Delta = no-wetland minus wetland
Delta S0 = -5,500 m3
```

The positive sign of Delta r follows from the stored definition `r=ST-S0-external_net`. Internal gate releases, wetland drainage/overflow, soil transfers and other transfers between compartments are not added to external flow. N is recovered from already recorded total external-net accounting and outlet volume; it is not an independent flux measurement or a second conservation test. Its paired contrast combines groundwater-boundary exchange and evapotranspiration differences after common rain/external inflow cancel. These components cannot be separated from the compact arrays used here, so no ET-only attribution is made.

Mean values across 64 pairs per fixed event are:

| Event | Initial term (m³) | Final term, wet−no (m³) | Net non-outlet term (m³) | Residual term (m³) | Outlet no−wet (m³) |
|---|---:|---:|---:|---:|---:|
| Sep 2022 * | -5,500.00 | +6,218.46 | +772.21 | +0.000495 | +1,490.67 |
| Aug 2023 | -5,500.00 | +4,218.18 | +758.11 | +0.000500 | -523.71 |
| Jul 2022 | -5,500.00 | +4,047.40 | +784.48 | +0.000499 | -668.12 |
| Aug 2019 | -5,500.00 | +3,889.92 | +701.13 | +0.000499 | -908.95 |
| Aug 2024 | -5,500.00 | +3,980.68 | +777.93 | +0.000500 | -741.40 |

The largest absolute paired-identity discrepancy is **1.24e-10 m³**. The largest case endpoint-balance discrepancy is **1.16e-10 m³**. Recorded outlet integrals match available saved discharge arrays within **0 m³**. These small arithmetic discrepancies check consistency of saved quantities; numerical timestep validation remains the separately documented original pilot.

The explicit initial-stock term is the same −5,500 m³ in every pair. Final storage and net external fluxes respond dynamically and differ by case. Therefore the four-new-event negative outlet-volume differences cannot be assigned entirely to release of the initial 5,500 m³, nor can the benchmark positive difference be called permanent storm-water removal. No exact causal percentage of a volume contrast is identified by rearranging the balance.

## Where the final system-storage contrast resides

The complete saved endpoint state vectors support a compartmentwise comparison. Mean final stocks below use wetland minus no-wetland; the no-wetland model has no wetland compartments, so their absent contribution to the total sum is zero, without inventing simulated zero-valued states.

| Event | Groundwater difference (m³) | Other shared-state difference (m³) | Wetland-unit stock (m³) | Total system difference (m³) |
|---|---:|---:|---:|---:|
| Sep 2022 * | 6,167.15 | 38.11 | 13.1944 | 6,218.46 |
| Aug 2023 | 4,192.22 | 25.97 | 0.00105174 | 4,218.18 |
| Jul 2022 | 4,022.47 | 24.93 | 1.89813e-05 | 4,047.40 |
| Aug 2019 | 3,865.95 | 23.97 | 3.64897e-05 | 3,889.92 |
| Aug 2024 | 3,956.16 | 24.51 | 0.000185837 | 3,980.68 |

The final total-system contrast is concentrated in groundwater, while the additional-event wetland compartments are numerically near empty. Final system storage therefore must not be described as water still held in the wetland units. This is a state-location accounting result, not tagged-water tracking: the final groundwater contrast divided by 5,500 is **not** a fraction of original wetland water remaining. Rainfall, seepage, external exchange and evapotranspiration act together during the event. The full [compartment table](tables/storage_compartment_differences.csv) preserves every pair, shared state and absent-compartment flag; its [summary](tables/storage_compartment_summary.csv) retains ranges.

## Wetland stock and late recession

The last positive-rain interval is defined by the saved forcing, not an outcome-selected hydrograph threshold. The new-event hourly rain is held through the end of its last positive hour; the benchmark retains its original linear 120 s forcing, including the final positive-to-zero segment. Stock is interpolated on the saved 15 s state grid. Net stock decline after rain ends is `S_w(t_rain_end)-S_w(T)`; it is not gross wetland drainage, because gate inflow, drainage, overflow, seepage and evaporation can continue together.

| Event | Timing pairs | Mean endpoint wetland stock across 64 (m³) | Mean post-rain net stock decline, timing subset (m³) | Timing pairs with more post-rain wetland-case outlet |
|---|---:|---:|---:|---:|
| Sep 2022 * | 2 | 13.1944 | 1,216.92 | 2/2 |
| Aug 2023 | 64 | 0.00105174 | 100.34 | 64/64 |
| Jul 2022 | 64 | 1.89813e-05 | 18.87 | 64/64 |
| Aug 2019 | 64 | 3.64897e-05 | 530.23 | 64/64 |
| Aug 2024 | 64 | 0.000185837 | 4,956.57 | 49/64 |

`storage_wetland_timing.csv` also retains the combined and branch storage peaks/times and, where achieved, the final downward crossing into combined stock≤55 m³. The 55 m³ mark is a declared descriptive threshold equal to 1% of the primary initial wetland stock. It is neither a safety threshold nor evidence that 99% of the tagged initial water has left. Very small endpoint volumes reflect the numerical model's near-empty asymptote, not measured water volumes at the displayed precision. The terminal storage and post-rain volume patterns describe late recession; they do not reconstruct unrecorded normal-drainage or overflow hydrographs.

## Proposed alternative initial-condition experiment — not executed

A defensible alternative is **common antecedent-weather conditioning**, preserving water through the unchanged equations. The existing [coverage audit](tables/antecedent_continuous_coverage.csv) finds 351, 93, 0, 254 and 109 consecutive supported nonfreezing hours immediately before the benchmark, August 2023, July 2022, August 2019 and August 2024 event starts, respectively. July 2022 has a missing rainfall hour immediately before its start; missing rain must not be treated as zero. The usable record combines station-derived/interpolated weather, grid-estimated solar and retained trace/QC flags; it is not uniformly direct observation.

A small specific proposed check is **24 versus 72 hours of common antecedent conditioning for August 2019 and August 2024**, using only the existing three pilot vectors `mc_00016`, `mc_00024` and `mc_00018`. This would require 12 pairs / 24 model trajectories including conditioning and the original event, with no new parameter draw. Both durations are supported by the existing [window audit](tables/antecedent_windows.csv). They are exploratory durations, not assumed equilibration times. This proposal is not executed, and it does not silently extend to a 7- or 14-day antecedent record. The primary comparisons remain unchanged.

At the start of that antecedent window, retain the accepted physical storage values for each model, including the wetland model's 5,500 m³ extra water. Apply the same previously saved parameter vector, hypothetical geometry, boundary conditions and feedback policy to both models throughout the common antecedent weather. Carry every resulting water state and the actual gate positions/latch memory continuously into the original 144-hour event; do not reset river level, soil storage, groundwater, wetland water or actuators at the event boundary. No water is deleted, moved arbitrarily between compartments or inserted to equalize totals. The 120 s controller clock must remain aligned through the boundary, and the event forcing remains the saved unchanged input.

This would produce explicit, model-conditioned alternative event-start states. Initial values at the antecedent-window start would still be assumptions, and finite conditioning is not evidence of equilibrium, spin-up convergence or observed antecedent reconstruction. All initial applicability checks, water balances, safety outcomes and numerical criteria remain unchanged. Invalid or unsafe conditioning/event runs must remain reported. Available longer durations could support a separately approved initial-memory sensitivity check; no duration or feasibility outcome is assumed here.

The two models can develop different shared event-start storages and actual gate positions under the same antecedent weather. This experiment therefore tests sensitivity to antecedent state and the persistence of the original initial-water difference; it does **not** isolate a pure causal effect of exactly 5,500 m³ or create matched event-start states by construction. Event metrics must include each resulting event-start stage and retain any initial peak. No such conditioning or alternative trajectory is executed by this saved-result analysis.

## Files

- [System balances](tables/storage_model_balances.csv): all 640 model endpoints and source checksums; full-waveform availability is explicit.
- [Paired decomposition](tables/storage_paired_decomposition.csv): all 320 outlet-volume identities.
- [Wetland timing](tables/storage_wetland_timing.csv): full-waveform diagnostics, with absent benchmark timings left blank.
- [Event summaries](tables/storage_event_summary.csv), [audit](tables/storage_audit.json), and [figure](plots/03_storage_balance.png).

Figure 03 shows the mean accounting terms, final system-storage contrasts and wetland-compartment endpoint stocks. It is an accounting illustration, not causal source-water tracing or empirical engineering validation.
