# Figure and table review

Reviewed on 2026-11-20. All 30 existing PNGs were opened and visually inspected, alongside their captions and associated saved tables: four original-weather figures, nine parameter-study figures (including seven official diagnostics), five historical-event figures, six saved-result robustness figures and six antecedent-condition figures. The [inventory](figure_inventory.csv) records dimensions and file hashes. No existing figure, table, model, forcing or result was overwritten, and no ODE or fitting run was started for this review.

No material discrepancy requiring a figure correction was found. Axes, units, comparison signs, legends and labeled time windows are consistent with the saved analyses. Dense multi-panel figures should be opened at full size; the index preserves links to the existing originals rather than producing smaller duplicate images.

| Numerical cross-check | Result |
|---|---|
| Original-weather primary ECDFs | At both 15 s and 7.5 s: 15/16 zero differences among warm pairs; 10/10 Copula/VAR and 5/6 bootstrap. The separate deterministic pair is excluded from these denominators. |
| Six official parameter marginal plots | Each uses the same 64 saved draws. Sorted values agree with the summary table within 2.86e-14 in displayed physical units; plotting positions agree with rank/65 within 1.12e-16. No ties occur in these six distributions. |
| Historical fixed-grid CDF table | Recomputed all 42,030 counts from the saved paired results: zero count mismatches and zero count/N discrepancies. Independently recomputed 95% Clopper–Pearson endpoints agree within 1.12e-16. |
| Historical zero-primary counts | September 2022: 0/64; August 2023: 9/64; July 2022: 11/64; August 2019: 37/64; August 2024: 30/64. No negative full-window peak difference was sampled. |
| Zero versus unavailable waveform | All 87 zero-primary cases have both full-window maxima at time zero and distinct subsequent hydrographs in the saved waveform audit. Full waveform coverage is 258 pairs, while all 320 have primary summaries/endpoints. Missing benchmark waveforms are not treated as zeros. |
| Antecedent accounting | Twelve new paired configurations and six reused no-prefix pairs; window-specific rows and numerical refinements are repeated summaries, not extra parameter draws. |

The numerical checks above read existing CSVs and evaluated counts/sorted values/binomial formulas. They do not constitute a new dynamic-model validation or a proof of numerical convergence.

## Interpretation checks

- Original weather ECDFs use tie-aware count/N. Official parameter plots retain rank/(N+1) and their pointwise order-statistic convention. Historical main plots distinguish the solid count/N CDF from dotted official rank positions; the shaded intervals come from fixed-threshold exact-binomial counts. A zero atom must not be replaced by the last tied rank/(N+1).
- August 2019 F(0)=37/64=0.578125 and August 2024 F(0)=30/64=0.46875. In this sample these equal the exact-zero fractions because no negative value occurs. Their population definitions remain different. Pointwise limits are not simultaneous confidence bands or estimates of historical storm frequency.
- Full-window maxima retain the initial state. Secondary wet-spell response and running-minimum drawup are explicitly separate; a nonzero secondary response does not turn a zero primary peak difference into a positive primary result.
- River stage is an elevation in metres; the comparison figures sometimes convert differences to millimetres. Outlet discharge is in m³/s, integrated volume in m³, rainfall intensity in mm/h, mapped hourly rainfall depth in mm, solar irradiance in W/m², temperature in °C and wind in m/s. The fixed rainfall-window totals are not presented as recurrence probabilities.
- Peak stage and peak discharge are connected by the uncapped outlet law in these retained results. Their agreement is not two independent pieces of validation. The coefficient of variation printed for absolute elevation in official diagnostics is datum-dependent and should not be interpreted as flood risk.
- Storage figures distinguish wetland-compartment storage from whole-system storage. The accepted additional 5,500 m³ is retained in the no-prefix initial conditions. Accounting differences do not identify the origin of individual water parcels or isolate a causal initial-water effect.
- Subset ribbons describe reuse of a finite 64-vector set. Their collapse at n=64 is automatic, not convergence evidence. Reweighting changes an assumed law on the same support; ESS is a weight-concentration diagnostic, not a new simulation count or a basis for an asserted weighted confidence band.
- Antecedent plots retain the negative-time prefix and mark the original event boundary. Original event weather is unchanged; carried states and controller memory are not reset. Event-only maxima and full extended-horizon maxima are separately reported. Conditioning for 24 h or 72 h does not establish equilibrium or reconstruct observed storage.

## Presentation notes

The existing captions resolve the few conventions that need care: official scatter-matrix cells display upper interval / correlation / lower interval vertically; hydrograph scales may differ between parameter columns or event rows; and the 72 h weather plot does not claim that the separate 24 h prefix has identical first-half-hour interpolation. These are documented reading cautions, not reasons to rerender or alter old results. All links in the master index were checked against local files.
