# Official parameter-uncertainty plots

These plots use the saved full-ODE results from the approved iid exploratory
parameter study. They do not use the 27-point factorial grid, the earlier
16-weather-scenario pilot, a surrogate fit or a flood threshold. Positive paired
differences mean no-wetland minus wetland. Negative differences remain results.

The six distribution plots show wetland, no-wetland and paired full-window peak
stage and outlet discharge. The unchanged official plotting utility uses
rank/(N+1) plotting positions and 95% pointwise Beta order-statistic limits,
not a ties-aggregated count/N ECDF or a simultaneous confidence band. Tied/atomic
outputs do not have the usual exact continuous-order-statistic coverage claim.
A conventional sample median is corrected for display only; returned plotting
positions and limits are unchanged. The coefficient-of-variation annotation is
not informative for a signed effect with a near-zero mean.

The five-variable scatter figure contains the three physical-parameter
multipliers and the two signed paired peak differences. Its upper panels give
Pearson correlations and approximate Fisher-z limits (upper limit displayed
above, lower limit below). These are association diagnostics, not Sobol indices,
causal effects or calibrated local engineering risk. Constant columns, if any,
are omitted only from correlation panels and remain in CSV/distribution outputs.

All statements are conditional on the declared exploratory log-uniform sampling
laws, fixed storm, geometry, feedback policy and initial conditions. The saved
metadata records the exact planned/valid/invalid denominator, any rejected-draw
causes, data/configuration hashes, official revision and plotting warnings.
Unsafe but physically valid runs are not removed. The distinct initial water
volumes and realization-specific outlet resistance remain part of the accepted
paired structural comparison.
