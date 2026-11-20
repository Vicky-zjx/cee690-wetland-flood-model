# Stage D.2 pilot metric definitions and interpretation

This analysis uses saved Stage D feedback trajectories. It does not alter initial states, remove an initial segment, run a spin-up, change gate operation, regenerate weather or replace the physical equations. The primary full-window peaks retain the initial sample. Additional response metrics describe features that those peaks can conceal; they do not replace the accepted objective or define a flood threshold.

## Saved population, clocks and references

The [scenario manifest](../g5_diagnostics/recovered/scenario_manifest.json) contains 20 original six-day weather windows: 10 Copula/VAR and 10 block-bootstrap cases. They represent ordinary September weather conditioned on the target calendar, not a storm-selected or flood-selected population. Every CSV has 144 hourly intervals. The UTC timestamp starts each interval: September 9, 2022 00:00 through September 14 23:00, with the analysis ending on September 15 00:00. The modeled duration is 518,400 s = 144 h. These calendar labels locate synthetic weather; they do not make it an observed September 2022 realization.

Precipitation is an hourly mean rate in mm/h; dividing by 3,600,000 converts it to m/s. Radiation is an hourly mean in W/m². The saved runner conserves the original hourly rain amount and solar energy using controller blocks within their original intervals. Response integrals use elapsed seconds from the saved ODE grid, including both event endpoints. The nominal pilot uses the uniform 15 s suite in [fine_nominal](../feedback_control/results/fine_nominal/); any finer comparison must remain a separately labeled repeated evaluation of the same weather IDs.

The initial river elevation is **H₀ = 97.35 m** in the model datum. The corresponding saved initial outlet discharge is **Q₀ = 0.44444444444441916 m³/s**, displayed approximately as 0.4444 m³/s. Calculations use the full-precision saved values and verify agreement across the compared cases. These are references derived from the accepted assumed initial state, not observed flood levels, channel capacities or approved safety thresholds. No active downstream flood threshold has been established; the retired 97.9 m value is not restored.

## Primary peaks and supplementary response metrics

For a saved scalar trajectory y(t), all metrics use the entire interval [0,T]. Piecewise-linear interpolation between saved output nodes is the explicit postprocessing convention; it is not a statement that the continuous ODE solution is exactly linear.

| Metric | Definition | Units and interpretation |
| --- | --- | --- |
| River peak | max Hᵣ(t), including t = 0 | m in model datum; original physical maximum |
| Outlet peak | max Qout(t), including t = 0 | m³/s; reported separately from river elevation |
| Peak time | First saved time attaining the maximum | h from event start; a maximum at t = 0 is labeled an initial-condition maximum |
| Maximum rise above initial level | max[Hᵣ(t) − H₀], including t = 0 | m; nonnegative because the initial difference is zero; equals the whole-window peak minus H₀ |
| River upward excursion (secondary) | max over t of [Hᵣ(t) − min over 0 ≤ s ≤ t of Hᵣ(s)] | m; largest later rise above an earlier running minimum |
| Outlet upward excursion | The same running-minimum construction for Qout | m³/s; largest later discharge rise |
| Initial-level exceedance duration | Integral of 1[Hᵣ(t) > H₀] dt | h; time above the initial reference, not flood duration |
| Initial-discharge excess volume | Integral of max(Qout(t) − Q₀, 0) dt | m³; excess above the initial-flow reference, not flood volume |
| Total outlet volume | Integral of Qout(t) dt | m³ over the complete window, including initial-water release |
| Signed paired level integral | Integral of [Hᵣ,no(t) − Hᵣ,wet(t)] dt | m·h; positive means lower time-integrated river level with wetlands |

The saved current outlet relation is Qout = max(Hᵣ − 97.15 m, 0)/(0.45 s/m²); the extracted table records zero discrepancy for the nominal saved outputs. Consequently, stage and discharge peaks are linked by the same monotone rating, and the initial-stage and initial-discharge reference metrics are also related. Their agreement is not independent validation of the model.

The requested maximum rise above the initial level uses the unchanged H₀ and retains the initial sample. It is exactly the full-window river peak translated by H₀, so it retains the same zero-versus-positive structure and is not an independent source of information. It must not be confused with the secondary running-minimum upward excursion.

The upward-excursion metric is evaluated over the full record and does not select a shorter analysis window. Its running minimum is part of the metric definition, not a replacement initial condition. It captures a later rebound even when that rebound stays below H₀. It is not necessarily rainfall-caused: the accepted dry `block_bootstrap_00002` has 0 mm precipitation but nonzero excursions because the rest of the dynamical system still evolves. No counterfactual dry-weather run has been performed to isolate a rainfall-only response.

Fifteen of the 16 valid nominal pairs attain their whole-window maximum at the common initial level. Their zero peak differences do not imply identical hydrographs or no wetland routing effect. Conversely, it would be misleading to discard these initial peaks merely to obtain larger apparent benefits. A difference in peak times is a descriptive comparison; if both maxima occur at initialization, the resulting zero is not a measured flood-wave travel time.

## Exact threshold crossing calculations for the reported interpolant

On an interval of length d = t₁ − t₀, let z₀ = y(t₀) − c and z₁ = y(t₁) − c. Define D as the time for which the linear interpolant is **strictly** above c, and A as its integrated positive excess.

- If z₀ ≤ 0 and z₁ ≤ 0: D = 0 and A = 0.
- If z₀ ≥ 0 and z₁ ≥ 0, with at least one positive: D = d and A = d(z₀ + z₁)/2. Isolated endpoint equality contributes no duration.
- If the signs differ: let p be the positive endpoint value and w = |z₁ − z₀|. Then D = d p/w and A = d p²/(2w). This splits the interval at its interpolated crossing.

A plateau exactly at the reference has zero strict-exceedance duration. Summing D and A gives exact results for this piecewise-linear interpolant. Simply applying trapezoidal integration to already-clipped endpoint excesses overestimates a crossing triangle; trapezoidal integration of a binary exceedance indicator also misplaces the crossing. Ordinary trapezoidal integration is exact for the unthresholded linear interpolant.

For stage, A has units m·s (or m·h), **not m³**. For discharge, A has units m³. Multiplying integrated stage excess by a surface area produces m³·s, not a discharge volume. If a physical flood threshold is approved later, any volume defined by discharge while stage exceeds that threshold requires integrating discharge over the stage-defined subintervals; that is a different metric from excess discharge above Q₀. No such flood metric is calculated here.

## Pairing, signs and water accounting

Pairs must share scenario ID, seed, source-weather hash, time grid, weather mapping, configuration and shared initial storages, and must both pass the existing trajectory checks. The feedback policy is common; actual gate histories must also be checked rather than assumed equal. The nominal 15 s pilot had identical realized gates within its 16 valid wetland/no-wetland pairs. Finer runs need their own check.

Reported differences use **no wetland minus wetland**. Positive differences mean a smaller metric with wetlands, but not every metric reduction constitutes a safety improvement or a flood benefit. Negative and zero differences remain valid outcomes. Different-case maxima and running-minimum excursions can occur at different times; a difference between their scalar summaries is not a pointwise hydrograph subtraction.

The accepted wetland initial condition contains **5,500 m³ more water** than the no-wetland system; common compartments match, and no redistribution is performed. Total outlet-volume differences therefore include initial-water release, differences in end storage, evaporation and groundwater exchange as well as routing. Internal gate releases, wetland drainage and overflow must not be added together as independent external inflows. A smaller outlet volume within six days does not by itself establish permanent flood-water removal: water can remain stored or leave after the analysis horizon.

## Preliminary ECDFs and denominators

For n available valid scenario metrics Yᵢ, the descriptive ECDF is Fₙ(y) = count(Yᵢ ≤ y)/n. It is right-continuous, with tied observations producing one jump equal to their combined frequency. Use the actual ties; do not jitter, smooth or fit a continuous distribution to make the plot more varied. This matches the step-function construction described in the [official SciPy ECDF documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ecdf.html).

The nominal denominator is **16 valid pairs**, comprising **10 Copula/VAR and 6 block-bootstrap** cases. The four cold block-bootstrap IDs, 00004–00007, remain visible in the 20-ID accounting as outside the unchanged liquid-water model's applicability. They are not zero outcomes, censored flood peaks or valid low-risk cases. Dividing the 16 observed ECDF counts by 20 would yield a subdistribution that ends at 0.8, not the conditional ECDF. Successful simulation of 16/20 original IDs is an accounting fraction, not a climatic probability.

The pooled ECDF is an equal-weight description of the retained 16-case pilot. It is no longer a 50/50 generator mixture: its weights are 10/16 and 6/16. Generator-specific curves with their own denominators make this visible. The two weather algorithms are alternative fitted models, not draws from an established model-probability distribution. Independently seeded simulations can be independent conditional on a fitted generator; that does not justify treating this pooled, heterogeneous, small pilot as identically distributed observations of regional floods. Within-trajectory timesteps and repeated numerical grids are never additional statistical samples. The wetland and no-wetland outcomes are paired, not independent groups.

For a strict positive difference, its empirical fraction is count(Δ > 0)/n = 1 − Fₙ(0); count(Δ ≥ 0)/n is different when there are ties. Report negative, zero and positive counts separately. No return period, annual exceedance probability, calibrated population distribution, inferential confidence band or proof of benefit follows from these pilot ECDFs.

## Initial-state sensitivity and possible future spin-up

The current study remains conditional on the accepted initial states. The weather generator's 336 h latent burn-in establishes neither hydrologic equilibration nor a reconstructed soil, groundwater, pool or river condition. The [USACE HEC-HMS Applications Guide](https://www.hec.usace.army.mil/confluence/hmsdocs/hmsag/introduction/how-should-hec-hms-be-used) emphasizes application-appropriate initial states, possible initialization periods for continuous simulation, and sensitivity of response to antecedent conditions. Its guidance supports examining initialization; it does not supply a validated spin-up duration or initial state for this hypothetical model.

For a separately approved future sensitivity study, use a documented antecedent weather period and the same boundaries and operating policy for each pair. Evolve and save every hydrologic state and the controller memory/actual gate position at the event boundary. Choose a lead-in length and initial-state sensitivity criterion before examining benefit outcomes, and compare plausible starting states or longer lead-ins. Do not reset river storage alone, repeat the target event backward as if observed, or choose a cutoff from the resulting flood peak. A longer run does not automatically establish equilibrium under changing weather. The feasibility of that future procedure depends on available antecedent forcing, treatment of freezing, and an explicitly approved initial-condition comparison; none is silently supplied here.
