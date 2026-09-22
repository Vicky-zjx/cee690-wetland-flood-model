# Eight-state wetland model — Stage A student summary

## Purpose and boundary

This model explores how two wetland–dam branches store and route water toward a downstream river. It retains the accepted eight ODEs and existing closures. The audit found consistent signs, units, state/input/output order, parameter names, and equation-to-code mapping; no equation change was needed. The example is **synthetic and uncalibrated**, not a prediction of real flood risk.

```text
Upland → Dam 1 → Gate 1 → Wetland 1 → River → Outside
Upland → Dam 2 → Gate 2 → Wetland 2 → River
Upland → River
Upland → Soil → Groundwater → River
Wetland 1 and Wetland 2 → Groundwater
External groundwater ↔ Groundwater; external inflow → River
Dam spills → corresponding wetlands; wetland overflows → River
```

Surface flow is one-way, although high receiving levels can suppress it. External groundwater exchange is signed. Rain falls on six disjoint surface footprints totaling 660,000 m²; soil and aquifer do not receive duplicate rainfall. Blockage affects only the final river outlet.

## States, decisions, and disturbances

Python indices are zero-based. All eight states are water volumes in m³; time is seconds.

| Index | State | Meaning |
|---:|---|---|
| 0 | S_u | Upland surface water |
| 1 | S_m | Upland soil water |
| 2 | S_w1 | Wetland 1 surface water |
| 3 | S_w2 | Wetland 2 surface water |
| 4 | S_d1 | Dam pool 1 water |
| 5 | S_d2 | Dam pool 2 water |
| 6 | S_g | Drainable groundwater above its lower reference |
| 7 | S_r | River-reach water |

Four conceptual decisions are Gate 1 opening, Gate 2 opening, Wetland 1 design, and Wetland 2 design. Openings g_1(t),g_2(t) vary between 0 and 1. Each wetland design combines categorical type τ_i and fractional cover ρ_i∈[0,1], fixed throughout a simulation. `type_a` and `type_b` are illustrative labels, not measured species traits.

The numerical-input order is `[g1,g2,P,T_a,R_s,v,Q_ext,h_d,h_g,B]`. The last eight remain external disturbances: precipitation [m/s], temperature [°C], radiation [W/m²], wind [m/s], external inflow [m³/s], downstream level [m], external groundwater level [m], and blockage [dimensionless]. Boundary signals are independently prescribed; Q_ext excludes runoff already generated inside the modeled footprints. External h_g differs from internal H_g.

## Eight water-balance ODEs

Every P A, Q, and E term below has units m³/s. Positive Q_gx enters groundwater. Each of the 16 internal transfers appears once negatively at its donor and once positively at its receiver.

1. dS_u/dt = P A_u − Q_inf − Q_ud1 − Q_ud2 − Q_ur.
2. dS_m/dt = Q_inf − Q_mg − E_m.
3. dS_w1/dt = P A_w1 + Q_G1 + Q_dsp1 − Q_wr1 − Q_wsp1 − Q_wg1 − E_w1.
4. dS_w2/dt = P A_w2 + Q_G2 + Q_dsp2 − Q_wr2 − Q_wsp2 − Q_wg2 − E_w2.
5. dS_d1/dt = P A_d1 + Q_ud1 − Q_G1 − Q_dsp1 − E_d1.
6. dS_d2/dt = P A_d2 + Q_ud2 − Q_G2 − Q_dsp2 − E_d2.
7. dS_g/dt = Q_mg + Q_wg1 + Q_wg2 + Q_gx − Q_gr.
8. dS_r/dt = P A_r + Q_ur + Q_wr1 + Q_wr2 + Q_wsp1 + Q_wsp2 + Q_gr + Q_ext − Q_out − E_r.

`system(t,x,u,c)` evaluates each flux once and returns `dxdt,y`. Levels and flows are algebraic outputs; peaks and event volumes are postprocessed metrics. The named channel lists document all eight states, ten inputs, and 37 outputs.

## Relationships and figures to study

Stage is H=z+S/A; groundwater uses H_g=z_g+S_g/(S_y A_g). Availability φ=S/(S+V*) reduces withdrawals near empty storage. Upland runoff uses capped kS rates; recharge uses k_mg[S_m−S_fc]+. Wetland seepage and groundwater discharge use positive head differences, conductance, and donor availability.

Design interpolates resistance and ET coefficients between bare and selected vegetation endpoints: R_w=(1−ρ)R_bare+ρR_veg,τ, with the same interpolation for K_ET,w. Baseline R_w=[2.1,2.2] s/m². These assumed mappings affect drainage and ET, not every wetland parameter.

The ten individually saved figures are listed below. Comparison values in A–G are teaching examples, not empirical ranges. Positive part means [z]+=max(z,0).

| Figure in `plots/` | Axes and varied settings | What to learn; ODE use |
|---|---|---|
| **A_gate_release.png** | Effective head ΔH=H_d−max(H_w,z_G) [m] versus Q_G [m³/s]; g=.25,.50,.75,1. | Q_G=g min(Q̄_G,[ΔH]+/R_G): zero against unfavorable head; opening increases slope and cap. Transfers pool water into its wetland (3–6). |
| **B_wetland_drainage.png** | ΔH=H_w−max(H_r,z_wo) [m] versus Q_wr [m³/s]; R_w=.8,2.10,2.88 s/m², representing bare/type_a/type_b at stated cover. | At the same effective head, higher resistance reduces uncapped discharge and delays the shared .65 m³/s cap. Both capped curves have equal discharge; this graph does not establish lower flood peaks. Transfers wetland water to river (3,4,8). |
| **C_wetland_overflow.png** | ΔH=H_w−max(z_wsp,H_r) [m] versus Q_wsp [m³/s]; c=1.5,3,4.5 m^(3/2)/s. | c[ΔH]+^(3/2) is zero before onset and steepens afterward; no explicit flow cap. An additional physical path removes wetland water (3,4,8). |
| **D_upland_infiltration.png** | Soil fraction S_m/C_m and upland storage S_u [m³] versus Q_inf [m³/s]; K_u=.5,1,2×10⁻⁶ m/s. | A_u K_u φ_u[1−S_m/C_m]+ vanishes at saturation or zero available water; bounded by infiltration capacity. Water moves upland→soil (1,2). |
| **E_evaporative_demand.png** | Temperature [°C] versus E_star [mm/h]; radiation=100,300,600 W/m², wind=2 m/s. | Radiation supplies energy; temperature changes vapor-pressure terms and weighting. Wind enters transfer conductance. No imposed cap; actual losses also require source water. ET affects (2–6,8). |
| **F_river_outlet_blockage.png** | ΔH=H_r−max(h_d,z_ro) [m] versus Q_out [m³/s]; B=0,.25,.50,.75. | Blockage scales slope and cap by 1−B; unfavorable head gives zero. Caps are 3,2.25,1.5,.75 m³/s. This removes river water (8). |
| **G_groundwater_exchange.png** | h_g−H_g [m] versus signed Q_gx [m³/s]; C_gx=.0175,.035,.07 m²/s. | Positive head gives inward flow; negative head gives outward flow reduced by φ_g. Zero crossing is head equality; no explicit cap. Boundary exchange enters (7). |
| **H_baseline_storages.png** | Time [hours] versus all eight storages [m³]. | Follow storage accumulation and recession across compartments. |
| **I_baseline_response.png** | Time versus river level, outlet discharge, gate releases, and wetland flows. | Connect downstream timing with upstream releases, drainage, and overflow. |
| **J_baseline_inputs.png** | Time versus both gates and eight disturbance channels, each with units. | Distinguish prescribed causes from calculated responses. |

Thresholds and caps prevent these particular flow laws from releasing water against their allowed direction or beyond rated normal capacity. Spill thresholds are onset levels, not hard storage bounds; overflow may coexist with normal drainage.

## Feedback mechanisms

Wetland storage rises → wetland level rises → normal drainage increases → this opposes further storage growth, **when downstream level and resistance are fixed, positive head permits drainage, and the cap is not active**. Below the outlet threshold or at the cap, that local drainage response is flat.

Soil storage approaches capacity → infiltration decreases → this opposes further soil-water accumulation, **holding surface-water supply and other infiltration parameters fixed**. At capacity infiltration is zero; when no surface water is available, infiltration is already zero.

## Parameters, workflow, and baseline

`parameters.md` gives symbols, Python names, units, current values, categories, admissible domains, and future data needs. No defensible site-specific numerical ranges are claimed. Physical constants, illustrative geometry/coefficients, initial conditions, designs, and disturbances are separated there.

Run from this folder:

```sh
python3 wetland_flood_model.py
```

The workflow is **define model → generate A–G → run one baseline with ode4u → generate H–J → interpret**. For algebraic figures alone, run `python3 algebraic_plots.py`; this requires NumPy and Matplotlib, not the solver. The baseline additionally requires NumPy 2+ and official `multivarious` with its dependencies, using `work/multivarious` when present. `load_ode4u()` loads the solver only when requested; missing dependencies produce a clear message without installation or substitution. Importing either script runs no simulation. The main retains `class c`, `system(t,x,u,c)`, and the professor's execution block. Solver arrays use `(channels,time)` and linearly interpolated inputs.

When editing geometry: **edit geometry → explicitly call `params.refresh_geometry_defaults()` if appropriate → construct `x_init` and inputs → validate → run**. This deliberately refreshes only the five linked defaults documented in `parameters.md`; independent empirical parameters remain manually editable.

The baseline lasts six days with 120 s steps. Two triangular storms deliver 60 and 100 mm. Both gates vary independently; vegetation remains fixed. Initial storage is `[600,90000,3000,2500,3000,2400,409600,7000]` m³. Exact forcing histories are in `parameters.md`.

The fresh run gives peak outlet discharge **1.357896 m³/s** and peak river level **97.879407 m**, below the illustrative 97.90 m reporting threshold. Both dam spills and wetland overflows activate. Stored water and release timing explain the delayed responses; this single baseline cannot establish a restoration benefit.

Numerical verification remained concise: execution and plots completed, states remained admissible, and maximum integrated water-balance residual was **0.076735 m³**, or **0.00004151%** of gross event inflow. `results/baseline_summary.json` records the current diagnostics. Numerical consistency does not establish physical validity; no scenario suite runs in the main workflow.

## Stage A checklist and remaining work

Stage A numerical workflow completed; the remaining-work checklist still applies.

- [x] Preserve and audit eight states, topology, balances, units, and indices.
- [x] Separate controls, designs, disturbances, parameters, and outputs.
- [x] Explain and plot A–G; retain one readable baseline and H–J.
- [x] Document illustrative parameters and numerical limitations.
- [ ] Evidence-supported site-specific parameter ranges: pending suitable evidence; mathematical domains do not satisfy this requirement.
- [ ] Your work: run the files, explain each curve and water balance in your own words, and make a small plot-parameter edit to demonstrate understanding.

Fixed areas, empirical ratings, one-way internal exchanges, synthetic vegetation endpoints, and simplified ET limit interpretation. No measured geometry or calibrated parameters support real flood predictions. Later course stages cover design optimization (Part B), data acquisition (Part C), and uncertainty propagation (Part D); none is implemented here.
