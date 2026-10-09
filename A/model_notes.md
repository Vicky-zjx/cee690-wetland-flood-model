# Eight-state wetland model — HW A.2

## Purpose and boundary

For **CEE690, Systems Modeling for Climate Response**, HW A.2 applies the ODE module to a project model. This eight-state model explores how two wetland–dam branches store and route water toward a downstream river. The example is **synthetic, illustrative, and uncalibrated**; it demonstrates model behavior rather than predicting real flood risk.

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

The controls are the two gate-opening histories g_1(t),g_2(t)∈[0,1]. Each fixed wetland design combines categorical type τ_i and fractional cover ρ_i∈[0,1], selected independently for each branch and held constant throughout a simulation. `type_a` and `type_b` are illustrative labels, not measured species traits.

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

Summing the eight equations cancels every internal transfer: dS_tot/dt = P A_rain + Q_ext + Q_gx − Q_out − E_tot. The baseline compares total storage change with the time integral of this external net flux.

## Implemented algebraic relationships

Positive part means [z]+=max(z,0). Stage is H=z+S/A; groundwater uses H_g=z_g+S_g/(S_y A_g), with specific yield S_y. Availability φ_j=S_j/(S_j+V_j*) reduces withdrawals near empty storage.

For branch i=1,2, upland routing is Q_udi=min(Q̄_udi,k_udi S_u), with Q_ur=min(Q̄_ur,k_ur S_u) for the direct river bypass. Infiltration is Q_inf=A_u K_u φ_u[1−S_m/C_m]+; recharge is Q_mg=k_mg[S_m−S_fc]+, starting above the assumed field-capacity storage S_fc.

The common capped rating is F(H_a,H_b,z_o,R,Q̄)=min(Q̄,[H_a−max(H_b,z_o)]+/R). Gate release is Q_Gi=g_i F(H_di,H_wi,z_Gi,R_Gi,Q̄_Gi); wetland drainage is Q_wri=F(H_wi,H_r,z_woi,R_wi,Q̄_wri); the final outlet is Q_out=(1−B)F(H_r,h_d,z_ro,R_r,Q̄_r). A, B, and F below illustrate these relationships.

Dam spill is Q_dspi=c_dspi[H_di−max(z_dspi,H_wi)]+^(3/2); wetland overflow is Q_wspi=c_wspi[H_wi−max(z_wspi,H_r)]+^(3/2). Spill thresholds mark flow onset, not hard storage bounds; overflow can coexist with normal drainage.

Wetland seepage is Q_wgi=C_wgi φ_wi[H_wi−H_g]+; groundwater discharge is Q_gr=C_gr φ_g[H_g−H_r]+. External exchange is Q_gx=C_gx[h_g−H_g]+−C_gx φ_g[H_g−h_g]+, so only its outward term is reduced by groundwater availability.

Design interpolates resistance and ET coefficients: R_wi=(1−ρ_i)R_bare,i+ρ_i R_veg,τi,i, with the same interpolation for K_ET,wi. Baseline R_w=[2.1,2.2] s/m². These assumed mappings affect drainage and ET, not every wetland parameter.

Reduced evaporative demand combines an energy term and an aerodynamic term: E_star = Δ/(Δ+γ) · ηR_s/(ρ_w λ_v) + γ/(Δ+γ) · k_a D_ρ/ρ_w [m/s]. Here Δ is the saturation-vapor-pressure slope, γ the psychrometric constant, D_ρ the vapor-density deficit, and k_a=k_0+k_v v the wind-dependent conductance; `evaporative_demand()` implements their thermodynamic expressions. Actual soil, wetland, dam, and river losses use E_j=A_j K_ET,j φ_j E_star [m³/s], with A_u for soil. This simplified demand is not standard reference ET₀.

## Figure guide

The ten existing figures are saved individually in `A/plots/`. Comparison values in A–G are teaching examples, not empirical ranges.

| Figure in `A/plots/` | Axes and varied settings | What to learn; ODE use |
|---|---|---|
| **A_gate_release.png** | Effective head ΔH=H_d−max(H_w,z_G) [m] versus Q_G [m³/s]; g=.25,.50,.75,1. | Q_G=g min(Q̄_G,[ΔH]+/R_G): zero against unfavorable head; opening increases slope and cap. Transfers pool water into its wetland (3–6). |
| **B_wetland_drainage.png** | ΔH=H_w−max(H_r,z_wo) [m] versus Q_wr [m³/s]; R_w=.8,2.10,2.88 s/m², representing bare/type_a/type_b at stated cover. | At the same effective head, higher resistance reduces uncapped discharge and delays the shared .65 m³/s cap. Once capped, the curves have equal discharge; this graph does not establish lower flood peaks. Transfers wetland water to river (3,4,8). |
| **C_wetland_overflow.png** | ΔH=H_w−max(z_wsp,H_r) [m] versus Q_wsp [m³/s]; c=1.5,3,4.5 m^(3/2)/s. | c[ΔH]+^(3/2) is zero before onset and steepens afterward; no explicit flow cap. An additional physical path removes wetland water (3,4,8). |
| **D_upland_infiltration.png** | Soil fraction S_m/C_m and upland storage S_u [m³] versus Q_inf [m³/s]; K_u=.5,1,2×10⁻⁶ m/s. | A_u K_u φ_u[1−S_m/C_m]+ vanishes at saturation or zero available water; bounded by infiltration capacity. Water moves upland→soil (1,2). |
| **E_evaporative_demand.png** | Temperature [°C] versus E_star [mm/h]; radiation=100,300,600 W/m², wind=2 m/s. | Radiation supplies energy; temperature changes vapor-pressure terms and weighting. Wind enters transfer conductance. No imposed cap; actual losses also require source water. ET affects (2–6,8). |
| **F_river_outlet_blockage.png** | ΔH=H_r−max(h_d,z_ro) [m] versus Q_out [m³/s]; B=0,.25,.50,.75. | Blockage scales slope and cap by 1−B; unfavorable head gives zero. Caps are 3,2.25,1.5,.75 m³/s. This removes river water (8). |
| **G_groundwater_exchange.png** | h_g−H_g [m] versus signed Q_gx [m³/s]; C_gx=.0175,.035,.07 m²/s. | Positive head gives inward flow; negative head gives outward flow reduced by φ_g. Zero crossing is head equality; no explicit cap. Boundary exchange enters (7). |
| **H_baseline_storages.png** | Time [hours] versus all eight storages [m³]. | Follow storage accumulation and recession across compartments. |
| **I_baseline_response.png** | Time versus river level, outlet discharge, gate releases, and wetland flows. | Connect downstream timing with upstream releases, drainage, and overflow. |
| **J_baseline_inputs.png** | Time versus both gates and eight disturbance channels, each with units. | Distinguish prescribed causes from calculated responses. |

## Feedback mechanisms

Wetland storage rises → wetland level rises → normal drainage increases → this opposes further storage growth, **when downstream level and resistance are fixed, positive head permits drainage, and the cap is not active**. Below the outlet threshold or at the cap, that local drainage response is flat.

Soil storage approaches capacity → infiltration decreases → this opposes further soil-water accumulation, **holding surface-water supply and other infiltration parameters fixed**. At capacity infiltration is zero; when no surface water is available, infiltration is already zero.

## Parameters, workflow, and baseline

`class c` in [wetland_flood_model.py](wetland_flood_model.py) defines parameter values, units, fixed designs, initial conditions, and forcing settings. `make_inputs()` constructs the prescribed input histories; `validate()` checks mathematical domains and assumed geometry limits. These checks do not establish evidence-supported site ranges.

After the environment setup in [A/README.md](README.md), run either command from the project root:

```sh
python3 -m A.algebraic_plots       # A–G only; NumPy and Matplotlib
python3 -m A.wetland_flood_model    # A–J and one baseline; official ode4u required
```

The second command follows **define model → generate A–G → run one baseline with official ode4u → generate H–J → interpret**. It additionally requires NumPy 2+ and official `multivarious` with its dependencies, using `A/work/multivarious` when present. `load_ode4u()` imports the solver only for simulation; importing either project script runs no simulation. The code follows the professor's `class c`, `system(t,x,u,c)`, and execution-block structure. Solver arrays use `(channels,time)` and linearly interpolated inputs.

When editing geometry: **edit geometry → explicitly call `params.refresh_geometry_defaults()` if appropriate → construct `x_init` and inputs → validate → run**. The helper updates only A_g, V_w, V_d, V_g, and V_r from their documented assumptions in `class c`; independent empirical parameters remain manually editable.

The baseline lasts six days with 120 s steps. Two triangular storms deliver 60 and 100 mm. Both gates vary independently; vegetation remains fixed. Initial storage is `[600,90000,3000,2500,3000,2400,409600,7000]` m³. Exact forcing knots and waveforms are in `class c` and `make_inputs()`.

The existing baseline gives peak outlet discharge **1.357896 m³/s** and peak river level **97.879407 m**, below the illustrative 97.90 m reporting threshold. Both dam spills and wetland overflows activate. Stored water and release timing explain the delayed responses; this single baseline cannot establish a restoration benefit.

Figures go to `A/plots/`; numerical arrays go to `A/results/baseline.npz` and diagnostics to `A/results/baseline_summary.json`. The existing run's maximum integrated water-balance residual was **0.076735 m³**, or **0.00004151%** of gross event inflow. This records numerical consistency, not physical validation. Rerunning overwrites the same output files.

## Limitations and A.2 remaining-work checklist

Fixed areas, empirical ratings, one-way internal exchanges, synthetic vegetation endpoints, and simplified ET limit interpretation. Geometry, site coefficients, and forcing histories are illustrative; evidence-supported site ranges remain unavailable. HW A.2 uses these relationships and one baseline to explain the ODE model. No calibration, physical validation, or design optimization is claimed.

- [ ] Explain the A–G curves, H–J baseline behavior, and total water balance in your own words.
- [ ] Identify the illustrative assumptions and explain how they limit interpretation.


## Complete numerical provenance for the accepted A example

This register covers every assigned input in `class c`, the derived geometry and initial states, and the additional algebraic-plot comparison settings. The authoritative values and formulas are in [the parameter definition](wetland_flood_model.py#L51), [the state/geometry properties](wetland_flood_model.py#L189), [the evaporation calculation](wetland_flood_model.py#L236), and [the synthetic input generator](wetland_flood_model.py#L311). These are **implementation sources**, not observational evidence. No A input is a local measurement, a reconstructed Milwaukee observation, or a calibrated coefficient. In particular, the 0.66 km² footprint is an invented demonstration installation, not the Underwood Creek gauge catchment.

The provenance categories used below are:

- **Assumed:** an explicit synthetic physical, engineering, or operating choice made for the course demonstration; empirical support for this particular value is unavailable.
- **Derived:** calculated by a stated formula from other listed values; derivation does not turn assumptions into observations.
- **Physical approximation / external formula:** a general thermodynamic value or relationship with the narrowly stated source below. It does not validate the complete ET method or any site-specific coefficient.
- **Numerical / reporting:** a unit conversion, time discretization, check tolerance, or illustrative reporting threshold; not a measured physical parameter or certified safety criterion.

Branch arrays always list **[branch 1, branch 2]**. All elevations share an **assumed common datum** and are not surveyed elevations. Historical A retains its empirical linear gate law and independent dam-spillway paths. The revised no-dam-spillway scope and gate applicability checks in B do not retroactively change these A values or equations.

### Geometry, soil, routing, and surface outlets

| Code field | Accepted A value and units | Provenance, derivation, and applicability |
|---|---|---|
| `A_u` | 600,000 m² | Assumed upland rainfall footprint. |
| `A_d` | [10,000, 8,000] m² | Assumed constant dam-pool surface areas. |
| `A_w` | [12,000, 10,000] m² | Assumed constant wetland surface areas. |
| `A_r` | 20,000 m² | Assumed river-reach rainfall and storage area. |
| `A_g` | 640,000 m² | Derived as `A_u + sum(A_d) + sum(A_w)`; assumed aquifer footprint. It receives no additional direct rainfall. |
| `A_rain` (property) | 660,000 m² = 0.66 km² | Derived sum of the six disjoint surface footprints, including `A_r`. |
| `S_y` | 0.20, dimensionless | Assumed drainable aquifer specific yield; no local aquifer test. |
| `z_d` | [100.0, 99.8] m | Assumed dam-pool bed elevations. |
| `z_w` | [98.0, 97.9] m | Assumed wetland bed elevations. |
| `z_r`, `z_g` | 97.0 m; 94.0 m | Assumed river bed and groundwater lower reference, respectively. |
| `C_m` | 120,000 m³ | Derived from `A_u × 0.20 m`; 0.20 m is an assumed equivalent soil-water capacity, not measured soil thickness. |
| `S_fc` | 60,000 m³ | Derived from `A_u × 0.10 m`; assumed recharge-onset/field-capacity storage. |
| `K_u` | 1 × 10⁻⁶ m/s | Assumed effective infiltration-rate scale; not a soil survey value or measured saturated conductivity. |
| `k_mg` | 1/(5 days) = 2.314815 × 10⁻⁶ s⁻¹ | Assumed five-day recharge time scale. |
| `k_ud` | [1/(8 h), 1/(10 h)] = [3.472222, 2.777778] × 10⁻⁵ s⁻¹ | Assumed upland-to-pool routing time scales. |
| `k_ur` | 1/(20 h) = 1.388889 × 10⁻⁵ s⁻¹ | Assumed direct upland-to-river routing time scale. |
| `Q_udcap`, `Q_urcap` | [1.2, 1.0] m³/s; 0.7 m³/s | Assumed capacities for the two pool branches and direct river bypass. |
| `z_G` | [100.05, 99.85] m | Assumed gate outlet elevations, 0.05 m above the corresponding pool beds. |
| `R_G` | [1.6, 1.8] s/m² | Assumed empirical linear gate resistance. A does not derive this from an orifice geometry or a measured rating curve. |
| `Q_Gcap` | [1.0, 0.8] m³/s | Assumed full-opening gate discharge caps; actual cap is multiplied by `g_i`. |
| `z_dsp` | [100.9, 100.7] m | Assumed independent dam-spillway onset elevations retained in A. These are not storage clips or B safety limits. |
| `c_dsp` | [4.0, 3.2] m^(3/2)/s | Assumed effective coefficients in the A dam-spillway 3/2-power ratings; not measured or independently geometry-derived. |
| `C_d` (property) | [9,000, 7,200] m³ | Derived as `A_d × (z_dsp − z_d)`; storage at nominal spill onset, not a hard capacity. |
| `z_wo` | [98.05, 97.95] m | Assumed normal wetland outlet elevations, 0.05 m above beds. |
| `z_wsp` | [98.6, 98.5] m | Assumed wetland overflow onset elevations, 0.60 m above beds. |
| `c_wsp` | [3.0, 2.5] m^(3/2)/s | Assumed effective wetland-overflow coefficients; separate from `c_dsp`. |
| `Q_wrcap` | [0.65, 0.55] m³/s | Assumed normal wetland drainage caps. Overflow remains an additional path. |
| `C_w` (property) | [7,200, 6,000] m³ | Derived as `A_w × (z_wsp − z_w)`; storage at nominal overflow onset, not a hard capacity. |
| `z_ro` | 97.1 m | Assumed river outlet elevation, 0.10 m above its bed. |
| `R_r`, `Q_rcap` | 0.45 s/m²; 3.0 m³/s | Assumed river outlet resistance and unblocked cap; no calibrated river rating. |

### Fixed wetland classes, groundwater exchange, and availability

`type_a` and `type_b` are synthetic hydraulic/ET response classes. The assumed linear interpolation in cover is part of the demonstration, not an empirical species relationship. Zero cover retains the wetland storage compartment and its routing paths.

| Code field | Accepted A value and units | Provenance, derivation, and applicability |
|---|---|---|
| `R_bare` | [0.8, 1.0] s/m² | Assumed zero-cover normal-drainage resistance endpoints. |
| `R_veg['type_a']` | [2.8, 3.0] s/m² | Assumed full-cover class-A resistance endpoints. |
| `R_veg['type_b']` | [4.0, 4.5] s/m² | Assumed full-cover class-B resistance endpoints; not an observed vegetation effect. |
| `K_bare` | [1.0, 1.0], dimensionless | Assumed zero-cover ET coefficients. |
| `K_veg['type_a']` | [1.1, 1.1], dimensionless | Assumed full-cover class-A ET coefficients. |
| `K_veg['type_b']` | [1.25, 1.25], dimensionless | Assumed full-cover class-B ET coefficients. |
| `tau` | [`type_a`, `type_a`] | Assumed fixed baseline class selection (A/A means class A in both wetlands). |
| `rho` | [0.65, 0.60], dimensionless | Assumed fixed cover fractions in the A baseline; no vegetation survey. |
| `design()` outputs | `R_w` = [2.1, 2.2] s/m²; `K_ETw` = [1.065, 1.06] | Derived independently per branch by `(1−rho) × bare + rho × selected-class endpoint`. |
| `C_wg` | [0.015, 0.012] m²/s | Assumed wetland-to-groundwater conductances; no permeability test or fitted flux. |
| `C_gr` | 0.08 m²/s | Assumed groundwater-to-river conductance. |
| `C_gx` | 0.035 m²/s | Assumed signed external-groundwater exchange conductance. |
| `V_u` | 600 m³ | Derived from the assumed 0.001 m availability depth times `A_u`. |
| `V_m` | 6,000 m³ | Derived from the assumed 0.010 m availability depth times `A_u`. |
| `V_w` | [120, 100] m³ | Derived from assumed 0.010 m surface depth times `A_w`. |
| `V_d` | [100, 80] m³ | Derived from assumed 0.010 m surface depth times `A_d`. |
| `V_g` | 6,400 m³ | Derived from assumed 0.05 m groundwater head times `S_y × A_g`. |
| `V_r` | 200 m³ | Derived from assumed 0.010 m surface depth times `A_r`. |

The availability volumes are constitutive assumptions in `S/(S+V*)`: withdrawal is half available at `S=V*`. They are not numerical corrections and do not add or remove stored water.

### Evaporative-demand constants and coefficients

The sources in this table support only the stated general constants or formula. They do not supply the project's assumed humidity, radiation fraction, conductance, wetland classes, or ET multipliers, and do not make this reduced demand identical to FAO reference ET₀.

| Code field or formula constant | Accepted A value and units | Provenance, derivation, and applicability |
|---|---|---|
| `lambda_v` | 2.45 × 10⁶ J/kg | General constant approximation, 2.45 MJ/kg, described in [FAO-56 Chapter 3](https://www.fao.org/4/x0490e/x0490e07.htm). Temperature dependence is omitted. |
| `c_p` | 1013 J/(kg K) | General constant approximation, 1.013 kJ/(kg °C), in [FAO-56 Chapter 3](https://www.fao.org/4/x0490e/x0490e07.htm). |
| `epsilon` | 0.622, dimensionless | Ratio of water-vapor and dry-air molecular weights used in the psychrometric constant; [FAO-56 Chapter 3](https://www.fao.org/4/x0490e/x0490e07.htm). |
| `rho_w` | 1000 kg/m³ | Constant-density approximation for water, adopted in the model; no local density measurement. |
| `R_v` | 461.5 J/(kg K) | General water-vapor gas constant listed in [the UCAR-hosted JRA-25 model constants](https://osdata.gdex.ucar.edu/web/datasets/d625000/docs/jra-25_model_constants.html); not a basin measurement. |
| `p_0` | 101,325 Pa | Assumed constant reference atmospheric pressure. The model does not use station pressure or elevation adjustment. |
| Saturation vapor-pressure constants in `evaporative_demand()` | 610.8 Pa; 17.27; 237.3 °C | `e_s = 610.8 exp(17.27 T/(T+237.3))` is the Pa form of the relationship in [FAO-56 Chapter 3](https://www.fao.org/4/x0490e/x0490e07.htm); its analytical derivative supplies `Delta`. |
| Temperature conversion | `T_K = T_a + 273.15` | Celsius-to-kelvin unit conversion for the ideal-gas vapor-density calculation. |
| `r_H` | 0.65, dimensionless | Assumed fixed relative humidity; no station observation. |
| `eta` | 0.60, dimensionless | Assumed fraction of incoming radiation available to the reduced energy term; not a measured net-radiation budget. |
| `k_0` | 0.001 m/s | Assumed zero-wind aerodynamic conductance component. |
| `k_v` | 0.002, dimensionless | Assumed slope in `k_a = k_0 + k_v v`; no aerodynamic calibration. |
| `K_ETm` | 0.85, dimensionless | Assumed soil ET multiplier. |
| `K_ETd` | [1.0, 1.0], dimensionless | Assumed dam-pool ET multipliers. |
| `K_ETr` | 1.0, dimensionless | Assumed river ET multiplier. |

### Initial conditions and their origins

All initial wetness and water levels are assumptions. They do not reconstruct a historical storm. The [initial-state property](wetland_flood_model.py#L210) converts them to the eight volumes without a separate state adjustment.

| Code field | Assumed condition | Derived initial state |
|---|---|---|
| `u_depth_init` | 0.001 m over `A_u` | `S_u(0)` = 600 m³. |
| `soil_fraction_init` | 0.75 of `C_m` | `S_m(0)` = 90,000 m³. |
| `w_depth_init` | [0.25, 0.25] m above wetland beds | [`S_w1(0)`, `S_w2(0)`] = [3,000, 2,500] m³; levels [98.25, 98.15] m. |
| `d_depth_init` | [0.30, 0.30] m above dam beds | [`S_d1(0)`, `S_d2(0)`] = [3,000, 2,400] m³; levels [100.30, 100.10] m. |
| `H_g_init` | 97.20 m on the assumed datum | `S_g(0) = (97.20−94.0) × 0.20 × 640,000 = 409,600` m³. |
| `r_depth_init` | 0.35 m above river bed | `S_r(0)` = 7,000 m³; river level 97.35 m. |

Thus `x_init = [600, 90000, 3000, 2500, 3000, 2400, 409600, 7000]` m³, totaling 518,100 m³. Small machine-roundoff differences in calculated groundwater storage do not represent different initial conditions.

### Exact synthetic forcing and operating schedules

Let `h` be elapsed hours, with no calendar-date or local-time interpretation. Define a triangular pulse `Tri(h; a,b,c,M)` that is zero outside `[a,c]`, rises linearly from zero at `a` to `M` at `b`, and falls linearly to zero at `c`. Every entry below is an **assumed synthetic signal**, not station, grid, gauge, or observed operational data. Inputs are generated on the 120 s grid and the official solver interpolates between those input nodes.

| Code field(s) | Values and implemented signal | Origin / applicability |
|---|---|---|
| `rain_pulses` | `[(18,22,26,15), (66,71,76,20)]`; start/peak/end in h, amplitude in mm/h. `P = (Tri_1+Tri_2) × 10⁻³/3600` m/s. | Synthetic storms with depths `0.5×8×15 = 60` mm and `0.5×10×20 = 100` mm. The total 160 mm gives 105,600 m³ of direct rain over `A_rain`. No return period is assigned. |
| `Q_ext_base`, `ext_pulses` | Base 0.08 m³/s plus pulses `(20,25,36,0.35)` and `(68,76,90,0.70)`, with h and m³/s units. | Synthetic external river inflow from outside the six rainfall footprints. Integrated six-day volume: 79,272 m³. It is not observed gauge discharge or internally generated runoff added twice. |
| `weather_period_h`, `temperature_phase_h`, `T_mean`, `T_amplitude` | 24 h; 9 h; 22 °C; 4 °C. `T_a(h)=22+4 sin[2π(h−9)/24]` °C. | Assumed daily temperature range 18–26 °C, peaking at elapsed 15 h plus multiples of 24 h. |
| `sunrise_phase_h`, `R_s_peak` | 6 h; 550 W/m². `R_s(h)=550 max(sin[2π(h−6)/24],0)` W/m², using `weather_period_h=24`. | Assumed half-sine daytime radiation; the phase is relative elapsed time, not an observed sunrise. |
| `wind_period_h`, `wind_mean`, `wind_amplitude` | 36 h; 2.0 m/s; 0.6 m/s. `v(h)=2+0.6 sin(2πh/36)` m/s. | Assumed wind range 1.4–2.6 m/s. |
| `h_d_base`, `h_d_pulse` | 97.15 m plus `Tri(h;66,77,94,0.10)` m. | Assumed downstream stage reaches 97.25 m; prescribed independently of the simulated river. |
| `groundwater_period_h`, `h_g_base`, `h_g_amplitude` | 144 h; 97.15 m; 0.025 m. `h_g(h)=97.15+0.025 sin(2πh/144)` m. | Assumed external groundwater signal in [97.125,97.175] m; distinct from simulated internal `H_g`. |
| `B_base` | 0.05, dimensionless, constant. | Assumed 5% blockage multiplier on the final river outlet only. |
| `gate1_hours`, `gate1_values` | Times `[0,17,18,26,28,50,65,66,76,80,110,144]` h; openings `[.35,.35,.20,.20,.50,.35,.35,.20,.20,.50,.35,.35]`. | Assumed operational schedule, linearly interpolated between knots; not optimized or observed. |
| `gate2_hours`, `gate2_values` | Times `[0,18,20,27,30,52,66,68,78,82,112,144]` h; openings `[.45,.45,.30,.30,.60,.45,.45,.30,.30,.60,.45,.45]`. | Independent assumed schedule, linearly interpolated; same interpretation as Gate 1. |

### Numerical settings, geometry checks, and reporting values

| Code field | Accepted value | Provenance and meaning |
|---|---|---|
| `hour`, `day` | 3,600 s; 86,400 s | Unit conversions (`day=24×hour`). |
| `duration` | 518,400 s = 144 h = 6 days | Numerical experiment window chosen for the two synthetic storms and recession. |
| `time_step` | 120 s | Numerical integration/input-sampling interval. The `time` property generates 4,321 nodes including both endpoints. |
| `H_flood` | 97.90 m | Illustrative reporting threshold for time above level, not a surveyed flood elevation or engineering safety standard. |
| `H_r_valid_max` | 99.0 m | Assumed river geometry-use limit; the diagnostic rejects a violating trajectory rather than clipping storage. |
| `H_g_valid_max` | 97.9 m | Assumed groundwater separation limit, equal to the lower wetland bed; admissible `H_g` is strictly below it. |
| `soil_bottom` | 100.0 m | Assumed soil-region lower reference checked against internal groundwater; not a borehole observation. |
| `balance_rtol` | 1 × 10⁻⁴, dimensionless | Numerical tolerance on maximum cumulative mass-balance discrepancy divided by gross event inflow (or initial total storage if no inflow). Not a calibration or physical-accuracy tolerance. |

The official `ode4u` implementation and pinned revision are documented in [the installation section](README.md#baseline-ode-simulation). This identifies the numerical method, not a source for the assumed installation parameters.

### Additional values used only in the A–G algebraic illustrations

These settings come from [algebraic_plots.py](algebraic_plots.py). They illustrate the equations and are neither additional baseline inputs nor empirical uncertainty ranges.

| Figure | Extra comparison or sampling values | Provenance |
|---|---|---|
| A | Effective head −0.2 to 3.0 m (401 points); openings 0.25, 0.50, 0.75, 1.00. | Chosen display domain and comparison openings; Branch 1 `R_G` and `Q_Gcap` from the tables above. |
| B | Effective head −0.2 to 3.0 m (401 points); bare cover 0 and class A/B cover 0.65. | Chosen comparison designs. Derived Branch 1 resistances 0.8, 2.10, 2.88 s/m²; shared cap 0.65 m³/s. |
| C | Overflow head −0.2 to 1.0 m (401 points); coefficient multipliers 0.5, 1.0, 1.5. | Chosen sensitivity illustrations around Branch 1 `c_wsp`, giving 1.5, 3.0, 4.5 m^(3/2)/s. |
| D | Soil fraction 0–1 and surface storage 0–`10 V_u`=6,000 m³ (301 points each); fixed `S_u=V_u` or soil fraction 0.50; `K_u` multipliers 0.5, 1, 2. | Chosen display settings; rates 0.5, 1, 2 × 10⁻⁶ m/s are not local measured ranges. |
| E | Temperature 5–35 °C (301 points); radiation 100, 300, 600 W/m²; wind 2 m/s. | Chosen evaporative-demand illustration; these radiation levels differ from the synthetic baseline peak by design. |
| F | Effective head −0.2 to 2.0 m (401 points); blockage 0, 0.25, 0.50, 0.75. | Chosen outlet comparisons; effective caps derived as 3, 2.25, 1.5, 0.75 m³/s. |
| G | External-minus-internal head −1 to 1 m (401 points); `C_gx` multipliers 0.5, 1, 2; initial `S_g` held fixed. | Chosen groundwater comparisons; conductances 0.0175, 0.035, 0.070 m²/s. |

### Recorded baseline results and calculation lineage

These values are **simulation outputs**, not observations. They were read from the existing local `A/results/baseline_summary.json`, generated by [baseline_diagnostics()](wetland_flood_model.py#L410), without rerunning or changing the model. The ignored result file need not be published for the recorded values and their definitions to be visible here; the source and run commands reproduce it.

| Quantity | Recorded value | How calculated |
|---|---:|---|
| Peak river stage | 97.87940653911606 m | Maximum saved `H_r = z_r + S_r/A_r`. |
| Peak outlet flow | 1.3578956737526984 m³/s | Maximum saved `Q_out`; this is not the same quantity as peak stage. |
| Total outlet volume | 171,870.64533543994 m³ | Trapezoidal time integral of `Q_out`. |
| Time above the illustrative 97.90 m threshold | 0 h | Piecewise-linear crossing calculation on the saved river-stage series. |
| Gross external inflow | 184,872 m³ | Trapezoidal integral of `P A_rain + Q_ext + Q_gx,in`, excluding all internal transfers. |
| Maximum cumulative water-balance residual | 0.07673464775143657 m³ | Maximum absolute `S_tot(t) − S_tot(0) − integral(external_net)`. |
| Final cumulative water-balance residual | 0.022189399605622384 m³ | Same expression at the final sample, signed. |
| Maximum relative balance residual | 4.1506906265652217 × 10⁻⁷ | Maximum absolute residual / gross event inflow; approximately 0.00004151%. |
| Maximum instantaneous derivative-balance discrepancy | 8.881784197001252 × 10⁻¹⁶ m³/s | Maximum absolute `sum(dxdt) − external_net`, showing algebraic cancellation to floating-point precision. |
| Maximum groundwater level | 97.4312180648551 m | Maximum simulated internal `H_g`; no measured groundwater record. |
| Maximum soil-capacity fraction | 0.7500005892726511 | Maximum `S_m/C_m`; no observed soil-moisture record. |

The minimum saved storages, in state order, are approximately `[3.900596e−38, 71954.401183, 202.940488, 104.274471, 523.718290, 390.169397, 409600.0, 3914.587143]` m³. They are the channel minima reported by the same diagnostic, not simultaneous states or imposed lower bounds. Numerical consistency and finite positive states do not establish physical calibration, engineering safety, or an optimization result.
