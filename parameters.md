# Parameter reference — accepted eight-state Stage A model

This reference supports the concise `model_notes.md` and the current `wetland_flood_model.py`. It lists the existing values; it does not introduce new equations or a parameter-estimation exercise. Branch arrays are ordered **[1, 2]**. Code names below are attributes of an instance of `class c`, unless otherwise stated.

**Value categories:** “illustrative geometry,” “illustrative empirical,” “illustrative model,” “illustrative initial,” “synthetic disturbance,” and “illustrative decision” all mean **uncalibrated assumptions**. “Physical approximation” identifies constants already used by accepted v3. Mathematical admissibility domains are not defensible site ranges. **Evidence-supported site-specific numerical ranges remain pending**: the supplied materials establish none. Admissibility domains alone do not satisfy that requirement. The last column identifies possible future data, without collecting it in Stage A.

## Geometry and reporting levels

| Symbol | Python name | Definition | Units | Current value | Value category | Admissible domain; uncertainty and future data |
|---|---|---|---|---|---|---|
| A_u | `A_u` | Upland footprint | m² | 600000 | Illustrative geometry | >0; uncertain; mapped contributing area. |
| A_d1,A_d2 | `A_d` | Dam-pool effective areas | m² | [10000,8000] | Illustrative geometry | >0; uncertain; pool surveys. |
| A_w1,A_w2 | `A_w` | Wetland effective areas | m² | [12000,10000] | Illustrative geometry | >0; uncertain; wetland surveys. |
| A_r | `A_r` | River effective storage/rainfall area | m² | 20000 | Illustrative geometry | >0; uncertain; reach geometry. |
| A_rain | `A_rain` property | Sum of six disjoint surface footprints | m² | 660000 | Derived geometry | No double counting; depends on footprint assumptions. |
| A_g | `A_g` | Effective aquifer area, excluding river footprint | m² | 640000 | Illustrative geometry | >0; uncertain; hydrogeologic mapping. |
| S_y | `S_y` | Aquifer specific yield | 1 | 0.20 | Illustrative empirical | 0<S_y≤1; uncertain; aquifer observations. |
| z_d1,z_d2 | `z_d` | Pool zero-storage elevations | m | [100.0,99.8] | Illustrative geometry | Common datum; uncertain; survey. |
| z_w1,z_w2 | `z_w` | Wetland zero-storage elevations | m | [98.0,97.9] | Illustrative geometry | Common datum; uncertain; survey. |
| z_r | `z_r` | River zero-storage elevation | m | 97.0 | Illustrative geometry | Common datum; uncertain; survey. |
| z_g | `z_g` | Lower reference for drainable groundwater | m | 94.0 | Illustrative geometry | Common datum; uncertain; aquifer interpretation. |
| H_flood | `H_flood` | River flood-reporting threshold | m | 97.90 | Illustrative model | Reporting only, not a flux or cap; needs site definition. |
| H_r,max | `H_r_valid_max` | Assumed river geometry-use limit | m | 99.0 | Illustrative model | Diagnostic limit, not clipping; needs geometry evidence. |
| H_g,max | `H_g_valid_max` | Conservative groundwater overlap limit | m | 97.9 | Illustrative model | H_g<H_g_valid_max; separate bed/soil checks also apply; needs hydrogeologic evidence. |
| z_soil,bottom | `soil_bottom` | Bottom of separate upland soil region | m | 100.0 | Illustrative geometry | H_g<soil_bottom; not another state or flow. |

The surface footprints are disjoint. Soil and groundwater receive infiltration/recharge rather than another rainfall input. The aquifer may underlie other footprints while representing a different vertical water-storage region. All levels use the same datum. Initial and computed groundwater levels must satisfy **H_g<H_g_valid_max, H_g<soil_bottom, and H_g<min(z_w)** separately. Relevant limits/elevations must be finite. These are applicability assumptions of this configuration, not universal hydraulic laws; a violation requires reviewing geometry, not clipping or assuming a smaller timestep fixes it.

**Editing order:** edit geometry → explicitly call `params.refresh_geometry_defaults()` when reapplying the existing defaults is intended → construct initial conditions and inputs → validate → run. The helper deliberately overwrites only `A_g=A_u+sum(A_d)+sum(A_w)`, `V_w=.01*A_w`, `V_d=.01*A_d`, `V_g=.05*S_y*A_g`, and `V_r=.01*A_r`. Simulation and plotting do not call it automatically. Independently supplied linked values can instead be retained by omitting the refresh or assigned after it. `C_m`, `S_fc`, `V_u`, `V_m`, conductances, ratings, and ET coefficients are not rescaled. To preserve the illustrative upland depth assumptions after changing A_u, manually set `C_m=.20*A_u`, `S_fc=.10*A_u`, `V_u=.001*A_u`, and `V_m=.01*A_u` if those assumptions remain appropriate.

## Soil, upland routing, and availability

| Symbol | Python name | Definition | Units | Current value | Value category | Admissible domain; uncertainty and future data |
|---|---|---|---|---|---|---|
| C_m | `C_m` | Maximum soil water storage | m³ | 120000 | Illustrative model | >0; uncertain; soil depth/storage data; 0.20 m water equivalent over A_u. |
| S_fc | `S_fc` | Soil drainage threshold | m³ | 60000 | Illustrative empirical | 0≤S_fc≤C_m; uncertain; soil retention/drainage data. |
| K_u | `K_u` | Effective infiltration coefficient | m/s | 1e−6 (3.6 mm/h) | Illustrative empirical | ≥0; uncertain; infiltration measurements. |
| k_mg | `k_mg` | Soil-to-groundwater drainage rate | s⁻¹ | 1/(5×86400) | Illustrative empirical | >0 in this implementation; uncertain; recession data. |
| k_ud1,k_ud2 | `k_ud` | Upland routing rates to pools | s⁻¹ | [1/(8×3600),1/(10×3600)] | Illustrative empirical | >0; uncertain; runoff timing. |
| k_ur | `k_ur` | Direct upland-to-river routing rate | s⁻¹ | 1/(20×3600) | Illustrative empirical | >0; uncertain; bypass runoff timing. |
| Q̄_ud1,Q̄_ud2 | `Q_udcap` | Pool-inflow routing caps | m³/s | [1.2,1.0] | Illustrative empirical | >0; uncertain; conveyance geometry/ratings. |
| Q̄_ur | `Q_urcap` | Direct-runoff cap | m³/s | 0.7 | Illustrative empirical | >0; uncertain; bypass rating. |
| V_u* | `V_u` | Upland half-availability storage | m³ | 600 | Illustrative empirical | >0; uncertain; equals 1 mm over A_u. |
| V_m* | `V_m` | Soil half-availability storage | m³ | 6000 | Illustrative empirical | >0; uncertain; equals 10 mm over A_u. |
| V_w1*,V_w2* | `V_w` | Wetland half-availability storages | m³ | [120,100] | Illustrative empirical | >0; uncertain; 1 cm over each wetland. |
| V_d1*,V_d2* | `V_d` | Pool half-availability storages | m³ | [100,80] | Illustrative empirical | >0; uncertain; 1 cm over each pool. |
| V_g* | `V_g` | Groundwater half-availability storage | m³ | 6400 | Illustrative empirical | >0; uncertain; 0.05 m head×S_y A_g. |
| V_r* | `V_r` | River half-availability storage | m³ | 200 | Illustrative empirical | >0; uncertain; 1 cm over A_r. |

Availability is phi_j=S_j/(S_j+V_j*), not a numerical clipping rule. The admissible state domain is S_j≥0, with S_m≤C_m. Assigned storage scales need future evidence before site use.

## Gates, wetlands, groundwater, and river ratings

| Symbol | Python name | Definition | Units | Current value | Value category | Admissible domain; uncertainty and future data |
|---|---|---|---|---|---|---|
| z_G1,z_G2 | `z_G` | Gate threshold elevations | m | [100.05,99.85] | Illustrative geometry | z_G≥z_d; uncertain; gate geometry. |
| R_G1,R_G2 | `R_G` | Full-opening gate resistances | s/m² | [1.6,1.8] | Illustrative empirical | >0; uncertain; head-discharge measurements. |
| Q̄_G1,Q̄_G2 | `Q_Gcap` | Full-opening gate flow caps | m³/s | [1.0,0.8] | Illustrative empirical | >0; actual cap is g_i Q̄_Gi; needs rating data. |
| z_dsp1,z_dsp2 | `z_dsp` | Pool spill crests | m | [100.9,100.7] | Illustrative geometry | z_dsp>z_d and z_dsp>z_G in each branch; configuration assumption; uncertain; survey. |
| c_dsp1,c_dsp2 | `c_dsp` | Dam spill coefficients | m^(3/2)/s | [4.0,3.2] | Illustrative empirical | >0; uncertain; spillway geometry/rating. |
| z_w1,o,z_w2,o | `z_wo` | Normal wetland outlet thresholds | m | [98.05,97.95] | Illustrative geometry | z_wo≥z_w; uncertain; survey. |
| Q̄_wr1,Q̄_wr2 | `Q_wrcap` | Normal wetland drainage caps | m³/s | [0.65,0.55] | Illustrative empirical | >0; uncertain; outlet ratings. |
| z_wsp1,z_wsp2 | `z_wsp` | Wetland overflow thresholds | m | [98.6,98.5] | Illustrative geometry | z_wsp>z_w and z_wsp>z_wo in each branch; configuration assumption; uncertain; survey. |
| c_wsp1,c_wsp2 | `c_wsp` | Wetland overflow coefficients | m^(3/2)/s | [3.0,2.5] | Illustrative empirical | >0; uncertain; overflow-path rating. |
| C_d1,C_d2 | `C_d` property | Pool spill-onset storages | m³ | [9000,7200] | Derived geometry | A_d(z_dsp−z_d); thresholds, not hard caps. |
| C_w1,C_w2 | `C_w` property | Wetland overflow-onset storages | m³ | [7200,6000] | Derived geometry | A_w(z_wsp−z_w); thresholds, not hard caps. |
| C_wg1,C_wg2 | `C_wg` | Wetland–groundwater conductances | m²/s | [0.015,0.012] | Illustrative empirical | ≥0; uncertain; seepage/head data. |
| C_gr | `C_gr` | Groundwater-to-river conductance | m²/s | 0.08 | Illustrative empirical | ≥0; uncertain; groundwater/baseflow data. |
| C_gx | `C_gx` | External groundwater conductance | m²/s | 0.035 | Illustrative empirical | ≥0; uncertain; boundary-head/flux data. |
| z_r,o | `z_ro` | Final river outlet threshold | m | 97.1 | Illustrative geometry | z_ro≥z_r; uncertain; outlet survey. |
| R_r | `R_r` | Final river outlet resistance | s/m² | 0.45 | Illustrative empirical | >0; uncertain; reach/outlet rating. |
| Q̄_r | `Q_rcap` | Unblocked river outlet cap | m³/s | 3.0 | Illustrative empirical | >0; actual cap (1−B)Q̄_r; needs rating data. |

## Wetland design choices and endpoint assumptions

| Symbol | Python name | Definition | Units | Current value | Value category | Admissible domain; uncertainty and future data |
|---|---|---|---|---|---|---|
| τ_1,τ_2 | `tau` | Vegetation categories | category | [type_a,type_a] | Illustrative decision | Choose a defined category; synthetic labels, not species identities. |
| ρ_1,ρ_2 | `rho` | Fractional vegetation cover | 1 | [0.65,0.60] | Illustrative decision | [0,1]; prescribed, not random or a growth state; future feasible designs. |
| R_bare,1,R_bare,2 | `R_bare` | Zero-cover resistance endpoints | s/m² | [0.8,1.0] | Illustrative empirical | >0; uncertain; drainage experiments. |
| R_veg,1,τ,R_veg,2,τ | `R_veg` | Full-cover resistance endpoints | s/m² | type_a: [2.8,3.0]; type_b: [4.0,4.5] | Illustrative empirical | >0; uncertain; actual vegetation/flow data. |
| K_bare,1,K_bare,2 | `K_bare` | Zero-cover ET endpoints | 1 | [1.0,1.0] | Illustrative empirical | ≥0; uncertain; ET evidence compatible with E_star. |
| K_veg,1,τ,K_veg,2,τ | `K_veg` | Full-cover ET endpoints | 1 | type_a: [1.1,1.1]; type_b: [1.25,1.25] | Illustrative empirical | ≥0; uncertain; vegetation/ET observations. |
| R_w1,R_w2 | First return from `design()` | Interpolated wetland resistance | s/m² | [2.1,2.2] | Derived design property | Positive from endpoints; inherits endpoint uncertainty. |
| K_ET,w1,K_ET,w2 | Second return from `design()` | Interpolated wetland ET coefficient | 1 | [1.065,1.06] | Derived design property | Nonnegative from endpoints; inherits endpoint uncertainty. |

`design()` uses (1−ρ)×bare+ρ×selected category endpoint independently for each wetland. `type_a` and `type_b` encode illustrative endpoint pairs only. Their ordering is not a universal vegetation law. Design changes normal drainage and wetland ET in this model; it does not change seepage, geometry, or overflow coefficients automatically. Designs remain fixed throughout a run.

## Reduced evaporation relationship

| Symbol | Python name | Definition | Units | Current value | Value category | Admissible domain; uncertainty and future data |
|---|---|---|---|---|---|---|
| λ | `lambda_v` | Latent heat of vaporization | J/kg | 2.45e6 | Physical approximation | >0; v3 constant approximation near 20 °C, not a site fit. |
| c_p | `c_p` | Air specific heat | J/(kg·K) | 1013 | Physical approximation | >0; accepted v3 approximation. |
| ε | `epsilon` | Molecular-mass ratio | 1 | 0.622 | Physical approximation | >0; accepted v3 constant. |
| ρ_w | `rho_w` | Liquid-water density | kg/m³ | 1000 | Physical approximation | >0; rounded value. |
| R_v | `R_v` | Water-vapor gas constant | J/(kg·K) | 461.5 | Physical approximation | >0; accepted v3 approximation. |
| p_0 | `p_0` | Fixed atmospheric pressure | Pa | 101325 | Illustrative model | >0; future site pressure/elevation. |
| r_H | `r_H` | Fixed effective relative humidity | 1 | 0.65 | Illustrative empirical | [0,1]; uncertain; humidity data, not a ninth disturbance here. |
| η | `eta` | Available-radiation fraction | 1 | 0.60 | Illustrative empirical | [0,1]; uncertain; radiation/energy evidence. |
| k_0 | `k_0` | Still-air transfer coefficient | m/s | 0.001 | Illustrative empirical | ≥0; uncertain; exchange/ET evidence. |
| k_v | `k_v` | Wind multiplier in air transfer | 1 | 0.002 | Illustrative empirical | ≥0; uncertain; wind-height/exchange evidence. |
| K_ET,m | `K_ETm` | Soil ET coefficient | 1 | 0.85 | Illustrative empirical | ≥0; uncertain; compatible ET observations. |
| K_ET,d1,K_ET,d2 | `K_ETd` | Pool evaporation coefficients | 1 | [1.0,1.0] | Illustrative empirical | ≥0; uncertain; evaporation observations. |
| K_ET,r | `K_ETr` | River evaporation coefficient | 1 | 1.0 | Illustrative empirical | ≥0; uncertain; evaporation observations. |
| e_s formula constants | Literals in `evaporative_demand()` | Saturation-pressure relationship retained from v3 | Pa; 1; °C | 610.8; 17.27; 237.3 | Accepted formula coefficients | Used as specified in v3; not new estimated parameters. |
| Celsius-to-Kelvin offset | Literal in `evaporative_demand()` | T_K=T_a+273.15 | K offset | 273.15 | Unit conversion | Exact temperature-scale conversion. |

The accepted closure is reduced wet-surface evaporative demand, not standardized FAO reference ET0. Changing temperature, radiation, or wind changes E_star; actual withdrawals also multiply source area, coefficient, and source availability. The liquid-water demonstration uses T_a>0 °C. Freezing, snow, condensation, variable humidity, and thermal storage are outside this setup.

## Initial conditions

| Symbol | Python name | Definition | Units | Current value | Value category | Admissible domain; uncertainty and future data |
|---|---|---|---|---|---|---|
| d_u(0) | `u_depth_init` | Initial upland surface depth | m | 0.001 | Illustrative initial | ≥0; uncertain; antecedent wetness. |
| S_m(0)/C_m | `soil_fraction_init` | Initial soil wetness fraction | 1 | 0.75 | Illustrative initial | [0,1]; uncertain; soil moisture. |
| d_w1(0),d_w2(0) | `w_depth_init` | Initial wetland depths | m | [0.25,0.25] | Illustrative initial | ≥0; uncertain; stage observations. |
| d_d1(0),d_d2(0) | `d_depth_init` | Initial pool depths | m | [0.30,0.30] | Illustrative initial | ≥0; uncertain; stage observations. |
| d_r(0) | `r_depth_init` | Initial river depth | m | 0.35 | Illustrative initial | ≥0; uncertain; stage observations. |
| H_g(0) | `H_g_init` | Initial internal groundwater level | m | 97.20 | Illustrative initial | ≥z_g and below all three geometry limits above; uncertain; piezometer data. |
| x(0) | `x_init` property | Storage vector in accepted state order | m³ | [600,90000,3000,2500,3000,2400,409600,7000] | Derived initial | Nonnegative; S_m≤C_m; inherits preceding uncertainty. |

Initial levels are H_w=[98.25,98.15], H_d=[100.30,100.10], H_r=97.35, H_g=97.20 m. Groundwater therefore begins below both wetlands and the river.

## Prescribed forcing and controllable gate histories

Let h=t/3600. A tuple `(a,p,b,A)` means a triangular pulse: zero at hour a, amplitude A at p, zero at b, and zero outside. Between knots it is linear. All eight disturbances are prescribed independently of the ODE states; the outside-watershed inflow is not recalculated from precipitation over the modeled footprints. These synthetic inputs need not imply statistical independence in a real watershed.

| Symbol | Python name | Definition | Units | Current value/history | Value category | Admissible domain; uncertainty and future data |
|---|---|---|---|---|---|---|
| Weather/boundary timing | `weather_period_h`, `temperature_phase_h`, `sunrise_phase_h`, `wind_period_h`, `groundwater_period_h` | Stored periods and phase offsets in the formulas below | h | 24; 9; 6; 36; 144 | Synthetic disturbance settings | Periods >0; phase offsets define timing; illustrative, not observed. |
| P(t) | `rain_pulses` | Rain triangles `(start,peak,end,amplitude)` | h; mm/h | [(18,22,26,15),(66,71,76,20)] | Synthetic disturbance | P≥0; illustrative event, no claimed return period; future rainfall. |
| T_a(t) | `T_mean`, `T_amplitude` | Mean and daily amplitude | °C | 22+4 sin(2π(h−9)/24) | Synthetic disturbance | This demonstration >0; future meteorological observations. |
| R_s(t) | `R_s_peak` | Daily daylight radiation | W/m² | 550 max(0,sin(2π(h−6)/24)) | Synthetic disturbance | ≥0; future radiation observations. |
| v(t) | `wind_mean`, `wind_amplitude` | Mean and 36-hour amplitude | m/s | 2+0.6 sin(2πh/36) | Synthetic disturbance | ≥0; future wind observations. |
| Q_ext(t) | `Q_ext_base`, `ext_pulses` | Independent outside-watershed river inflow | m³/s; h | Base 0.08; pulses (20,25,36,0.35), (68,76,90,0.70) | Synthetic disturbance | ≥0 inward; future boundary hydrograph; exclude modeled runoff. |
| h_d(t) | `h_d_base`, `h_d_pulse` | Downstream outlet boundary level | m; h | Base 97.15; pulse (66,77,94,0.10) | Synthetic disturbance | Common datum; future boundary stages. |
| h_g(t) | `h_g_base`, `h_g_amplitude` | External groundwater boundary level | m | 97.15+0.025 sin(2πh/144) | Synthetic disturbance | Common datum, distinct from H_g; future groundwater monitoring. |
| B(t) | `B_base` | Final-outlet blockage factor | 1 | 0.05 | Synthetic disturbance | [0,1]; uncertain; outlet condition evidence. |
| g_1(t) | `gate1_hours`, `gate1_values` | Gate 1 knot interpolation | h; 1 | Hours [0,17,18,26,28,50,65,66,76,80,110,144]; values [.35,.35,.20,.20,.50,.35,.35,.20,.20,.50,.35,.35] | Illustrative decision | [0,1]; prescribed operating choice; future feasible operation rules. |
| g_2(t) | `gate2_hours`, `gate2_values` | Gate 2 knot interpolation | h; 1 | Hours [0,18,20,27,30,52,66,68,78,82,112,144]; values [.45,.45,.30,.30,.60,.45,.45,.30,.30,.60,.45,.45] | Illustrative decision | [0,1]; independently prescribed operating choice. |

Rainfall conversion is P[m/s]=P[mm/h]×1e−3/3600. The two storms have depths 60 and 100 mm. Gate values are linearly interpolated between listed knots. The solver linearly interpolates all sampled numerical-input channels between time samples, including sampled smooth weather signals.

## Simulation settings and indexing reference

| Symbol | Python name | Definition | Units | Current value | Value category | Admissible domain; uncertainty and future data |
|---|---|---|---|---|---|---|
| hour; day | `hour`, `day` | Unit conversions | s | 3600; 86400 | Unit conversion | Exact. |
| T | `duration` | Baseline event window | s | 518400 (6 days) | Demonstration setting | >0; chosen teaching example. |
| Δt | `time_step` | Fixed integration/input-sampling step | s | 120 | Numerical setting | >0; duration is an integer number of steps. |
| ε_balance | `balance_rtol` | Relative water-balance reporting criterion | 1 | 1e−4 | Numerical setting | Illustrative numerical criterion, not physical accuracy. |

The baseline uses 4,321 time samples. It does not automatically run a scenario suite or timestep comparison. Exact channel names are defined by `STATE_NAMES`, `INPUT_NAMES`, and `OUTPUT_NAMES` in the main script. Array shape is **(channels, time)**; all indices below are zero-based.

- States: `0 S_u, 1 S_m, 2 S_w1, 3 S_w2, 4 S_d1, 5 S_d2, 6 S_g, 7 S_r`.
- Inputs: `0 g1, 1 g2, 2 P, 3 T_a, 4 R_s, 5 v, 6 Q_ext, 7 h_d, 8 h_g, 9 B`.
- Outputs: `0 H_w1, 1 H_w2, 2 H_d1, 3 H_d2, 4 H_g, 5 H_r, 6 Q_G1, 7 Q_G2, 8 Q_wr1, 9 Q_wr2, 10 Q_dsp1, 11 Q_dsp2, 12 Q_wsp1, 13 Q_wsp2, 14 Q_inf, 15 Q_mg, 16 Q_wg1, 17 Q_wg2, 18 Q_gr, 19 Q_gx, 20 Q_out, 21 E_tot, 22 S_tot, 23 E_star, 24 Q_gx_in, 25 Q_gx_out, 26 Q_ud1, 27 Q_ud2, 28 Q_ur, 29 E_m, 30 E_w1, 31 E_w2, 32 E_d1, 33 E_d2, 34 E_r, 35 P_total, 36 external_net`.

State units are m³. Output indices 0–5 are m; 22 is m³; 23 is m/s; all other outputs are m³/s. Input units are in the forcing table. Algebraic-plot comparison values are local teaching settings in `algebraic_plots.py`, not extra baseline parameters or new whole-system scenarios.
