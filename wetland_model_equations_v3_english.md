# Wetland Restoration and Flood Mitigation: Eight-State Equation Draft v3

Date: 2026-09-18

**Status: equation-review draft; no simulation code has been written or executed.**

## 1. Accepted structure

Based on the latest user specification, subsequent topology corrections, and the concept–equations–code workflow in `Project-ODE-how-to.md`. The earlier placement of wetlands upstream of the gates has been superseded.

```text
             +--> Dam pool 1 -- Gate 1 --> Wetland 1 --+
             |                                        |
Upland ------+--> Dam pool 2 -- Gate 2 --> Wetland 2 --+--> River --> Outside
             |                                        |
             +--------------- direct runoff ----------+
             |
             +--> Soil --> Groundwater --> River

Wetland 1 --> Groundwater
Wetland 2 --> Groundwater
External groundwater <--> Groundwater
External inflow --> River

Dam-pool spillway i --> Wetland i   (bypasses Gate i)
Wetland overflow i --> River
```

- Surface routing is one-way for now; high tailwater can still suppress drainage.
- Blockage B acts only on the final river outlet.
- Dam spills enter the associated wetlands; wetland overflows enter the river.
- The direct upland-to-river bypass, soil state, and groundwater state are retained.
- The no-reverse-flow choice applies to the discussed surface connections. The original signed external groundwater-boundary exchange is retained separately.

**Source distinction:** The topology and decision categories come from the user; the eight balance equations follow conservation; the stage, flow-rating, ET, and design-interpolation laws below are proposed closures for review.

## 2. Categories, ordering, and units

Time t is in seconds, storage in m³, and volumetric fluxes in m³/s. All levels share one elevation datum.

\[
x=[S_u,S_m,S_{w1},S_{w2},S_{d1},S_{d2},S_g,S_r]^T.
\]

| Index | State | Definition | Units |
|---|---|---|---|
| 1 | S_u | Upland surface water | m³ |
| 2 | S_m | Water in the lumped upland soil compartment | m³ |
| 3 | S_w1 | Wetland 1 surface water | m³ |
| 4 | S_w2 | Wetland 2 surface water | m³ |
| 5 | S_d1 | Dam pool 1 water | m³ |
| 6 | S_d2 | Dam pool 2 water | m³ |
| 7 | S_g | Drainable groundwater above a fixed lower reference | m³ |
| 8 | S_r | Water in the modeled downstream reach | m³ |

Operational controls:
\[
u_c(t)=[g_1(t),g_2(t)]^T,\qquad 0\le g_i\le1.
\]
g_i is normalized opening: 0 closed and 1 fully open.

Designs fixed within a simulation:
\[
D_i=(\tau_i,\rho_i),\quad i=1,2,\qquad 0\le\rho_i\le1.
\]
τ_i is categorical vegetation type, not a continuous variable. ρ_i is fractional vegetation cover used as a density descriptor; it is not stem density or biomass. Four conceptual decision groups do not imply four numerical scalars.

External disturbances:
\[
w(t)=[P,T_a,R_s,v,Q_{\rm ext},h_d,h_g,B]^T.
\]

| Disturbance | Units and sign | Explicit entry |
|---|---|---|
| P | m/s, P ≥ 0 | P A_j: surface precipitation inputs |
| T_a | °C | Vapor pressure, ET weights, vapor density |
| R_s | W/m², ≥ 0 | Effective radiative energy: η R_s |
| v | m/s, ≥ 0 | Air-transfer conductance: k_0 + k_v v |
| Q_ext | m³/s, ≥ 0 inward | Direct inflow to S_r |
| h_d | m | River-outlet tailwater |
| h_g | m | Groundwater-boundary exchange |
| B | dimensionless, 0–1 | Multiplier 1−B on Q_out |

Parameters and initial conditions are not states; levels and fluxes are algebraic outputs; peaks are postprocessed metrics. Future code storage conventions do not change these categories.

## 3. Common algebraic definitions

### 3.1 Positive part and availability

\[
[z]_+=\max(z,0),\qquad
\phi_j(S_j)=\frac{S_j}{S_j+V_j^*},\quad S_j\ge0,\quad V_j^*>0.
\]

V_j* is the storage in m³ at which the availability factor reaches one half. The factor limits infiltration, seepage, groundwater withdrawals, or ET near dry conditions. It is a constitutive storage scale, not a timestep or a numerical state-clipping device. Different compartments can have different values.

### 3.2 Storage–stage relationships

The initial approximation uses fixed areas:
\[
H_j=z_j+\frac{S_j}{A_j},\qquad j\in\{d1,d2,w1,w2,r\},
\]
\[
H_g=z_g+\frac{S_g}{S_y A_g}.
\]

A_j is effective horizontal area and z_j is the zero-storage level. A_g and S_y are aquifer area and specific yield. H_g is internal groundwater level, distinct from external h_g. Soil wetness is represented by S_m/C_m without another state.

### 3.3 Rainfall areas

\[
P_j=P(t)A_j,\quad j\in\{u,d1,d2,w1,w2,r\}.
\]

The six surface footprints must be disjoint and together cover the modeled area whose rainfall is included here. Soil and aquifer areas are not added again to rainfall area. Under the fixed-area approximation, precipitation over each assigned footprint enters that compartment.

\[
P\,[\mathrm{m/s}]=P\,[\mathrm{mm/h}]\frac{10^{-3}}{3600}.
\]

Q_ext excludes runoff already generated within these modeled footprints.

## 4. All eight water-balance ODEs

Every Q, E, and P_j below is in m³/s and is defined in Sections 5–6. Each internal transfer uses exactly the same flux at its two ends.

\[
\boxed{\dot S_u=P_u-Q_{\rm inf}-Q_{ud,1}-Q_{ud,2}-Q_{ur}} \tag{1}
\]
\[
\boxed{\dot S_m=Q_{\rm inf}-Q_{mg}-E_m} \tag{2}
\]
\[
\boxed{\dot S_{w1}=P_{w1}+Q_{G1}+Q_{dsp,1}-Q_{wr,1}-Q_{wsp,1}-Q_{wg,1}-E_{w1}} \tag{3}
\]
\[
\boxed{\dot S_{w2}=P_{w2}+Q_{G2}+Q_{dsp,2}-Q_{wr,2}-Q_{wsp,2}-Q_{wg,2}-E_{w2}} \tag{4}
\]
\[
\boxed{\dot S_{d1}=P_{d1}+Q_{ud,1}-Q_{G1}-Q_{dsp,1}-E_{d1}} \tag{5}
\]
\[
\boxed{\dot S_{d2}=P_{d2}+Q_{ud,2}-Q_{G2}-Q_{dsp,2}-E_{d2}} \tag{6}
\]
\[
\boxed{\dot S_g=Q_{mg}+Q_{wg,1}+Q_{wg,2}+Q_{gx}-Q_{gr}} \tag{7}
\]
\[
\boxed{\dot S_r=P_r+Q_{ur}+Q_{wr,1}+Q_{wr,2}+Q_{wsp,1}+Q_{wsp,2}+Q_{gr}+Q_{\rm ext}-Q_{\rm out}-E_r} \tag{8}
\]

**Interpretation of each ODE**

| Eq. | Gains and losses | Controls, disturbances, interactions |
|---|---|---|
| 1: Upland | Rain adds water; infiltration and three runoff routes remove it. | P acts directly; soil saturation suppresses infiltration; pools and bypass receive runoff. |
| 2: Soil | Infiltration adds water; recharge and ET remove it. | Connects upland to groundwater; weather acts through E_m. |
| 3: Wetland 1 | Rain, Gate 1, and dam spill add water; drainage, overflow, seepage, and ET remove it. | g_1 affects inflow; D_1 affects drainage and ET; river and groundwater heads affect drainage. |
| 4: Wetland 2 | Same categories as Wetland 1, with independent fluxes and parameters. | Uses g_2, D_2, and branch-2 parameters. |
| 5: Pool 1 | Rain and runoff add water; gate release, spill, and evaporation remove it. | g_1 controls normal release; H_w1 is tailwater. |
| 6: Pool 2 | Same categories as Pool 1, with independent storage. | g_2 controls release; H_w2 is tailwater. |
| 7: Groundwater | Recharge and boundary inflow add water; boundary outflow and river discharge remove it. | h_g acts through Q_gx; H_r affects Q_gr. |
| 8: River | Internal arrivals, external inflow, and local rain add water; outlet flow and evaporation remove it. | h_d and B affect the outlet; gates and designs act through upstream routing. |

Initial simplifications omit evaporation from transient upland ponding and seepage through dam-pool beds, while retaining pool and river evaporation. Wetland ET is withdrawn from S_wi only. Root uptake from subsurface water after a wetland dries is not represented.

## 5. Explicit transfer relationships

### 5.1 Runoff, infiltration, and recharge

\[
Q_{ud,i}=\min(\bar Q_{ud,i},k_{ud,i}S_u),\quad i=1,2,
\]
\[
Q_{ur}=\min(\bar Q_{ur},k_{ur}S_u),
\]
\[
Q_{\rm inf}=A_uK_u\phi_u(S_u)\left[1-\frac{S_m}{C_m}\right]_+,
\]
\[
Q_{mg}=k_{mg}[S_m-S_{fc}]_+.
\]

The k coefficients are in s⁻¹, flow caps in m³/s, and K_u in m/s. C_m is maximum soil-water storage, and S_fc is the gravitational-drainage threshold, with 0 ≤ S_fc ≤ C_m. Both pool inflows are withdrawn from S_u. Water rejected by saturated soil remains in S_u.

These are empirical lumped closures, not a Richards-equation model. The soil and aquifer must remain distinguishable storage regions in the modeled regime; groundwater intrusion into the soil compartment requires further treatment.

### 5.2 One-way head–discharge rating

Define an explicit helper relationship:
\[
\mathcal F(H_a,H_b,z_o;\mathcal R,\bar Q)
=\min\left(\bar Q,\frac{[H_a-\max(H_b,z_o)]_+}{\mathcal R}\right).
\]

H_a and H_b are donor and receiver levels; z_o is the outlet threshold elevation. R has units s/m² and Q̄ is rated maximum discharge. The threshold must be at or above the donor's zero-storage elevation. Unfavorable head produces zero flow, preserving backwater without reverse discharge.

**This is a proposed piecewise-linear empirical rating, not an exact orifice law.** It explicitly represents resistance, discharge capacity, tailwater, and one-way flow. Its validity range requires later evidence. A more detailed orifice relationship, including partial submergence, could replace it without adding states. EPA describes nonlinear orifice ratings and user-specified outlet rating curves; this draft does not implement the manual’s orifice equations. [EPA Hydraulics Manual, §§2.1.7, 6.2](https://downloads.tuflow.com/SWMM/SWMM5_Reference_Manual_Volume2_Hydaulics_P100S9AS.pdf)

Gates:
\[
Q_{Gi}=g_i(t)\mathcal F(H_{di},H_{wi},z_{Gi};\mathcal R_{Gi},\bar Q_{Gi}),\quad i=1,2.
\]

Normal wetland drainage:
\[
Q_{wr,i}=\mathcal F(H_{wi},H_r,z_{wi,o};\mathcal R_{wi}(D_i),\bar Q_{wr,i}).
\]

Final river outlet:
\[
Q_{\rm out}=(1-B)\mathcal F(H_r,h_d,z_{r,o};\mathcal R_r,\bar Q_r).
\]

Gate opening and B affect only their own links. A closed gate leaves its spillway available; full blockage closes only the final outlet. Gate tailwater is wetland level, wetland tailwater is internal river level, and river-outlet tailwater is external h_d.

### 5.3 Dam spills and wetland overflow

\[
Q_{dsp,i}=c_{dsp,i}[H_{di}-\max(z_{dsp,i},H_{wi})]_+^{3/2},
\]
\[
Q_{wsp,i}=c_{wsp,i}[H_{wi}-\max(z_{wsp,i},H_r)]_+^{3/2}.
\]

c_sp has units m^(3/2)/s and combines effective width and discharge factors. The spill threshold z_sp is above the corresponding donor bed. These are proposed tailwater-limited power-law ratings, not exact submerged-weir formulas. Normal drainage and overflow represent distinct physical paths, so the same conveyance capacity is not counted twice.

Normal-storage thresholds:
\[
C_{di}=A_{di}(z_{dsp,i}-z_{di}),\qquad C_{wi}=A_{wi}(z_{wsp,i}-z_{wi}).
\]

These are thresholds, not hard clipping limits. Water above them remains stored until discharged. High tailwater may suppress overflow. Geometry and ratings still have finite domains of physical validity.

### 5.4 Groundwater connections

\[
Q_{wg,i}=C_{wg,i}\phi_{wi}(S_{wi})[H_{wi}-H_g]_+,
\]
\[
Q_{gr}=C_{gr}\phi_g(S_g)[H_g-H_r]_+,
\]
\[
Q_{gx}^{in}=C_{gx}[h_g-H_g]_+,\qquad
Q_{gx}^{out}=C_{gx}\phi_g(S_g)[H_g-h_g]_+,
\]
\[
Q_{gx}=Q_{gx}^{in}-Q_{gx}^{out}.
\]

The C coefficients are effective conductances in m²/s. Wetland recharge and groundwater discharge to the river retain the original one-way approximation; unfavorable receiving heads can stop them. Availability factors prevent withdrawals from empty stores merely because their reference elevations remain high.

Positive Q_gx means external groundwater enters; negative means it leaves. Outward withdrawal stops at the lower drainable-storage limit. Groundwater exfiltration into wetlands and river-bank infiltration are not represented in this draft.

## 6. Explicit temperature, radiation, and wind roles

### 6.1 Common evaporative demand

Use a simplified Penman-type wet-surface demand E_*, not standardized FAO grass-reference ET0. Fixed relative humidity r_H, radiation fraction η, and transfer coefficients are proposed uncertain parameters, not observations or new disturbance signals.

With T_a expressed in °C:
\[
T_K=T_a+273.15,
\]
\[
e_s(T_a)=610.8\exp\left(\frac{17.27T_a}{T_a+237.3}\right)\quad[\mathrm{Pa}],
\]
\[
\Delta(T_a)=\frac{17.27\times237.3\,e_s(T_a)}{(T_a+237.3)^2}\quad[\mathrm{Pa/K}],
\]
\[
D_\rho=\frac{(1-r_H)e_s(T_a)}{R_vT_K}\quad[\mathrm{kg/m^3}],
\]
\[
k_a=k_0+k_vv,\qquad R_{avail}=\eta R_s,\qquad
\gamma=\frac{c_pp_0}{\epsilon\lambda},
\]
\[
\boxed{
E_*=
\frac{\Delta}{\Delta+\gamma}\frac{R_{avail}}{\rho_w\lambda}
+\frac{\gamma}{\Delta+\gamma}\frac{k_aD_\rho}{\rho_w}
}\quad[\mathrm{m/s}].
\]

The saturation-pressure relation is documented by FAO; Δ is its derivative. The specific reduced closure remains a project assumption. Temperature affects vapor pressure, density deficit, and weighting; radiation supplies energy; wind changes air transfer. No thermal-storage state is introduced. [FAO meteorological relationships](https://www.fao.org/4/x0490e/x0490e07.htm)

FAO describes the combination of energy balance and mass transfer and explains that the standard method requires humidity. This supports the combination-method basis, not the particular simplifications adopted here. Fixed humidity and ηR_s omit changing humidity, explicit longwave exchange, and thermal storage. The result is not standardized reference ET or an established accurate hourly prediction. Condensation, snow, and freezing are outside scope. [FAO combination-method discussion](https://www.fao.org/4/x0490e/x0490e06.htm)

### 6.2 Actual withdrawals from each store

\[
E_m=A_u K_{ET,m}\phi_m(S_m)E_*,
\]
\[
E_{wi}=A_{wi}K_{ET,wi}(D_i)\phi_{wi}(S_{wi})E_*,
\]
\[
E_{di}=A_{di}K_{ET,di}\phi_{di}(S_{di})E_*,\qquad
E_r=A_rK_{ET,r}\phi_r(S_r)E_*.
\]

Every K_ET is nonnegative and dimensionless. Each loss is withdrawn once from its supplying state and vanishes at zero storage. Open-water and vegetation coefficients are distinct; coefficients calibrated to standard FAO ET0 cannot be transferred to E_* without checking compatibility.

## 7. Separate design–parameter mappings

As an explicit initial interpolation, propose:
\[
\mathcal R_{wi}(\tau_i,\rho_i)
=(1-\rho_i)\mathcal R_{bare,i}+\rho_i\mathcal R_{veg,i,\tau_i},
\]
\[
K_{ET,wi}(\tau_i,\rho_i)
=(1-\rho_i)K_{bare,i}+\rho_iK_{veg,i,\tau_i}.
\]

Each vegetation category selects its own parameter pair. Resistance endpoints are strictly positive and ET endpoints nonnegative. This is scenario interpolation of effective properties, not a universal vegetation law; endpoint ordering is not assumed.

| Wetland | Design → property → parameter → flux → ODEs |
|---|---|
| 1 | (τ1,ρ1) → hydraulic resistance → R_w1 → Q_wr,1 → Eq. 3 and Eq. 8 |
| 1 | (τ1,ρ1) → ET coefficient → K_ET,w1 → E_w1 → Eq. 3 |
| 2 | (τ2,ρ2) → hydraulic resistance → R_w2 → Q_wr,2 → Eq. 4 and Eq. 8 |
| 2 | (τ2,ρ2) → ET coefficient → K_ET,w2 → E_w2 → Eq. 4 |

Vegetation does not automatically alter seepage conductance, area, bed elevation, or overflow thresholds. Soil amendment or earthwork could later be added as explicit design components with justified mappings, without adding vegetation-growth states.

## 8. Parameters, evidence, bounds, and uncertainty

No site dimensions or observations were supplied, so defensible site-specific numerical ranges cannot be manufactured. The table gives admissible domains and value status; these domains are not empirical site ranges. Equations are evaluable once parameter values, initial conditions, and forcing are supplied.

| Symbol | Meaning | Units | Value status and admissible range | Use; uncertainty |
|---|---|---|---|---|
| A_u,A_di,A_wi,A_r | Surface and effective storage areas | m² | TBD; >0 | Rain, stage, ET; site geometry |
| z_di,z_wi,z_r,z_g | Zero-storage elevations | m | TBD; common datum | Stage; datum consistency critical |
| A_g,S_y | Aquifer area and specific yield | m²; 1 | TBD; A_g>0, 0<S_y≤1 | H_g; uncertain |
| C_m | Soil capacity | m³ | TBD; >0 | Infiltration bound; uncertain |
| S_fc | Drainage threshold | m³ | TBD; 0≤S_fc≤C_m | Recharge; uncertain |
| V_j* | Half-availability storage | m³ | TBD; >0, compartment-specific | φ_j; empirical and uncertain |
| K_u | Effective infiltration rate | m/s | TBD; ≥0 | Q_inf; highly uncertain |
| k_ud,1,k_ud,2,k_ur,k_mg | Routing and drainage rates | s⁻¹ | TBD; >0 | Runoff and recharge; uncertain |
| Q̄_ud,i,Q̄_ur | Runoff-link capacities | m³/s | TBD; >0 | Upland flow limits; geometry-dependent |
| R_Gi,Q̄_Gi,z_Gi | Gate rating | s/m²; m³/s; m | TBD; R,Q̄>0; z_Gi≥z_di | Gates; rating needs evidence |
| R_bare,i,R_veg,i,τ | Wetland resistance endpoints | s/m² | TBD; >0 | Design interpolation; uncertain |
| Q̄_wr,i,z_wi,o | Wetland normal outlet | m³/s; m | TBD; Q̄>0; z_wi,o≥z_wi | Q_wr; geometry-dependent |
| R_r,Q̄_r,z_r,o | River outlet rating | s/m²; m³/s; m | TBD; R,Q̄>0; z_r,o≥z_r | Q_out; uncertain |
| c_dsp,i,c_wsp,i | Overflow rating coefficients | m^(3/2)/s | TBD; >0 | Spills; uncertain |
| z_dsp,i,z_wsp,i | Spill thresholds | m | TBD; above donor bed; normally above normal-outlet onset | Capacities and spills |
| C_wg,i,C_gr,C_gx | Groundwater conductances | m²/s | TBD; ≥0 | Seepage and exchange; highly uncertain |
| K_ET,m,K_ET,di,K_ET,r | Other surface-loss coefficients | 1 | TBD; ≥0 | E_m,E_di,E_r; uncertain |
| K_bare,i,K_veg,i,τ | Wetland ET endpoints | 1 | TBD; ≥0 | Design → E_wi; uncertain |
| r_H | Assumed effective relative humidity | 1 | TBD; 0≤r_H≤1 | ET closure; not observed |
| η | Available-energy fraction | 1 | TBD; 0≤η≤1 | ET closure; empirical |
| k_0,k_v | Still-air and wind transfer | m/s; 1 | TBD; ≥0 | k_a; height-dependent, uncertain |
| p_0 | Fixed atmospheric pressure | Pa | TBD for elevation; >0 | γ; site-dependent |
| λ | Latent heat | J/kg | 2.45×10^6, FAO approximation near 20°C | ET physical approximation |
| c_p | Air specific heat | J/(kg·K) | 1013, FAO approximation | γ |
| ε | Molecular-mass ratio | 1 | 0.622, FAO convention | γ |
| ρ_w | Water density | kg/m³ | 1000, proposed rounded liquid-water approximation | ET conversion |
| R_v | Water-vapor gas constant | J/(kg·K) | approximately 461.5, physical constant | Vapor-density conversion |
| S_j(0) | Initial storage | m³ | TBD; ≥0 and S_m(0)≤C_m | Initial conditions; often highly uncertain |
| H_flood | River flood-reporting threshold | m | TBD from intended site/definition | Metrics; not a hard storage cap |

The stated λ, c_p, and ε conventions are documented in [FAO meteorological data](https://www.fao.org/4/x0490e/x0490e07.htm). No site measurements, synthetic storm, or unannounced design values have been assigned.

## 9. Dimensions, bounds, and conservation

### 9.1 Dimensions

| Expression | Unit result |
|---|---|
| P A | (m/s)(m²)=m³/s |
| k S | (1/s)(m³)=m³/s |
| A K_u φ | (m²)(m/s)=m³/s |
| ΔH/R | m/(s/m²)=m³/s |
| c_sp ΔH^(3/2) | (m^(3/2)/s)(m^(3/2))=m³/s |
| C_g ΔH φ | (m²/s)(m)=m³/s |
| R_avail/(ρ_w λ) | (J/s/m²)/(J/m³)=m/s |
| k_a D_ρ/ρ_w | (m/s)(kg/m³)/(kg/m³)=m/s |
| A K_ET φ E_* | (m²)(m/s)=m³/s |

Δ and γ both have units Pa/K, so the ET weights are dimensionless. Every ODE has units m³/s on both sides.

### 9.2 Nonnegativity, saturation, and capacity

- At S_u = 0, runoff and infiltration vanish; the remaining rainfall input is nonnegative.
- At S_m = 0, recharge and soil ET vanish. At S_m = C_m, infiltration vanishes while losses remain nonnegative. Therefore [0,C_m] is invariant for the continuous soil equation with admissible initial conditions.
- At S_di = 0, S_wi = 0, or S_r = 0, stage equals the corresponding bed elevation. Outlet and spill thresholds are no lower than that bed, so their withdrawals vanish. Availability factors reduce ET and seepage to zero.
- At S_g = 0, groundwater discharge to the river and the external boundary vanish; recharge and external boundary inflow remain nonnegative.
- Surface fluxes are nonnegative while tailwater still reduces drainage.
- Spill thresholds are not hard caps; storage is never clipped to them.

These are continuous-ODE checks, not guarantees for arbitrary numerical timesteps. Timestep sensitivity, physical admissibility, and residuals must be checked during implementation.

With full blockage, outlet discharge is zero and arrivals remain in S_r. No extra floodplain state or external overbank-loss path is introduced. Levels outside the river rating's valid domain, or groundwater intruding into unmodeled soil/surface storage, exceed this model's applicability.

### 9.3 Total water balance

\[
S_{tot}=S_u+S_m+S_{w1}+S_{w2}+S_{d1}+S_{d2}+S_g+S_r,
\]
\[
A_{rain}=A_u+A_{d1}+A_{d2}+A_{w1}+A_{w2}+A_r,
\]
\[
E_{tot}=E_m+E_{w1}+E_{w2}+E_{d1}+E_{d2}+E_r.
\]

\[
\boxed{\dot S_{tot}=P A_{rain}+Q_{ext}+Q_{gx}-Q_{out}-E_{tot}}.
\]

Infiltration, runoff routes, recharge, gates, dam spills, wetland drainage/overflow/seepage, and groundwater discharge to the river all cancel internally. Signed Q_gx remains a boundary term.

Future integrated residual:
\[
\mathcal E(t)=S_{tot}(t)-S_{tot}(0)
-\int_0^t\left(P A_{rain}+Q_{ext}+Q_{gx}-Q_{out}-E_{tot}\right)\,d\tau.
\]

The analytical residual is zero; numerical residuals must be measured later and have not been verified here.

## 10. Feedback, outputs, and later plots

### 10.1 Feedback present in these equations

- **Negative:** Higher pool stage increases release/spill and drains the pool; closure, caps, and tailwater can weaken this response.
- **Negative:** Wetter soil suppresses infiltration and increases drainage above its threshold.
- **Negative:** Higher wetland stage increases drainage/overflow and reduces storage.
- **Reinforcing backwater coupling:** Wetland rise can raise river inflow and stage; higher river stage suppresses wetland drainage and retains water upstream. This coupling does not imply instability or create water.
- River-to-wetland-to-gate-to-pool backwater is an upstream propagation chain, not by itself proof of unstable feedback.

### 10.2 Instantaneous outputs and event metrics

Suggested ordered output groups; final numerical-array indexing belongs to the implementation stage:

1. H_w1,H_w2,H_d1,H_d2,H_g,H_r [m].
2. Q_G1,Q_G2,Q_wr,1,Q_wr,2,Q_out [m³/s].
3. Q_dsp,1,Q_dsp,2,Q_wsp,1,Q_wsp,2 [m³/s].
4. Q_inf,Q_mg,Q_wg,1,Q_wg,2,Q_gr,Q_gx [m³/s].
5. E_tot [m³/s], S_tot [m³].

Postprocessing:
\[
Q_{peak}=\max_{t\in[0,T]}Q_{out}(t),\quad
H_{peak}=\max_{t\in[0,T]}H_r(t),
\]
\[
D_{flood}=\int_0^T\mathbf1\{H_r(t)>H_{flood}\}\,dt,\quad
V_{out}=\int_0^TQ_{out}(t)\,dt.
\]

Plot both gate histories, all eight disturbances, and states directly as well. Flood reduction must be evaluated from results, not assumed from the presence of wetland storage.

### 10.3 Algebraic plots requested by the professor

| x-axis | y-axis | Varied parameter | What it reveals |
|---|---|---|---|
| Effective gate head | Q_Gi | g_i, R_Gi | Closure, slope, capacity cap |
| H_wi | Q_wr,i and Q_wsp,i separately | Wetland design, H_r | Resistance, spill onset, tailwater |
| S_m/C_m | Q_inf at fixed S_u | K_u | Saturation shuts off infiltration |
| B | Q_out at fixed H_r,h_d | h_d | Blockage and downstream head |
| T_a | E_* | R_s,v; fixed r_H | Weather dependence |
| h_g−H_g | Q_gx at fixed S_g | C_gx,S_g | Exchange direction and availability |

## 11. Review boundary

The accepted topology is retained. Review concerns the proposed closures: fixed areas, capped linear hydraulic ratings, approximate tailwater-limited spill ratings, lumped soil/groundwater relations, availability factors, reduced ET, and vegetation interpolation.

Completed: eight explicit balance equations, explicit algebraic closures, dimensional and continuous-boundary checks, and cancellation of internal water transfers.

Not performed or claimed: calibration, defensible site-specific numerical ranges, numerical simulation, plot generation, timestep tests, hydrologic validation, optimization, or uncertainty propagation. Any later example parameters must be labeled illustrative and uncalibrated.

**Stop at equation review. After acceptance, implement using the structure of `earth_heat_2.py` and verify the actual `ode4u` interface.**
