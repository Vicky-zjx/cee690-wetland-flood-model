# Stage D upstream-level feedback policy: basis and assumptions

## Purpose and preserved model

This is a separate, explicitly hypothetical operating-policy experiment. It replaces the time-only gate commands with local upstream-water-level feedback, while reusing the existing Stage D hydraulic candidate unchanged. A/B/C, previous D studies, physical initial storages, meteorology, boundaries, wetland designs, geometry, coefficients and safety limits remain unchanged. No policy optimization or scenario-specific tuning is part of this case.

The controller selects the commanded opening. The hydraulic function receives the actual opening and computes discharge with its existing applicability checks. No discharge cap, state clipping, head-to-opening clipping or relaxed hydraulic domain is introduced. Successful operation within a supported domain would demonstrate the tested controller/model combination; it would not validate the hydraulic transition or eliminate every possible G5-related failure.

## Primary operational references

The [USACE HEC-RAS 7.0 elevation-controlled-gate guidance](https://www.hec.usace.army.mil/confluence/cwmsdocs/rasum/7.0/boundary-conditions-150733912.html), under “Elevation Controlled Gate,” supports upstream-level operation with separately specified opening and closing elevations, movement rates, maximum/minimum openings and an initial opening. Using separate thresholds with retained operating mode provides the hysteresis used here. The guidance does not supply this project's numerical thresholds.

The [USACE rule-operation setup](https://www.hec.usace.army.mil/confluence/rasdocs/rasum/7.0/advanced-features-for-unsteady-flow-routing/user-defined-rules-for-hydraulic-structures-and-pumps/entering-rule-operations) distinguishes a target setting from finite actuator travel, with explicit movement rates and position limits. Its illustrative 0.1 ft/min rate equals 0.000508 m/s. The present 0.0005 m/s assumption is a rounded analogue of that example, not a manufacturer specification or a recommended speed for a real installation. USACE also evaluates rule decisions once per computational step; this implementation instead uses a fixed 120 s controller clock, held constant during ODE refinement.

The [USACE operation-rules reference](https://www.hec.usace.army.mil/confluence/rasdocs/rasum/7.0/advanced-features-for-unsteady-flow-routing/user-defined-rules-for-hydraulic-structures-and-pumps/operation-rules) supports retained rule variables, current model-state inputs and target gate positions. The present implementation changes opening commands only; it does not use the separate capability to impose a desired discharge independently of hydraulics.

## Declared policy, identical for both branches

Let H_i be the current dam-pool level above gate sill i, a_i the actual opening in metres, and m_i the retained ON/OFF operating mode. Both branches have the same allowable pool depth: 1.20 m above their respective sills.

| Setting | Value | Classification and reason |
|---|---:|---|
| Controller decision interval | 120 s | Numerical/control assumption, independent of ODE step |
| Close threshold | H_i <= 0.20 m | Hypothetical operating trigger, above the 0.05 m shallow-crest applicability boundary |
| Reopen threshold | H_i >= 0.30 m | Hypothetical operating trigger; 0.10 m hysteresis band |
| Mode inside 0.20 < H_i < 0.30 m | Retain previous mode | Prevents immediate open/close reversal at a single threshold |
| Initial operating modes | ON / ON | Explicit initialization consistent with initial H_i = 0.30 m |
| Initial actual opening | 0.0015 m each | Preserved initial normalized opening 0.005 |
| OFF target | 0 m | Closed-gate command; not instantaneous closure |
| Lowest ON target | 0.10 m | Engineering command level; not a lower bound on actual opening |
| Full-opening target reached at | H_i >= 0.80 m | Engineering trigger, leaving 0.40 m to the unchanged dam safety depth |
| Maximum actual opening | 0.30 m | Preserved physical maximum |
| Maximum opening/closing speed | 0.0005 m/s | Hypothetical actuator assumption; approximately the USACE illustrative rate |

The absolute close/reopen elevations are 100.20/100.30 m for Gate 1 and 100.00/100.10 m for Gate 2. Full opening is commanded at 100.80 and 100.60 m. The hard dam safety limits remain 101.20 and 101.00 m; wetland safety limits also remain unchanged. Operating triggers are not replacement safety limits.

At each controller time, update the mode using only the current accepted state:

    if H_i <= 0.20: m_i = OFF
    else if H_i >= 0.30: m_i = ON
    else: retain m_i

Then set the target opening:

    OFF:                 a_target = 0
    ON and H_i <= 0.30:   a_target = 0.10
    ON and 0.30 < H_i < 0.80:
                         a_target = 0.10 + 0.20 (H_i - 0.30)/0.50
    ON and H_i >= 0.80:   a_target = 0.30

The rising target is a transparent engineering choice to release more water as the pool rises. Its breakpoints and slopes are not fitted to a weather scenario, river-peak improvement or completion outcome. There is no separate emergency override: the target is already fully open above 0.80 m. Whether available discharge and actuator response keep the pool safe requires actual simulation.

## Actuator trajectory and decision timing

For control interval Delta = 120 s and speed limit v = 0.0005 m/s, calculate

    a_next = a_current + sign(a_target - a_current)
             min(abs(a_target - a_current), v Delta).

Actual opening is linear from a_current to a_next over the entire interval. If the remaining travel is less than v Delta, the controller deliberately uses a slower constant speed and reaches the endpoint at the next decision time. It does not move at maximum speed and then stop early. This distinction defines the reproducible input trajectory.

The maximum change per decision interval is 0.06 m; full 0.30 m travel takes at least 600 s. Position bounds arise from the actuator target and movement rule. They do not modify the hydraulic equation or depend on a computed discharge cap. The actual opening can be below 0.10 m during startup or closure; the ON target floor must not be reported as a physical minimum opening.

Commands use accepted upstream levels at controller times, not future weather, subsequent states or intermediate RK trial evaluations. ODE refinements retain the same controller clock, hysteresis and rate limits. Perfect upstream-level measurement, no communication delay beyond the decision clock, reliable power, and a speed-controlled actuator with no acceleration/inertia model are additional idealizations. Mechanical force, torque, wear and sensor accuracy have not been established.

## Pre-run hydraulic rationale and its limits

For an ON target, the piecewise rule gives H_i/a_target >= 2 before closure. In the ramp region it decreases from 3 to 8/3; above 0.80 m it is at least 8/3. This concerns the commanded target only. Actual actuator lag can produce a different ratio, so the hydraulic checks remain active at every evaluation. The controller never sets a_target = min(command, H_i/constant) inside the discharge calculation.

The close threshold leaves 0.15 m above the existing shallow-crest bound. A conservative hydraulic-release-only screening calculation uses the ideal critical-flow capacity at H = 0.20 m, a full 120 s detection delay and 600 s closure:

    Delta H_release <= Ccrit b (0.20)^(3/2) (120 + 600) / A_d
                    = 0.010978 m,
    Ccrit = 1.704603788 m^(1/2)/s.

Both branches give the same value because b/A_d = 1/10,000 = 0.8/8,000 = 0.0001 m^-1. The screen assumes head is at most 0.20 m while draining and uses the full ideal-capacity bound throughout closure. It omits evaporation and is not a full invariant proof or a flood-safety guarantee. It explains the selected buffer independently of scenario outcomes; actual low-head, rate and safety behavior must still be checked.

The completed [policy checks](results/policy_checks.json) also retain an evaporation-inclusive version of this screen. The unchanged ET equation is bounded componentwise over the saved deterministic forcing and the 16 non-freezing weather scenarios, using maximum radiation, wind and vapor-pressure terms and minimum absolute temperature. The resulting demand envelope is 5.878783827e-7 m/s (above the largest directly sampled demand, 1.876186450e-7 m/s). Including the existing dam evaporation multipliers over 720 s gives a release-plus-ET drawdown bound of **0.011400709 m** in either branch and leaves **0.138599291 m** above the 0.05 m boundary. This is a deterministic screen conditional on those saved input ranges and the existing ET equation, not a guarantee for arbitrary future weather or physical actuator behavior; it does not alter the policy or hydraulic law.

The preserved initial 1.5 mm opening and brief small-opening passages during closure remain outside any demonstrated local actuator/hydraulic calibration. Closed gates may allow pool levels to decline through evaporation; this is distinct from evaluating an unsupported positive gate flow. The existing dry/closed behavior and strict hydraulic exceptions remain unchanged.

## Comparison and reporting requirements

Use this same policy, controller initialization and actuator assumptions in both wetland and no-wetland models. Save actual opening histories and mode changes. A common feedback policy does not automatically imply identical realized gate schedules; compare the histories before describing a downstream difference as a wetland effect under identical operation.

Report hydraulic applicability, completion and safety separately. A complete simulation that exceeds a safety limit remains an unsafe result, not a hydraulic failure. A failure to complete within the unchanged hydraulic domain is retained with its first diagnostic. Any success applies only to the tested weather, geometry and hypothetical policy; it is neither proof of safety under all disturbances nor validation of the existing Stage D hydraulic rating.
