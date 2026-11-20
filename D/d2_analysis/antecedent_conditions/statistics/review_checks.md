# Independent method and saved-array review

This review reads the new antecedent adapter/provider/analysis and saved outputs. It does not execute a physical model, change parameters, or alter any original file. The importance-weight audit in this directory concerns the earlier 64-vector ensemble and is separate from the designed antecedent cases.

## Method checks

- The adapter reuses the unchanged physical right-hand side, official solver, feedback policy, actuator rates and 120 s controller clock. Its additional wrappers observe ordinary solver evaluations; they do not evaluate a second physical trajectory.
- Accepted initial states are imposed at the beginning of the 24 h or 72 h prefix. The event begins inside one continuous integration. State, actual gate opening and pre-update latch memory must match across that boundary exactly; the legitimate controller update at the boundary is not mistaken for a reset.
- The provider maps the prefix and original event separately, shifting the event clock without changing its saved weather or boundary values. The resulting temperature/wind seam jumps are explicitly documented. Rainfall and solar energy are audited by completed control blocks.
- Phase-relative mass error uses the phase start as its storage reference and its own external-inflow scale. A phase error is not allowed to hide behind a larger full-run denominator. Dam and wetland safety checks retain both accepted-node and successful-RHS evidence.
- Full extended-horizon and original 144 h event-window metrics are separate. The event-start sample remains eligible for every maximum. Different evolved event-start conditions in the wetland/no-wetland systems are legitimate consequences of continuous conditioning, not a reason to discard a pair.
- Pairwise differences require matching clocks, weather and shared physical/policy provenance. Actual gate-history equality is a reported outcome, not an extra pairing criterion.

The accounting sign is correct. With the reported convention Δ = no wetland minus wetland,

`ΔV_out = ΔS_start − ΔS_end + ΔN_nonoutlet + Δmass_residual`.

The near-zero error of this reconstructed identity is an arithmetic consistency check. Independent evidence for numerical conservation comes from direct block-integrated fluxes, phase mass residuals, and successful-RHS physical mass identities. It must not be described as a second independent physical validation.

## Saved first-pilot check

For August 2024, a 72 h prefix and saved vector `mc_00018`, both completed 15 s model outputs were independently read from NPZ/JSON/CSV without importing or running the model:

| Check | Eight-state wetland | Six-state no wetland |
|---|---:|---:|
| Original event-relative timestamps bit-identical | Yes | Yes |
| Original event weather and boundary arrays bit-identical | Yes | Yes |
| Prefix-start physical states equal accepted original states | Yes | Yes |
| Event external-net volume: reconstructed minus direct blocks (m³) | −3.09×10⁻¹¹ | −1.46×10⁻¹¹ |
| Event-start river stage (m) | 97.15960716473992 | 97.15811104944115 |
| Event peak including event start (m) | 97.34204450968306 | 97.47059563232564 |
| Minimum successful open-gate H/a during event | 1.9943988751 | 1.9943988751 |

Both 15 s trajectories are classified as completed, valid and safety-feasible. Their state/actuator/latch handoff flags are exact. These checks cover the first saved coarse pilot only; they do not establish a completed 12-pair study or a passed 15/7.5 s comparison.

## Review corrections and reporting cautions

The new runner's obsolete inherited pilot-selection text was removed; its declared numerical pilot is the specified single long-prefix pair. The analysis now checks pairing provenance before reporting differences and labels full-horizon gate diagnostics explicitly. Direct phase flux and residual fields are retained in the analysis.

One minor summary-count issue was reported: full-horizon latch mismatch counts should be aggregated once per physical pair, not summed again across its event/antecedent window rows. The first pilot's zero mismatches are unaffected. No further physical-method blocker was found in this source review. Final saved-result accounting remains to be checked after the intended runs finish.
