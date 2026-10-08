"""Documented B.2 sluice/orifice branches; no project geometry defaults.

Sources and deliberately unsupported regimes are in docs/07_gate_hydraulics.md.
This module uses only the standard library and performs no simulation or I/O.
It does not modify the shared A.2 head_flow() rating or its separate spillways.
"""
from dataclasses import dataclass
from math import isfinite, sqrt


class GateConfigurationError(ValueError):
    """Missing/invalid fixed gate evidence or an unsupported head assumption."""


class GateInputError(ValueError):
    """A nonfinite level or normalized opening outside [0, 1]."""


class UnsupportedGateRegime(RuntimeError):
    """Hydraulic state needs an unimplemented law, not a clipped discharge.

    Rejecting this state is a model-coverage failure. It does not authorize
    shrinking the user's optimization domain to avoid that state.
    """
    def __init__(self, regime, head_to_opening_ratio, submergence):
        self.regime = regime
        self.head_to_opening_ratio = head_to_opening_ratio
        self.submergence = submergence
        super().__init__(
            f"Unsupported gate regime: {regime}; H_energy/a="
            f"{head_to_opening_ratio:.12g}, tailwater_depth/H_energy="
            f"{submergence:.12g}. Supply a documented applicable closure or "
            "obtain explicit approval for any operating-domain restriction; "
            "do not substitute a minimum, switch, cap, or clipped opening."
        )


def _finite_number(value, name, error_type):
    if value is None or isinstance(value, bool):
        raise error_type(f"{name} requires an explicit finite numerical value.")
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise error_type(f"{name} requires an explicit finite numerical value.") from exc
    if not isfinite(number):
        raise error_type(f"{name} must be finite.")
    return number


@dataclass(frozen=True)
class GateParameters:
    """One branch's fixed configuration, in SI; None is intentionally missing.

    Coefficients must match the head and geometry conventions in the sources.
    No coefficient, sill, width, opening, or gravity value is inferred from A.2.
    The caller must verify gate_sill_z >= the donor pool's zero-storage bed.
    """
    gravity: float | None = None                 # m/s2
    gate_width: float | None = None              # m
    gate_max_opening: float | None = None        # m
    gate_sill_z: float | None = None             # m, common vertical datum
    Cd_gate_free: float | None = None            # dimensionless
    Cd_gate_submerged: float | None = None       # dimensionless, independently set
    neglect_approach_velocity_head: bool | None = None

    def __post_init__(self):
        for name in ("gravity", "gate_width", "gate_max_opening",
                     "gate_sill_z", "Cd_gate_free", "Cd_gate_submerged"):
            value = _finite_number(getattr(self, name), name, GateConfigurationError)
            if name != "gate_sill_z" and value <= 0:
                raise GateConfigurationError(f"{name} must be strictly positive.")
            object.__setattr__(self, name, value)
        if self.neglect_approach_velocity_head is not True:
            raise GateConfigurationError(
                "Set neglect_approach_velocity_head=True only when explicitly "
                "accepting pool water-surface elevation as upstream energy "
                "elevation. Non-negligible approach velocity is not implemented."
            )


@dataclass(frozen=True)
class GateFlow:
    discharge: float                  # m3/s, positive pool -> associated wetland
    regime: str
    opening_height: float             # m, exactly gate_max_opening * normalized opening
    energy_head: float                # m above sill, NOT absolute elevation
    differential_head: float          # m: upstream energy elevation - downstream WSE
    submergence: float | None         # downstream depth above sill / energy head


def sluice_gate_flow(upstream_level, downstream_level, normalized_opening,
                     params: GateParameters) -> GateFlow:
    """Evaluate one gate once, returning its discharge and hydraulic regime.

    Implemented positive-flow domain: H_energy/a >= 1.25. Submergence
    <= 0.67 uses free sluice flow, >= 0.80 uses submerged orifice flow,
    and the intervening interval uses the explicit USACE 7.0 weighted law.
    Low-depth/weir and weir-gate transition states remain unsupported.

    HEC-RAS 6.6 sluice Eq. 1: Q = Cd_free * width * a * sqrt(2*gravity*H).
    Linked radial-gates Eq. 3: Q = Cd_sub * width * a * sqrt(2*gravity*dH).
    The intervening transition follows the explicit USACE 7.0 overview.
    See documentation for head conventions and the 6.6 overview typo.

    The one-way project assumption returns zero at equal/adverse head. That
    routing assumption is distinct from either hydraulic discharge equation.
    """
    if not isinstance(params, GateParameters):
        raise GateConfigurationError("params must be an explicit GateParameters instance.")
    upstream = _finite_number(upstream_level, "upstream_level", GateInputError)
    downstream = _finite_number(downstream_level, "downstream_level", GateInputError)
    opening = _finite_number(normalized_opening, "normalized_opening", GateInputError)
    if not 0 <= opening <= 1:
        raise GateInputError("normalized_opening must lie in [0, 1]; it is not clipped.")

    # Explicit modeling assumption: linear opening map, negligible approach head.
    height = params.gate_max_opening * opening
    energy_head = upstream - params.gate_sill_z
    differential_head = upstream - downstream
    submergence = ((downstream - params.gate_sill_z) / energy_head
                   if energy_head > 0 else None)

    def result(discharge, regime):
        if not isfinite(discharge):
            raise GateInputError("Nonfinite gate discharge; inspect the supplied scales.")
        return GateFlow(discharge, regime, height, energy_head,
                        differential_head, submergence)

    if opening == 0:
        return result(0.0, "closed")
    if energy_head <= 0:
        return result(0.0, "dry")
    if differential_head <= 0:
        return result(0.0, "adverse_head")

    head_to_opening_ratio = energy_head / height
    if head_to_opening_ratio <= 1:
        raise UnsupportedGateRegime("weir_or_partially_wetted_opening",
                                    head_to_opening_ratio, submergence)
    if head_to_opening_ratio < 1.25:
        raise UnsupportedGateRegime("weir_gate_transition",
                                    head_to_opening_ratio, submergence)

    area = params.gate_width * height
    if submergence <= 0.67:
        discharge = params.Cd_gate_free * area * sqrt(2 * params.gravity * energy_head)
        return result(discharge, "free")
    if submergence < 0.80:
        # USACE 7.0 overview explicitly specifies both weights and heads.
        # Its rounded 0.67 boundary produces a small documented jump; do not
        # silently replace 0.67 with 2/3 or smooth the hydraulic relationship.
        fraction = (submergence - 0.67) / (0.80 - 0.67)
        free_component = (params.Cd_gate_free * area
                          * sqrt(2 * params.gravity * 3 * differential_head))
        submerged_component = (params.Cd_gate_submerged * area
                               * sqrt(2 * params.gravity * differential_head))
        discharge = ((1 - fraction) * free_component
                     + fraction * submerged_component)
        return result(discharge, "transition")
    discharge = (params.Cd_gate_submerged * area
                 * sqrt(2 * params.gravity * differential_head))
    return result(discharge, "submerged")
