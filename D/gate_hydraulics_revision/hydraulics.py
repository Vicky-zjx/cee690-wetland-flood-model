"""Bounded experimental Stage D closure; never a change to the accepted B law.

See geometry_and_method.md for sources, derivation and the deliberately exposed
empirical transition assumption. Unsupported crest/tailwater states raise; no
opening, storage, discharge or safety-limit clipping is performed.
"""
from dataclasses import dataclass
from math import isfinite, sqrt
from B.gate_hydraulics import GateFlow, GateInputError


class ApplicabilityError(RuntimeError):
    def __init__(self, reason, **details):
        self.details = dict(reason=reason, **details)
        super().__init__(str(self.details))


@dataclass(frozen=True)
class Assumptions:
    crest_length_m: float = 1.0
    Cw_US: float = 2.8
    contraction_at_zero_relative_opening: float = 0.617
    transition_end_ratio: float = 1.5
    minimum_crest_head_length_ratio: float = 0.05
    maximum_crest_head_length_ratio: float = 0.5

    def __post_init__(self):
        if not all(isfinite(v) for v in self.__dict__.values()):
            raise ValueError('All hydraulic assumptions must be finite')
        if self.Cw_US not in (2.6, 2.8, 3.0):
            raise ValueError('Use the declared, pre-run coefficient sensitivity set')
        if (self.crest_length_m != 1.0 or self.contraction_at_zero_relative_opening != 0.617 or
                self.transition_end_ratio != 1.5 or
                self.minimum_crest_head_length_ratio != 0.05 or
                self.maximum_crest_head_length_ratio != 0.5):
            raise ValueError('Geometry and domain may not be tuned by the run configuration')

    def coefficients(self, gravity):
        ideal = sqrt(gravity) * (2.0/3.0)**1.5
        cw = self.Cw_US * sqrt(0.3048)
        loss = 3.0*((ideal/cw)**(2.0/3.0)-1.0)
        if loss < 0:
            raise ValueError('Negative loss is outside this energy model')
        return ideal, cw, loss, (3.0+loss)/2.0


def contraction(relative_opening):
    """USBR PAP-1094 Eq. 2, vertical gate; only 0 <= a/H <= 2/3."""
    s = float(relative_opening)
    if not 0 <= s <= 2/3+1e-14:
        raise ValueError('Contraction relation outside declared free-orifice range')
    return 0.617-0.06*s+0.026*s*s+0.026*s*s*s


def gate_flow(upstream, downstream, opening, gate, assumptions=None):
    assumptions = assumptions or Assumptions()
    upstream, downstream, opening = float(upstream), float(downstream), float(opening)
    if not all(isfinite(v) for v in (upstream, downstream, opening)) or not 0 <= opening <= 1:
        raise GateInputError('Finite elevations and an unchanged opening in [0,1] are required')
    a = gate.gate_max_opening*opening
    H = upstream-gate.gate_sill_z
    dH = upstream-downstream
    sub = (downstream-gate.gate_sill_z)/H if H > 0 else None
    def out(q, regime):
        return GateFlow(q, regime, a, H, dH, sub)
    if a == 0: return out(0.0, 'closed')
    if H <= 0: return out(0.0, 'dry')
    if dH <= 0: return out(0.0, 'adverse_head_one_way')
    details = dict(upstream_stage_m=upstream, tailwater_stage_m=downstream,
                   water_depth_m=H, opening_height_m=a, head_opening_ratio=H/a,
                   legacy_G5_residual_m=1.25*a-H)
    # Safety-feasible eight-state tails and applicable six-state river levels
    # lie below the sill. Do not invent a drowned-crest or reverse-flow law.
    if downstream >= gate.gate_sill_z:
        raise ApplicabilityError('tailwater_at_or_above_sill', **details)
    ideal, cw, loss, contact = assumptions.coefficients(gate.gravity)
    r = H/a
    if r < assumptions.transition_end_ratio:
        h_over_l = H/assumptions.crest_length_m
        if not assumptions.minimum_crest_head_length_ratio <= h_over_l <= assumptions.maximum_crest_head_length_ratio:
            raise ApplicabilityError('broad_crest_head_length_outside_supported_range',
                head_length_ratio=h_over_l, required_range=[0.05, 0.5], **details)
    b = gate.gate_width
    qw = cw*b*H**1.5
    if r <= 1.0:
        return out(qw, 'free_weir')
    if r >= assumptions.transition_end_ratio:
        jet = contraction(1/r)*a
        qg = b*jet*sqrt(2*gate.gravity*(H-jet)/(1+loss))
        return out(qg, 'energy_gate')
    # Similarity-based endpoint interpolation. The gate equation is evaluated
    # ONLY at its supported endpoint r=1.5, never extrapolated to a/H>2/3.
    # This is a declared engineering closure inspired by HEC-RAS interpolation;
    # it is NOT HEC-RAS's exact implementation or a measured gate rating.
    s_end = 1/assumptions.transition_end_ratio
    d_end = contraction(s_end)
    c_end = d_end*s_end*sqrt(2*gate.gravity*(1-d_end*s_end)/(1+loss))
    weight = (r-1)/(assumptions.transition_end_ratio-1)
    q = b*H**1.5*((1-weight)*cw+weight*c_end)
    # Diagnostic assertions expose an incompatible closure; they do not cap it.
    # In this transition a >= 2H/3, so the complete critical-flow depth fits
    # below the opening; do not incorrectly impose a full-opening flow depth.
    capacity = ideal*b*H**1.5
    if q > capacity*(1+1e-12):
        raise ApplicabilityError('transition_energy_inconsistent', discharge_m3_s=q,
                                 opening_capacity_m3_s=capacity, **details)
    return out(q, 'assumed_weir_gate_transition')
