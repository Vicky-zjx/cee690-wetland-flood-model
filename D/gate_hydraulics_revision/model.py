"""D-only gate replacement; all other accepted B constitutive laws are reused."""
import numpy as np
from A import wetland_flood_model as a2
from B import no_spillway_model as wet
from B import wetland_counterfactual as dry
from .hydraulics import Assumptions, ApplicabilityError, gate_flow


def make_rhs(case, gates, variant=None, diagnostics=None):
    if case not in ('wetland', 'no_wetland'):
        raise ValueError('Unknown physical model')
    settings = Assumptions(**(variant or {}))
    base = wet.make_b2_rhs(gates) if case == 'wetland' else dry.make_rhs(gates)
    Y = a2.Y if case == 'wetland' else dry.Y
    d_indices = (4, 5) if case == 'wetland' else (2, 3)
    def rhs(t, x, u, p):
        closed = np.array(u, copy=True)
        closed[:2] = 0.0
        dx, y = base(t, x, closed, p)
        for i in range(2):
            head = float(y[Y[f'H_d{i+1}']])
            tail = float(y[Y[f'H_w{i+1}'] if case == 'wetland' else Y['H_r']])
            try:
                flow = gate_flow(head, tail, u[i], gates[i], settings)
            except ApplicabilityError as exc:
                record = dict(exc.details, branch=i+1, time_s=float(t),
                    normalized_opening=float(u[i]), regime=exc.details['reason'],
                    residual_m=1.25*gates[i].gate_max_opening*float(u[i])-(head-gates[i].gate_sill_z),
                    rule='Stage D physical domain; legacy G5 is reported only')
                if diagnostics is not None:
                    diagnostics.setdefault('first_applicability_failure', record)
                raise wet.HydraulicApplicabilityFailure(record) from exc
            q = flow.discharge
            dx[d_indices[i]] -= q
            dx[2+i if case == 'wetland' else 5] += q
            y[Y[f'Q_G{i+1}']] = q
            if diagnostics is not None:
                counts = diagnostics.setdefault('regime_rhs_counts', {})
                key = f'gate_{i+1}:{flow.regime}'
                counts[key] = counts.get(key, 0)+1
                ranges = diagnostics.setdefault('energy_head_range_m', {})
                old = ranges.setdefault(str(i+1), [flow.energy_head, flow.energy_head])
                old[0] = min(old[0], flow.energy_head)
                old[1] = max(old[1], flow.energy_head)
        return dx, y
    return rhs
