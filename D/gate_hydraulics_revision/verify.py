"""Focused hydraulic, dimensional and full-RHS regression checks; no fitting."""
import json
from math import sqrt
from pathlib import Path
import numpy as np
from B import gate_hydraulics as original
from B import wetland_counterfactual as dry
from B import no_spillway_model as wet
from D.g5_diagnostics import runner
from .hydraulics import gate_flow, Assumptions, ApplicabilityError, contraction
from .model import make_rhs

HERE = Path(__file__).resolve().parent


def run():
    records = []
    def check(name, condition, **evidence):
        records.append(dict(check=name, passed=bool(condition), **evidence))
    fixed = runner.load_fixed()
    p, x8, gates = runner.adapter()._parameters(fixed)
    g = gates[0].gravity
    critical = sqrt(g)*(2/3)**1.5
    old_c = .6*sqrt(2*g)/1.25
    discrepancy = old_c/critical-1
    check('original endpoint discrepancy independently recomputed',
          abs(discrepancy-.247076581449591)<1e-14,
          old_coefficient_SI=old_c, ideal_coefficient_SI=critical,
          excess_percent=100*discrepancy, energy_head='pool WSE minus sill')
    for cu in (2.6,2.8,3.0):
        settings=Assumptions(Cw_US=cu)
        ideal,cw,loss,contact=settings.coefficients(g)
        check(f'{cu}: SI conversion and nonnegative loss',loss>=0,
              Cw_SI=cw, loss=loss, critical_depth_head_ratio=2/(3+loss))
        check(f'{cu}: critical-flow derivation',abs(sqrt(g)*(2/(3+loss))**1.5-cw)<1e-14)
        for gate in gates:
            z=gate.gate_sill_z
            for label,up,tail,opening in [('closed',z+.1,z-1,0),
                    ('dry',z,z-1,.5),('below_sill',z-.1,z-1,.5),
                    ('equal',z+.1,z+.1,.5),('adverse',z+.1,z+.2,.5)]:
                check(f'{cu}/{gate.gate_width}: {label}',gate_flow(up,tail,opening,gate,settings).discharge==0)
            # Fixed H sweep tests all requested old and newly revised boundaries.
            H=.12
            rs=np.unique(np.r_[np.linspace(.5,20,3000),1.,1.25,1.5])
            fs=[gate_flow(z+H,z-1,H/r/gate.gate_max_opening,gate,settings) for r in rs]
            qs=np.array([f.discharge for f in fs])
            check(f'{cu}/{gate.gate_width}: energy bound/nonnegative',
                  np.all(qs>=0) and max(qs)<=ideal*gate.gate_width*H**1.5*(1+1e-10))
            check(f'{cu}/{gate.gate_width}: opening monotonicity',np.max(np.diff(qs))<=1e-12)
            for r in (1.,1.25,1.5):
                eps=1e-8
                values=[gate_flow(z+H,z-1,H/(r+e)/gate.gate_max_opening,gate,settings).discharge for e in (-eps,0,eps)]
                check(f'{cu}/{gate.gate_width}: continuity at {r}',max(values)-min(values)<1e-7)
            # Fixed opening, increase H across every active regime. No smoothing.
            a=.1;heads=np.linspace(.05,.5,5000)
            q=np.array([gate_flow(z+h,z-1,a/gate.gate_max_opening,gate,settings).discharge for h in heads[1:]])
            check(f'{cu}/{gate.gate_width}: head monotonicity',np.min(np.diff(q))>=-1e-12,
                  minimum_discharge_increment_m3s=float(np.min(np.diff(q))))
            # Froude similarity Q scales as length^(5/2) when b,a,H scale.
            f1=gate_flow(z+.10,z-1,.2,gate,settings).discharge
            gate2=original.GateParameters(gravity=g,gate_width=2*gate.gate_width,
                gate_max_opening=2*gate.gate_max_opening,gate_sill_z=z,
                Cd_gate_free=.6,Cd_gate_submerged=.8,neglect_approach_velocity_head=True)
            f2=gate_flow(z+.20,z-1,.2,gate2,settings).discharge
            check(f'{cu}/{gate.gate_width}: dimensional similarity',abs(f2/f1-2**2.5)<1e-10)
            for label,h,tail,opening in [('shallow',.049,z-1,.8),('drowned',.2,z+.01,.8)]:
                try:
                    gate_flow(z+h,tail,opening,gate,settings)
                    rejected=False
                except ApplicabilityError:
                    rejected=True
                check(f'{cu}/{gate.gate_width}: explicit {label} rejection',rejected)
    for case in ('wetland','no_wetland'):
        x=x8.copy() if case=='wetland' else x8[dry.SHARED].copy()
        Y=runner.adapter().a2.Y if case=='wetland' else dry.Y
        base=wet.make_b2_rhs(gates) if case=='wetland' else dry.make_rhs(gates)
        u=np.array([.005,.005,1e-7,20,100,2,0,97.15,97.15,0])
        for cu in (2.6,2.8,3.0):
            d,y=make_rhs(case,gates,{'Cw_US':cu})(0,x,u,p)
            db,yb=base(0,x,u,p)
            unaffected=[j for j in range(len(y)) if j not in (Y['Q_G1'],Y['Q_G2'])]
            check(f'{case}/{cu}: every non-gate output unchanged',np.array_equal(y[unaffected],yb[unaffected]))
            expected=db.copy()
            for i,di in enumerate((4,5) if case=='wetland' else (2,3)):
                delta=y[Y[f'Q_G{i+1}']]-yb[Y[f'Q_G{i+1}']]
                expected[di]-=delta
                expected[2+i if case=='wetland' else 5]+=delta
            check(f'{case}/{cu}: only internal gate transfer changes',np.max(np.abs(expected-d))<1e-14)
            check(f'{case}/{cu}: full ODE physical mass identity',abs(sum(d)-y[Y['external_net']])<1e-13,
                  residual_m3_s=float(sum(d)-y[Y['external_net']]))
        # Closed-gate RHS exactly reuses original behavior.
        u[:2]=0
        d,y=make_rhs(case,gates)(0,x,u,p);db,yb=base(0,x,u,p)
        check(f'{case}: closed-gate regression',np.array_equal(d,db) and np.array_equal(y,yb))
        # Exercise the complete balances in the newly entered flow regimes,
        # rather than checking conservation only at the high-head initial state.
        for ratio in (.75,1.,1.1,1.25,1.5,3.):
            state=x.copy();uin=u.copy();uin[:2]=.12/ratio/.3
            for i,di in enumerate((4,5) if case=='wetland' else (2,3)):
                state[di]=p.A_d[i]*(gates[i].gate_sill_z+.12-p.z_d[i])
            d,y=make_rhs(case,gates)(0,state,uin,p)
            closed=uin.copy();closed[:2]=0
            dc,yc=base(0,state,closed,p)
            expected=dc.copy()
            for i,di in enumerate((4,5) if case=='wetland' else (2,3)):
                q=y[Y[f'Q_G{i+1}']]
                expected[di]-=q
                expected[2+i if case=='wetland' else 5]+=q
            check(f'{case}: low-head full balance r={ratio}',
                  abs(sum(d)-y[Y['external_net']])<1e-13 and np.max(np.abs(d-expected))<1e-14)
    report=dict(passed=all(r['passed'] for r in records),checks=records,
                check_count=len(records),scope='Numerical/physical consistency checks, not empirical validation')
    (HERE/'results').mkdir(exist_ok=True)
    (HERE/'results'/'hydraulic_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(passed=report['passed'],checks=len(records),failed=[r['check'] for r in records if not r['passed']])))
    return report


if __name__=='__main__':
    raise SystemExit(0 if run()['passed'] else 1)
