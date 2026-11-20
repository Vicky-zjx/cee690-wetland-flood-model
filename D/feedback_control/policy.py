"""Explicit sampled controller and actuator plan, separate from gate hydraulics.

Two memories are retained outside the water-storage ODE: enabled/off latch and
actual actuator opening. Commands use only the accepted upstream level at the
fixed decision clock. No RHS trial state can mutate control memory.
"""
from pathlib import Path
import json
from math import isfinite

CONTROL_PERIOD_S = 120.0
POLICY_PATH = Path(__file__).with_name('policy.json')


class Controller:
    def __init__(self, actual_opening_m, config=None):
        self.config = json.loads(POLICY_PATH.read_text()) if config is None else dict(config)
        c = self.config
        if c['control_period_s'] != CONTROL_PERIOD_S:
            raise ValueError('Keep the controller clock fixed independently of the ODE grid')
        if not (0 < c['close_head_m'] < c['reopen_head_m'] < c['full_open_head_m'] < min(c['unchanged_dam_safe_head_m'])):
            raise ValueError('Invalid or non-conservative ordering of controller thresholds')
        if not (0 < c['enabled_low_head_target_m'] < c['maximum_opening_m'] == .3):
            raise ValueError('Retain existing maximum opening and valid enabled target')
        if min(c['opening_rate_limit_m_per_s'],c['closing_rate_limit_m_per_s']) <= 0:
            raise ValueError('Actuator speeds must be finite and positive')
        self.actual = [float(v) for v in actual_opening_m]
        if len(self.actual) != 2 or any(not isfinite(v) or not 0 <= v <= c['maximum_opening_m'] for v in self.actual):
            raise ValueError('Supply two actual initial openings; values are not clipped')
        self.enabled = list(c['initial_enabled'])
        if len(self.enabled) != 2 or any(type(v) is not bool for v in self.enabled):
            raise ValueError('Supply the two initial latch states explicitly')
        self.next_time = 0.0

    @classmethod
    def from_initial(cls, actual_opening_m, config=None):
        return cls(actual_opening_m, config)

    def update(self, heads_m, time_s):
        c = self.config
        if not isfinite(time_s) or abs(time_s-self.next_time)>1e-9:
            raise ValueError('Controller updates must occur once per fixed 120 s interval')
        heads = [float(v) for v in heads_m]
        if len(heads)!=2 or not all(isfinite(v) for v in heads):
            raise ValueError('Two finite upstream energy heads are required')
        start = self.actual.copy(); before = self.enabled.copy()
        targets=[]; ends=[]; rates=[]; modes=[]
        for i,H in enumerate(heads):
            if H <= c['close_head_m']:
                self.enabled[i] = False
            elif H >= c['reopen_head_m']:
                self.enabled[i] = True
            if not self.enabled[i]:
                target=0.0; mode='closing_or_closed_low_water'
            elif H <= c['reopen_head_m']:
                target=c['enabled_low_head_target_m']; mode='enabled_low_head_release'
            elif H >= c['full_open_head_m']:
                target=c['maximum_opening_m']; mode='full_open_requested'
            else:
                fraction=(H-c['reopen_head_m'])/(c['full_open_head_m']-c['reopen_head_m'])
                target=c['enabled_low_head_target_m']+fraction*(c['maximum_opening_m']-c['enabled_low_head_target_m'])
                mode='level_proportional_release'
            delta=target-start[i]
            speed=c['opening_rate_limit_m_per_s'] if delta>=0 else c['closing_rate_limit_m_per_s']
            travel=speed*CONTROL_PERIOD_S
            end=target if abs(delta)<=travel else start[i]+(travel if delta>0 else -travel)
            targets.append(target);ends.append(end);rates.append((end-start[i])/CONTROL_PERIOD_S);modes.append(mode)
        self.actual=ends.copy();self.next_time+=CONTROL_PERIOD_S
        return dict(time_s=float(time_s),period_s=CONTROL_PERIOD_S,heads_m=heads,
            start_opening_m=start,target_opening_m=targets,end_opening_m=ends,rate_m_per_s=rates,
            enabled=self.enabled.copy(),previous_enabled=before,modes=modes,
            latch_changed=[a!=b for a,b in zip(before,self.enabled)],
            dam_limit_exceeded_at_update=[h>s for h,s in zip(heads,c['unchanged_dam_safe_head_m'])],
            endpoint_is_plan_until_ODE_block_completes=True)
