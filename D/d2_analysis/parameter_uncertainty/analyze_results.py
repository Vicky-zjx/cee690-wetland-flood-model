"""Summarize saved ODE evaluations; never fit or run a physical model.

The deterministic factorial and exploratory iid Monte Carlo remain separate.
Positive paired differences mean no-wetland minus wetland. No flood threshold
is supplied or inferred by this analysis.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/cee690-parameter-analysis-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
FACTORS = ['m_K_u', 'm_routing', 'm_R_r']
LABELS = ['Soil infiltration multiplier', 'Common routing-rate multiplier', 'River resistance multiplier']
METRICS = {
    'peak_stage_m': 'm',
    'peak_discharge_m3_s': 'm3/s',
    'maximum_rise_above_initial_m': 'm',
    'initial_level_exceedance_h': 'h',
    'initial_discharge_excess_m3': 'm3',
    'total_outlet_volume_m3': 'm3',
}


def boolean(series):
    values = series.astype(str).str.lower()
    if not values.isin(['true', 'false']).all():
        raise ValueError('Missing or ambiguous Boolean classification')
    return values.eq('true')


def clean(value):
    if isinstance(value, dict): return {k: clean(v) for k, v in value.items()}
    if isinstance(value, list): return [clean(v) for v in value]
    if isinstance(value, np.generic): value = value.item()
    if isinstance(value, float) and not np.isfinite(value): return None
    return value


def load(study, count):
    path = HERE / 'results' / study / 'paired_results.csv'
    frame = pd.read_csv(path)
    if len(frame) != count or frame.sample_id.duplicated().any():
        raise AssertionError(f'{study}: expected exactly {count} uniquely identified planned rows')
    if not np.allclose(frame.dt_s, 15):
        raise AssertionError('Do not pool numerical refinement repeats with study samples')
    if not np.allclose(frame.K_u, frame.m_K_u * 1e-6, rtol=1e-12, atol=0):
        raise AssertionError('Infiltration mapping changed')
    for key, base in [('k_ud_1', 1/28800), ('k_ud_2', 1/36000), ('k_ur', 1/72000)]:
        if not np.allclose(frame[key], frame.m_routing * base, rtol=1e-12, atol=0):
            raise AssertionError('Shared routing mapping changed: '+key)
    if not np.allclose(frame.R_r, frame.m_R_r * .45, rtol=1e-12, atol=0):
        raise AssertionError('Resistance mapping changed')
    if (boolean(frame.paired_valid) & ~(boolean(frame.wetland_valid) & boolean(frame.no_wetland_valid))).any():
        raise AssertionError('A valid pair includes an invalid model outcome')
    valid = boolean(frame.paired_valid)
    for metric in METRICS:
        name = 'delta_' + metric
        frame[name] = frame['no_wetland_'+metric] - frame['wetland_'+metric]
        frame.loc[~valid, name] = np.nan
    for model in ('wetland', 'no_wetland'):
        eligible = boolean(frame[model+'_valid'])
        h = frame.loc[eligible, model+'_peak_stage_m']
        h0 = frame.loc[eligible, model+'_initial_stage_m']
        if not np.isfinite(h).all() or np.any(h < h0-1e-12):
            raise AssertionError('Full-horizon peak missing or initial state omitted')
        # Fixed accepted boundary h_d=97.15 m, B=0; outlet cap=3 m3/s.
        expected_q = np.minimum(3.0, np.maximum(h-97.15, 0)/frame.loc[eligible, 'R_r'])
        if not np.allclose(frame.loc[eligible, model+'_peak_discharge_m3_s'], expected_q, rtol=0, atol=1e-10):
            raise AssertionError('Peak discharge disagrees with the unchanged sampled-resistance rating law')
    return frame


def summarize(frame, study):
    valid = boolean(frame.paired_valid)
    rows = []
    outcomes = []
    for model in ('wetland', 'no_wetland'):
        for category, group in frame.groupby(model+'_classification', dropna=False):
            outcomes.append(dict(study=study, model=model, classification=category, count=len(group), planned=len(frame)))
    columns = [prefix+'_'+metric for metric in METRICS for prefix in ('wetland', 'no_wetland', 'delta')]
    for column in columns:
        # A common complete-pair denominator is used even for model marginals.
        values = frame.loc[valid, column].to_numpy(float)
        if not np.isfinite(values).all():
            raise AssertionError('A valid pair has a missing outcome: '+column)
        row = dict(study=study, metric=column, planned_pairs=len(frame), valid_pairs=int(valid.sum()),
                   missing_pairs=int((~valid).sum()), minimum=float(values.min()) if len(values) else None,
                   maximum=float(values.max()) if len(values) else None,
                   denominator='All complete valid pairs, including unsafe and negative-effect pairs')
        if study == 'monte_carlo' and len(values):
            row.update(mean=float(values.mean()), median=float(np.median(values)),
                       sample_sd=float(values.std(ddof=1)), mean_mc_standard_error=float(values.std(ddof=1)/np.sqrt(len(values))),
                       q05=float(np.quantile(values, .05)), q95=float(np.quantile(values, .95)))
        if column.startswith('delta_'):
            row.update(positive_count=int((values>0).sum()), zero_count=int((values==0).sum()), negative_count=int((values<0).sum()))
        rows.append(row)
    audit = dict(planned_pairs=len(frame), valid_pairs=int(valid.sum()), invalid_pairs=int((~valid).sum()),
                 unsafe_valid_pairs=int((valid & (boolean(frame.wetland_unsafe.fillna(False)) | boolean(frame.no_wetland_unsafe.fillna(False)))).sum()),
                 failure_sample_ids=frame.loc[~valid, 'sample_id'].tolist(),
                 unsafe_sample_ids=frame.loc[valid & (boolean(frame.wetland_unsafe.fillna(False)) | boolean(frame.no_wetland_unsafe.fillna(False))), 'sample_id'].tolist(),
                 maximum_paired_gate_difference_m=float(frame.maximum_paired_gate_opening_difference_m.max()),
                 paired_latch_mismatch_count=int(frame.paired_latch_mismatch_count.sum()),
                 all_valid_pairs_have_identical_weather=bool(boolean(frame.loc[valid, 'paired_weather_inputs_exact']).all()),
                 maximum_relative_mass_residual=max(float(frame.wetland_balance_relative.max()),float(frame.no_wetland_balance_relative.max())),
                 minimum_wetland_safety_margin_m=float(frame.wetland_minimum_safety_margin_m.min()),
                 minimum_no_wetland_safety_margin_m=float(frame.no_wetland_minimum_safety_margin_m.min()),
                 wetland_initial_peak_count=int((valid & frame.wetland_peak_stage_time_s.eq(0)).sum()),
                 no_wetland_initial_peak_count=int((valid & frame.no_wetland_peak_stage_time_s.eq(0)).sum()),
                 wetland_outlet_cap_peak_count=int((valid & frame.wetland_peak_discharge_m3_s.ge(3.0)).sum()),
                 no_wetland_outlet_cap_peak_count=int((valid & frame.no_wetland_peak_discharge_m3_s.ge(3.0)).sum()),
                 flood_threshold_status='Undefined: no accepted downstream flood threshold',
                 initial_reference_warning='Initial-stage/discharge exceedance is a response benchmark, not a flood exceedance statistic')
    return rows, outcomes, audit


def sensitivity_figures(frame, plots, tables):
    rows = []
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), sharey=True, layout='constrained')
    colors = ['#005A9C', '#E17C05', '#2C8C68']
    for ax, (resistance, group) in zip(axes, frame.groupby('m_R_r')):
        for color, (infiltration, part) in zip(colors, group.groupby('m_K_u')):
            part = part.sort_values('m_routing')
            valid = boolean(part.paired_valid)
            y = np.where(valid, part.delta_peak_stage_m * 1000, np.nan)
            ax.plot(part.m_routing, y, 'o-', color=color, label=f'Infiltration ×{infiltration:g}')
        ax.set(title=f'River resistance ×{resistance:g}', xlabel='Common routing-rate multiplier', xticks=[.5,1,2])
        ax.axhline(0, color='0.5', lw=.8); ax.grid(alpha=.2)
    axes[0].set_ylabel('Peak-stage reduction (mm)\nno wetland − wetland')
    axes[-1].legend(fontsize=8)
    fig.suptitle('27 deterministic sensitivity configurations — no probability weights')
    fig.savefig(plots/'factorial_peak_reduction.png', dpi=180); plt.close(fig)
    for factor in FACTORS:
        other = [x for x in FACTORS if x != factor]
        part = frame.loc[np.isclose(frame[other[0]], 1) & np.isclose(frame[other[1]], 1)].sort_values(factor)
        for _, row in part.iterrows():
            rows.append(dict(varied_factor=factor, factor=row[factor], sample_id=row.sample_id,
                paired_valid=row.paired_valid, wetland_peak_m=row.wetland_peak_stage_m,
                no_wetland_peak_m=row.no_wetland_peak_stage_m, reduction_mm=1000*row.delta_peak_stage_m,
                discharge_reduction_m3_s=row.delta_peak_discharge_m3_s))
    pd.DataFrame(rows).to_csv(tables/'one_factor_at_a_time.csv', index=False)


def mc_figure(frame, plots):
    valid = boolean(frame.paired_valid)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), sharey=True, layout='constrained')
    unsafe = boolean(frame.wetland_unsafe.fillna(False)) | boolean(frame.no_wetland_unsafe.fillna(False))
    for ax, column, label in zip(axes, FACTORS, LABELS):
        ax.scatter(frame.loc[valid & ~unsafe, column], 1000*frame.loc[valid & ~unsafe, 'delta_peak_stage_m'],
                   s=24, alpha=.8, color='#005A9C', label='Both safety feasible')
        if (valid & unsafe).any():
            ax.scatter(frame.loc[valid & unsafe, column], 1000*frame.loc[valid & unsafe, 'delta_peak_stage_m'],
                       s=35, marker='x', color='#C44E52', label='At least one unsafe')
        ax.set(xlabel=label, xlim=(.45,2.05)); ax.axhline(0, color='0.5', lw=.8); ax.grid(alpha=.2)
    axes[0].set_ylabel('Peak-stage reduction (mm)\nno wetland − wetland')
    axes[-1].legend(fontsize=8)
    fig.suptitle(f'{int(valid.sum())}/{len(frame)} complete pairs — exploratory parameter law, fixed 2022 storm')
    fig.savefig(plots/'mc_parameter_effects.png', dpi=180); plt.close(fig)


def detailed_audit(study, frame):
    """Retain every model diagnosis and every applicable upstream safety check."""
    limits = {'H_d1': 101.2, 'H_d2': 101.0, 'H_w1': 98.9, 'H_w2': 98.8}
    safety, trajectories = [], []
    for sample in frame.sample_id:
        for model in ('wetland', 'no_wetland'):
            path = HERE/'results'/study/f'{sample}__{model}__15s.json'
            d = json.loads(path.read_text())
            obs = d.get('all_rhs_observations', {})
            nodes = {r['variable']: r for r in d.get('safety', [])}
            rhs = {r['variable']: r for r in obs.get('safety_successful_rhs', [])}
            for variable, threshold in limits.items():
                applicable = not (model == 'no_wetland' and variable.startswith('H_w'))
                records = [source[variable] for source in (nodes,rhs) if variable in source]
                if any(r['threshold_m'] != threshold for r in records):
                    raise AssertionError('A physical safety threshold changed')
                peak = max((r['peak_m'] for r in records), default=np.nan)
                if applicable and d.get('full_valid_trajectory_available') and len(records) != 2:
                    raise AssertionError('A valid run lacks a node or RHS safety check')
                safety.append(dict(study=study, sample_id=sample, model=model, variable=variable,
                    applicable=applicable, threshold_m=threshold if applicable else None,
                    peak_over_nodes_and_successful_RHS_m=peak if applicable else None,
                    margin_m=threshold-peak if applicable else None,
                    exceeded=bool(peak>threshold) if applicable and records else None,
                    full_horizon_completed=d.get('full_integration_horizon_completed'),
                    classification=d['classification']))
            ratios = [r['minimum_open_gate_H_over_a'] for r in d.get('actual_gate_successful_rhs_observations', [])
                      if r.get('minimum_open_gate_H_over_a') is not None]
            trajectories.append(dict(study=study, sample_id=sample, model=model,
                classification=d['classification'], full_horizon_completed=d.get('full_integration_horizon_completed'),
                valid=d.get('full_valid_trajectory_available'), unsafe=d.get('safety_violation'),
                failure_reason=d.get('failure_reason'), failure_message=d.get('failure_message'),
                failure_stage=d.get('failure_stage'), first_failure_time_s=d.get('first_failure_time_s'),
                minimum_open_gate_H_over_a=min(ratios) if ratios else None,
                regimes=json.dumps(d.get('hydraulic_diagnostics',{}).get('regime_rhs_counts',{}), sort_keys=True),
                relative_mass_residual=d.get('balance_relative'),
                relative_mass_tolerance=d.get('balance_relative_tolerance'),
                maximum_derivative_mass_identity_m3_s=obs.get('derivative_physical_mass_identity_max_abs_m3_s'),
                solver_termination=d.get('solver_termination'), source=str(path.relative_to(HERE))))
    return safety, trajectories


def main():
    tables, plots = HERE/'tables', HERE/'plots'
    tables.mkdir(exist_ok=True); plots.mkdir(exist_ok=True)
    all_stats, all_outcomes, all_safety, all_trajectories, audit = [], [], [], [], {}
    for study, count in [('sensitivity', 27), ('monte_carlo', 64)]:
        frame = load(study, count)
        stats, outcomes, audit[study] = summarize(frame, study)
        all_stats.extend(stats); all_outcomes.extend(outcomes)
        safety, trajectories = detailed_audit(study, frame)
        all_safety.extend(safety); all_trajectories.extend(trajectories)
        frame.to_csv(tables/(study+'_paired_metrics.csv'), index=False)
        if study == 'sensitivity': sensitivity_figures(frame, plots, tables)
        else: mc_figure(frame, plots)
    pd.DataFrame(all_stats).to_csv(tables/'outcome_summaries.csv', index=False)
    pd.DataFrame(all_outcomes).to_csv(tables/'classification_counts.csv', index=False)
    pd.DataFrame(all_safety).to_csv(tables/'safety_checks.csv', index=False)
    pd.DataFrame(all_trajectories).to_csv(tables/'trajectory_audit.csv', index=False)
    (tables/'analysis_audit.json').write_text(json.dumps(clean(audit), indent=2, allow_nan=False)+'\n')
    print(json.dumps(clean(audit), indent=2, allow_nan=False))


if __name__ == '__main__': main()
