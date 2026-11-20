"""Plot completed paired parameter MC results with the inspected official utilities.

This script only reads saved results. It never samples, solves an ODE or fits a
surrogate. The deterministic factorial grid cannot be used as an iid ensemble.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from D.d2_analysis.parameter_uncertainty.official_plots import render_official_plots

FACTORS = ['m_K_u', 'm_routing', 'm_R_r']
OUTCOMES = [
    'wetland_peak_stage_m', 'no_wetland_peak_stage_m', 'delta_peak_stage_m',
    'wetland_peak_discharge_m3_s', 'no_wetland_peak_discharge_m3_s',
    'delta_peak_discharge_m3_s',
]
LABELS = {
    'm_K_u': r'$m_{K_u}$',
    'm_routing': r'$m_{route}$',
    'm_R_r': r'$m_{R_r}$',
    'wetland_peak_stage_m': r'Wetland river peak $H_{r,max}$ (m)',
    'no_wetland_peak_stage_m': r'No-wetland river peak $H_{r,max}$ (m)',
    'delta_peak_stage_m': r'$\Delta H_{r,max}$ (m)',
    'wetland_peak_discharge_m3_s': r'Wetland outlet peak $Q_{out,max}$ (m$^3$/s)',
    'no_wetland_peak_discharge_m3_s': r'No-wetland outlet peak $Q_{out,max}$ (m$^3$/s)',
    'delta_peak_discharge_m3_s': r'$\Delta Q_{out,max}$ (m$^3$/s)',
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def inspect_saved_inputs(results_csv, config_path, plan_path):
    """Verify fixed approved draw identity before invoking any plotting routine."""
    config = json.loads(Path(config_path).read_text())
    plan = json.loads(Path(plan_path).read_text())
    digest = sha(config_path)
    if digest != plan['config_sha256'] or config != plan['config']:
        raise ValueError('Configuration differs from the immutable pre-run sampling plan')
    settings = config['monte_carlo']
    if settings['sampling_kind'] != 'iid_monte_carlo':
        raise ValueError('Only the explicitly configured iid parameter ensemble is eligible')
    frame = pd.read_csv(results_csv)
    required = ['sample_id', 'study_kind', 'dt_s', 'paired_valid', 'failure_reason',
                *FACTORS, *OUTCOMES, 'K_u', 'k_ud_1', 'k_ud_2', 'k_ur', 'R_r']
    missing = [name for name in required if name not in frame]
    if missing:
        raise ValueError('Missing saved-result columns: '+', '.join(missing))
    planned = pd.DataFrame(plan['mc']).set_index('sample_id')
    if len(frame) != settings['sample_count'] or len(planned) != settings['sample_count']:
        raise ValueError('All planned MC draws must have a saved result row before plotting')
    if frame.sample_id.duplicated().any() or set(frame.sample_id) != set(planned.index):
        raise ValueError('Saved sample IDs differ from the configured MC draw IDs')
    if not frame.study_kind.eq(settings['study_kind']).all():
        raise ValueError('Factorial, numerical pilot or weather rows are not MC parameter draws')
    dt = frame.dt_s.to_numpy(float)
    if not np.isfinite(dt).all() or not np.all(dt == dt[0]) or dt[0] <= 0:
        raise ValueError('Each sample must use the same recorded integration resolution; do not pool repeats')
    for name in FACTORS:
        expected = planned.loc[frame.sample_id, name].to_numpy(float)
        if not np.allclose(frame[name].to_numpy(float), expected, rtol=2e-14, atol=0):
            raise ValueError('Saved multiplier differs from its immutable draw: '+name)
    for row in frame.itertuples(index=False):
        physical = planned.loc[row.sample_id, 'physical_values']
        expected = dict(K_u=physical['K_u'], k_ud_1=physical['k_ud'][0],
                        k_ud_2=physical['k_ud'][1], k_ur=physical['k_ur'], R_r=physical['R_r'])
        for name, value in expected.items():
            if not np.isclose(getattr(row, name), value, rtol=2e-14, atol=0):
                raise ValueError(f'{row.sample_id}: saved physical mapping changed for {name}')
    valid = frame.paired_valid
    if not valid.isin([True, False]).all():
        raise ValueError('Every draw needs an explicit paired_valid result')
    for suffix in ['peak_stage_m', 'peak_discharge_m3_s']:
        expected = frame.loc[valid, 'no_wetland_'+suffix] - frame.loc[valid, 'wetland_'+suffix]
        actual = frame.loc[valid, 'delta_'+suffix]
        if not np.isfinite(actual).all() or not np.allclose(actual, expected, rtol=1e-10, atol=1e-13):
            raise ValueError('Signed paired difference must equal no-wetland minus wetland: '+suffix)
    provenance = {
        'study_kind': config['study'],
        'sampling_kind': settings['sampling_kind'],
        'distribution_status': config['distribution_status'],
        'distribution_config_sha256': digest,
        'seed': settings['seed'],
        'sample_count': settings['sample_count'],
        'sampling_plan_sha256': sha(plan_path),
        'sampling_plan_verified': True,
        'dt_s': float(dt[0]),
        'distribution': settings['multiplier_distribution'],
        'support': [settings['lower_multiplier'], settings['upper_multiplier']],
        'interpretation': config['interpretation'],
        'paired_difference': 'no-wetland minus wetland; negative values are retained',
        'fixed_weather': config['fixed_weather'],
        'flood_threshold_m': config['outcomes']['flood_threshold_m'],
        'surrogate_fit_performed': False,
    }
    return config, provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=HERE/'results/monte_carlo/paired_results.csv')
    parser.add_argument('--config', type=Path, default=HERE/'config.json')
    parser.add_argument('--plan', type=Path, default=HERE/'results/study_plan.json')
    parser.add_argument('--output', type=Path, default=HERE/'plots/official_mc')
    parser.add_argument('--allow-conditional-valid-subset', action='store_true',
                        help='Explicitly expose invalid-draw accounting and label all plots as conditional')
    args = parser.parse_args()
    config, provenance = inspect_saved_inputs(args.results, args.config, args.plan)
    settings = config['monte_carlo']
    metadata = render_official_plots(
        args.results, FACTORS, OUTCOMES, args.output,
        iid_parameter_draws_confirmed=True, provenance=provenance,
        allow_conditional_valid_subset=args.allow_conditional_valid_subset,
        ecdf_confidence_percent=settings['ecdf_confidence_percent'],
        correlation_confidence=settings['correlation_confidence'],
        variable_labels=LABELS,
        scatter_output_columns=['delta_peak_stage_m', 'delta_peak_discharge_m3_s'],
    )
    captions = '''# Official parameter-uncertainty plots

These plots use the saved full-ODE results from the approved iid exploratory
parameter study. They do not use the 27-point factorial grid, the earlier
16-weather-scenario pilot, a surrogate fit or a flood threshold. Positive paired
differences mean no-wetland minus wetland. Negative differences remain results.

The six distribution plots show wetland, no-wetland and paired full-window peak
stage and outlet discharge. The unchanged official plotting utility uses
rank/(N+1) plotting positions and 95% pointwise Beta order-statistic limits,
not a ties-aggregated count/N ECDF or a simultaneous confidence band. Tied/atomic
outputs do not have the usual exact continuous-order-statistic coverage claim.
A conventional sample median is corrected for display only; returned plotting
positions and limits are unchanged. The coefficient-of-variation annotation is
not informative for a signed effect with a near-zero mean.

The five-variable scatter figure contains the three physical-parameter
multipliers and the two signed paired peak differences. Its upper panels give
Pearson correlations and approximate Fisher-z limits (upper limit displayed
above, lower limit below). These are association diagnostics, not Sobol indices,
causal effects or calibrated local engineering risk. Constant columns, if any,
are omitted only from correlation panels and remain in CSV/distribution outputs.

All statements are conditional on the declared exploratory log-uniform sampling
laws, fixed storm, geometry, feedback policy and initial conditions. The saved
metadata records the exact planned/valid/invalid denominator, any rejected-draw
causes, data/configuration hashes, official revision and plotting warnings.
Unsafe but physically valid runs are not removed. The distinct initial water
volumes and realization-specific outlet resistance remain part of the accepted
paired structural comparison.
'''
    (args.output/'figure_captions.md').write_text(captions)
    print(json.dumps({key:metadata[key] for key in ('total_draws','valid_draws','invalid_draws','plots')}, indent=2))


if __name__ == '__main__':
    main()
