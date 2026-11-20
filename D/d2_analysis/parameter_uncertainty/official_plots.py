"""Guarded access to unchanged official plots for a configured parameter MCS.

No physical model, sampler or fitting routine is executed by this module.
A deterministic factorial sensitivity grid is not an iid Monte Carlo sample.
"""
from __future__ import annotations

import hashlib
import importlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import sys
import warnings

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OFFICIAL_ROOT=ROOT/'A/work/multivarious'
OFFICIAL_REVISION='06d7aed5b541ec7fe26a6dc69be6aa7fb610aa02'
OFFICIAL_HASHES={
    'multivarious/utl/plot_ECDF_ci.py':'3fcfcf204fb9b01793c637247fe0dc4bbb3ccb1bf7ef97fbcd0311e82ce64b35',
    'multivarious/utl/plot_scatter_dist.py':'e9117b3766c6b59419f4eb1488b7cc84f34802f70a73742bb5f2cf7fb2947791',
    'multivarious/utl/format_plot.py':'83eb7c26465ba076057b7d3f0cf17330c99080e5c3532dbd65fefb6670bd8db2',
}


def _sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _clean(value):
    import numpy as np
    if isinstance(value,dict):return {str(k):_clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [_clean(v) for v in value]
    if isinstance(value,np.generic):return _clean(value.item())
    if isinstance(value,float) and not np.isfinite(value):return None
    return value


def _official_functions():
    for name,expected in OFFICIAL_HASHES.items():
        if _sha(OFFICIAL_ROOT/name)!=expected:raise ValueError(f'Official source differs from reviewed revision: {name}')
    if str(OFFICIAL_ROOT) not in sys.path:sys.path.insert(0,str(OFFICIAL_ROOT))
    a=importlib.import_module('multivarious.utl.plot_ECDF_ci')
    b=importlib.import_module('multivarious.utl.plot_scatter_dist')
    for m in (a,b):
        if not Path(m.__file__).resolve().is_relative_to(OFFICIAL_ROOT.resolve()):
            raise ValueError('Another multivarious installation is already loaded; use a fresh process with the pinned checkout')
    return a.plot_ECDF_ci,b.plot_scatter_dist


def render_official_plots(
    results_csv, parameter_columns, output_columns, output_dir=None, *,
    iid_parameter_draws_confirmed=False, provenance=None,
    valid_column='paired_valid', sample_id_column='sample_id',
    failure_reason_columns=('failure_reason',), allow_conditional_valid_subset=False,
    ecdf_confidence_percent=95.0, correlation_confidence=0.95,
    variable_labels=None, scatter_output_columns=None,
):
    """Render official plots from one real paired-result row per parameter draw.

    Required provenance fields:
      study_kind='fixed_storm_parameter_uncertainty'
      sampling_kind='iid_monte_carlo'
      distribution_status in {'explicit_exploratory_assumptions','documented'}
      distribution_config_sha256, seed, sample_count

    Every planned draw, including invalid ones, must have a unique CSV row.
    Invalid rows require explicit permission to plot a conditional valid subset
    and nonempty diagnostics in the supplied failure-reason columns. Unsafe
    but valid results and negative effects are retained without filtering.

    Returns metadata and writes it to official_plot_metadata.json. No CLI auto-
    execution is provided because real configured physical results are required.
    """
    import numpy as np
    import pandas as pd
    p=dict(provenance or {})
    required=dict(study_kind='fixed_storm_parameter_uncertainty',sampling_kind='iid_monte_carlo')
    if iid_parameter_draws_confirmed is not True:raise ValueError('Configured iid parameter draws must be explicitly confirmed')
    for key,value in required.items():
        if p.get(key)!=value:raise ValueError(f'{key} must equal {value!r}; weather pilots and factorial grids are not parameter MCS data')
    if p.get('distribution_status') not in {'explicit_exploratory_assumptions','documented'}:
        raise ValueError('A resolved, explicitly labeled distribution specification is required')
    if not re.fullmatch(r'[0-9a-fA-F]{64}',str(p.get('distribution_config_sha256',''))):
        raise ValueError('Provide the SHA-256 of the resolved distribution configuration')
    if p.get('seed') is None:raise ValueError('Provide the recorded parameter-sampling seed provenance')
    n=p.get('sample_count')
    if isinstance(n,bool) or not isinstance(n,int) or n<4:raise ValueError('sample_count must be an integer >=4 for official correlation intervals')
    if not (0<float(ecdf_confidence_percent)<100) or not (0<float(correlation_confidence)<1):
        raise ValueError('ECDF confidence uses percent; correlation confidence uses a fraction')
    parameters=list(parameter_columns);outputs=list(output_columns)
    if len(parameters)!=3 or len(set(parameters))!=3 or not outputs or len(set(outputs))!=len(outputs):
        raise ValueError('Supply exactly three distinct parameter columns and distinct nonempty output columns')
    if set(parameters)&set(outputs):raise ValueError('Parameter and output column names must be distinct')
    scatter_outputs=outputs if scatter_output_columns is None else list(scatter_output_columns)
    if len(set(scatter_outputs))!=len(scatter_outputs) or not set(scatter_outputs)<=set(outputs):
        raise ValueError('Scatter output columns must be a distinct subset of the ECDF output columns')
    source=Path(results_csv).resolve();frame=pd.read_csv(source)
    needed=[sample_id_column,valid_column,*parameters,*outputs]
    missing=[name for name in needed if name not in frame]
    if missing:raise ValueError(f'Missing required result columns: {missing}')
    if len(frame)!=n or frame[sample_id_column].isna().any() or frame[sample_id_column].duplicated().any():
        raise ValueError('Retain exactly one unique record for every planned parameter draw')
    status=frame[valid_column]
    if not status.isin([True,False]).all():raise ValueError('Validity column must contain explicit Boolean values for all draws')
    valid=status.eq(True);rejected=frame[~valid]
    if len(rejected):
        if not allow_conditional_valid_subset:raise ValueError('Invalid draws exist: explicit conditional-subset accounting is required before plotting')
        reasons=[c for c in failure_reason_columns if c in frame]
        if not reasons:raise ValueError('Invalid draws require recorded failure reasons')
        has_reason=rejected[reasons].fillna('').astype(str).apply(lambda c:c.str.strip().ne('')).any(axis=1)
        if not has_reason.all():raise ValueError('Every invalid draw needs a recorded cause')
    else:reasons=[]
    chosen=frame[valid]
    if len(chosen)<4:raise ValueError('Fewer than four valid paired parameter draws; official Fisher intervals are unavailable')
    if not np.isfinite(frame[parameters].to_numpy(float)).all():raise ValueError('Every sampled parameter vector must remain finite and recorded')
    if not np.isfinite(chosen[outputs].to_numpy(float)).all():raise ValueError('A valid row contains missing/nonfinite physical outputs')
    destination=Path(output_dir).resolve() if output_dir is not None else HERE/'plots'
    destination.mkdir(parents=True,exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR','/private/tmp/cee690-parameter-official-plots-mpl')
    sys.dont_write_bytecode=True
    import matplotlib
    matplotlib.use('Agg',force=True)
    import matplotlib.pyplot as plt
    official_ecdf,official_scatter=_official_functions()
    labels=dict(variable_labels or {})
    metadata=dict(provenance=p,input_csv=str(source),input_csv_sha256=_sha(source),
        official_revision=OFFICIAL_REVISION,official_sources=OFFICIAL_HASHES,
        total_draws=n,valid_draws=len(chosen),invalid_draws=len(rejected),
        invalid_records=rejected[[sample_id_column,*reasons]].to_dict('records'),
        conditioning='Complete valid paired physical trajectories' if len(rejected) else 'All configured parameter draws',
        safety_or_benefit_sign_filter_applied=False,
        ecdf_method='Official rank/(N+1) plotting positions; active pointwise Beta order-statistic limits; not simultaneous, not ties-aggregated count/N ECDF',
        ci_scope='Conditional on declared sampling laws; not calibrated parameter uncertainty or flood probability',
        ecdf_confidence_percent=float(ecdf_confidence_percent),correlation_confidence=float(correlation_confidence),
        scatter_method='Pearson correlation with approximate Fisher-z intervals; not causal or variance-based sensitivity',
        dependencies={name:importlib.metadata.version(name) for name in ('numpy','scipy','matplotlib')},
        scatter_output_columns=scatter_outputs,outputs={},plots=[],warnings=[])
    was_interactive=plt.isinteractive()
    try:
        with plt.rc_context({'text.usetex':False}):
            for column in outputs:
                values=chosen[column].to_numpy(float);unique=np.unique(values)
                stem=re.sub(r'[^A-Za-z0-9_.-]+','_',column).strip('_')
                fig=plt.figure();number=fig.number
                try:
                    with warnings.catch_warnings(record=True) as caught:
                        warnings.simplefilter('always')
                        x,F,lower,upper=official_ecdf(values,float(ecdf_confidence_percent),number,
                                                     x_label=labels.get(column,column),save_plots=False)
                    metadata['warnings'].extend(str(w.message) for w in caught)
                    fig=plt.figure(number)
                    for annotation in fig.axes[0].texts:
                        if annotation.get_text().startswith('$x_{med}$'):
                            annotation.set_text(rf'$x_{{med}}$ = {np.median(values):.3f}')
                    fig.axes[0].set_title(f'Official plotting positions with {ecdf_confidence_percent:g}% pointwise Beta limits',fontsize=12)
                    tie_note='Ties present: exact continuous-order-statistic coverage is not asserted.' if len(unique)<len(values) else 'Continuous i.i.d. interpretation is conditional on the assumed sampling model.'
                    note=f'Rank/(N+1); pointwise, not simultaneous. Valid n={len(values)} of {n} draws.\n{tie_note}\nConventional sample median={np.median(values):.6g}; index-based median annotation corrected for display only.'
                    if len(rejected):note+='\nConditional valid subset; rejected-draw causes remain in the metadata.'
                    fig.text(.5,.015,note,ha='center',va='bottom',fontsize=8)
                    fig.subplots_adjust(bottom=.28,top=.86)
                    image_path=destination/f'official_ecdf_{stem}.png';fig.savefig(image_path,dpi=180,bbox_inches='tight')
                    table_path=destination/f'official_ecdf_{stem}.csv'
                    pd.DataFrame({'sorted_value':x,'official_plotting_position':F,'pointwise_lower':lower,'pointwise_upper':upper}).to_csv(table_path,index=False)
                    metadata['plots'].append(str(image_path))
                    metadata['outputs'][column]=dict(n=len(values),unique_values=len(unique),
                        negative_count=int((values<0).sum()),zero_count=int((values==0).sum()),positive_count=int((values>0).sum()),
                        conventional_median=float(np.median(values)),official_index_median=float(np.sort(values)[int(round(len(values)/2))]),median_annotation_corrected=True,
                        ties_present=len(unique)<len(values),continuous_exact_beta_coverage_asserted=False,
                        returned_array_csv=str(table_path))
                finally:plt.close(number)
            varying=[name for name in parameters+scatter_outputs if np.ptp(chosen[name].to_numpy(float))>0]
            excluded=[name for name in parameters+scatter_outputs if name not in varying]
            xnames=[name for name in parameters if name in varying];ynames=[name for name in scatter_outputs if name in varying]
            metadata['constant_columns_omitted_from_correlation_only']=excluded
            if len(varying)>=2:
                X=chosen[xnames].to_numpy(float).T if xnames else None
                Y=chosen[ynames].to_numpy(float).T if ynames else None
                names={'X':[labels.get(c,c) for c in xnames],'Y':[labels.get(c,c) for c in ynames]}
                fig=plt.figure();number=fig.number;plt.close(number)
                try:
                    with warnings.catch_warnings(record=True) as caught:
                        warnings.simplefilter('always')
                        fig=official_scatter(X,Y,fig_num=number,var_names=names,n_bins=20,font_size=10,ci=float(correlation_confidence))
                    metadata['warnings'].extend(str(w.message) for w in caught)
                    caption=f'Official Pearson/Fisher-z diagnostic; valid n={len(chosen)} of {n}.\nIntervals are approximate; not simultaneous or causal.'
                    if excluded:caption+='\nConstant columns omitted only from correlations: '+', '.join(excluded)+'. Their results remain in the CSV/ECDF outputs.'
                    if len(rejected):caption+='\nConditional valid subset; rejected-draw causes remain in metadata.'
                    fig.text(.5,.015,caption,ha='center',va='bottom',fontsize=9)
                    fig.subplots_adjust(left=.10,right=.98,bottom=.17,top=.98,hspace=.14,wspace=.14)
                    path=destination/'official_parameter_output_scatter.png';fig.savefig(path,dpi=150,bbox_inches='tight');metadata['plots'].append(str(path))
                finally:plt.close(number)
            else:metadata['scatter_skipped_reason']='Fewer than two nonconstant variables; correlation is undefined'
    finally:
        if not was_interactive:plt.ioff()
    metadata=_clean(metadata)
    (destination/'official_plot_metadata.json').write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n')
    return metadata
