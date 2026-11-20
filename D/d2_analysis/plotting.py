"""Plot existing Stage D.2 tables and saved trajectories; never run an ODE."""
from __future__ import annotations
import hashlib
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/cee690-d2-report-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
RESULTS=HERE/'results'
PLOTS=HERE/'plots'
DET='deterministic_september_2022'
COLORS={'wetland':'#177e73','no_wetland':'#a64b36'}
LABELS={'wetland':'Wetland B/B','no_wetland':'No wetland'}


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(fig,name):
    fig.savefig(PLOTS/name,dpi=180,bbox_inches='tight')
    plt.close(fig)
    print(name)


def ecdf(ax,values,*,label,color,ls='-',marker='o',bounds=None,lw=1.8):
    values=np.asarray(values,dtype=float)
    if not len(values) or not np.isfinite(values).all():raise ValueError('ECDF requires finite, nonempty saved metrics')
    unique,count=np.unique(values,return_counts=True)
    cumulative=np.cumsum(count)/len(values)
    if bounds is None:
        span=np.ptp(unique);pad=.06*span if span else max(abs(unique[0])*.05,.01)
        bounds=(unique[0]-pad,unique[-1]+pad)
    ax.step(np.r_[bounds[0],unique,bounds[1]],np.r_[0.,cumulative,1.],where='post',
            color=color,ls=ls,lw=lw,label=label)
    ax.plot(unique,cumulative,linestyle='none',marker=marker,ms=4,color=color,markerfacecolor='white')
    ax.set_xlim(*bounds);ax.set_ylim(-.025,1.025);ax.set_yticks([0,.25,.5,.75,1.])
    ax.set_ylabel('Empirical fraction')


def shared_bounds(*values):
    v=np.concatenate([np.asarray(a,dtype=float) for a in values]);span=np.ptp(v)
    pad=.06*span if span else max(abs(v[0])*.05,.01)
    return float(v.min()-pad),float(v.max()+pad)


def comparison_row(axes,pairs,metric,label,units):
    a,b=axes
    bounds=shared_bounds(pairs['wetland_'+metric],pairs['no_wetland_'+metric])
    for kind in ('wetland','no_wetland'):
        ecdf(a,pairs[kind+'_'+metric],label=f'{LABELS[kind]} (n={len(pairs)})',color=COLORS[kind],
             ls='-' if kind=='wetland' else '--',marker='o' if kind=='wetland' else 's',bounds=bounds)
    a.set_title(label);a.set_xlabel(f'{label} ({units})');a.legend(loc='lower right',fontsize=8)
    bounds=shared_bounds(pairs['difference_'+metric])
    populations=[('All valid',pairs,'#303030','-',2.1),
                 ('Copula/VAR',pairs[pairs.method=='copula_var'],'#5873b8','--',1.6),
                 ('Block bootstrap',pairs[pairs.method=='block_bootstrap'],'#b17a25',':',1.8)]
    for name,p,color,ls,lw in populations:
        if len(p):ecdf(b,p['difference_'+metric],label=f'{name} (n={len(p)})',color=color,ls=ls,bounds=bounds,lw=lw)
    difference_units=units.split(',')[0]
    b.set_title('Paired difference, including ties');b.set_xlabel(f'No wetland − wetland ({difference_units})')
    b.legend(loc='lower right',fontsize=8)


def primary_ecdf(pairs,cold):
    fig,ax=plt.subplots(2,2,figsize=(12,8.5))
    comparison_row(ax[0],pairs,'river_peak_m','Full-window river peak','m, model datum')
    comparison_row(ax[1],pairs,'outlet_peak_m3_s','Full-window outlet peak','m³/s')
    fig.suptitle('Primary peak ECDFs: complete 144 h window, including initialization',fontsize=14)
    fig.text(.5,.015,f'15 s output grid; {len(pairs)} valid warm pairs; {cold} original cold IDs remain outside the liquid-water domain.\n'
             'Right-continuous empirical steps; tied values retained. Pooled curves describe this pilot, not a flood probability.',
             ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.065,1,.965]);save(fig,'02_primary_peak_ecdfs.png')


def response_ecdf(pairs,cold):
    fig,ax=plt.subplots(3,2,figsize=(12,12))
    rows=[('maximum_rise_above_initial_m','Maximum rise above H₀','m'),
          ('initial_level_exceedance_h','Time strictly above H₀','h'),
          ('initial_discharge_excess_m3','Excess outlet volume above Q₀','m³')]
    for i,(name,label,unit) in enumerate(rows):comparison_row(ax[i],pairs,name,label,unit)
    fig.suptitle('Initial-reference response ECDFs — these references are not flood thresholds',fontsize=14)
    fig.text(.5,.012,'H₀ = 97.35 m; Q₀ ≈ 0.4444 m³/s (calculations use the exact saved initial values).\n'
             f'15 s grid; {len(pairs)} valid warm pairs; {cold} cold exclusions retained in accounting. No initial samples removed.',
             ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.055,1,.97]);save(fig,'03_initial_reference_ecdfs.png')


def load_trajectory(metrics,sid,kind):
    row=metrics[(metrics.dt_s==15)&(metrics.scenario_id==sid)&(metrics.case_name==kind)&(metrics.design=='best')]
    if len(row)!=1 or not bool(row.iloc[0]['valid']):raise ValueError(f'Expected one valid saved primary trajectory: {sid}, {kind}')
    source=ROOT/row.iloc[0]['source']
    if digest(source)!=row.iloc[0]['source_sha256']:raise ValueError('Diagnostic source JSON changed after metric extraction')
    with np.load(source.with_suffix('.npz'),allow_pickle=False) as z:return {name:z[name] for name in z.files}


def initial_diagnosis(metrics,weather):
    chosen=[('copula_var_00000','First Copula/VAR case'),
            ('block_bootstrap_00000','Noninitial peak example'),
            ('block_bootstrap_00002','Saved zero-rain case')]
    fig,ax=plt.subplots(3,3,figsize=(15,10),sharex=True)
    for j,(sid,title) in enumerate(chosen):
        source=ROOT/'D/g5_diagnostics/recovered/scenarios'/f'{sid}.csv'
        row=weather[weather.scenario_id==sid].iloc[0]
        if digest(source)!=row['forcing_sha256']:raise ValueError('Saved weather source hash changed')
        w=pd.read_csv(source);times=pd.to_datetime(w.timestamp_utc,utc=True)
        hours=(times-times.iloc[0]).dt.total_seconds()/3600
        if len(w)!=144 or not np.array_equal(hours.to_numpy(),np.arange(144)):raise ValueError('Expected original 144 hourly interval starts')
        ax[0,j].bar(hours,w.precipitation_mm_per_h,width=1,align='edge',color='#5c89b0',alpha=.85)
        ax[0,j].set_title(f'{title}\n{sid}',fontsize=11)
        ax[0,j].set_ylabel('Hourly precipitation (mm/h)');ax[0,j].set_ylim(bottom=0)
        ax[0,j].text(.98,.92,f"Total: {row['total_rain_mm']:.3g} mm",transform=ax[0,j].transAxes,ha='right',va='top',fontsize=9)
        for kind in ('wetland','no_wetland'):
            z=load_trajectory(metrics,sid,kind);t=z['time']/3600;h=z['downstream_stage']
            ax[1,j].plot(t,h,color=COLORS[kind],ls='-' if kind=='wetland' else '--',lw=1.5,label=LABELS[kind])
            drawup=h-np.minimum.accumulate(h)
            ax[2,j].plot(t,drawup,color=COLORS[kind],ls='-' if kind=='wetland' else '--',lw=1.5,label=LABELS[kind])
        ax[1,j].axhline(h[0],color='#666666',ls=':',label='Initial level H₀')
        ax[1,j].set_ylabel('River elevation (m)');ax[2,j].set_ylabel('Rise above running minimum (m)')
        ax[2,j].set_xlabel('Hours since synthetic-window start (UTC)')
        ax[2,j].set_xlim(0,144);ax[2,j].set_xticks([0,24,48,72,96,120,144])
    ax[1,0].legend(fontsize=8,loc='upper right')
    fig.suptitle('Why an initial peak can hide later response: all initial samples are retained',fontsize=14)
    fig.text(.5,.012,'Illustrative existing IDs, not a selected population for inference. Lower panels show the secondary running-minimum excursion,\n'
             'which differs from the requested rise above the initial level. Nonzero dry-case response is not rainfall-only runoff.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.06,1,.965]);save(fig,'01_initial_peak_diagnosis.png')


def convergence_plot(frame):
    d=frame[frame.comparable.eq(True)].copy()
    if d.empty:raise ValueError('No comparable completed resolution runs for a convergence figure')
    ids=([DET] if DET in set(d.scenario_id) else [])+sorted(set(d.scenario_id)-{DET})
    position={sid:i for i,sid in enumerate(ids)}
    labels=['Det.' if s==DET else s.replace('copula_var_','CV').replace('block_bootstrap_','BB').replace('0000','0') for s in ids]
    fig,ax=plt.subplots(2,2,figsize=(14,9.5),sharex=True)
    groups=[('wetland','reference','A/A reference','#777777','^'),
            ('wetland','best','Wetland B/B',COLORS['wetland'],'o'),
            ('no_wetland','best','No wetland',COLORS['no_wetland'],'s')]
    for kind,design,label,color,marker in groups:
        g=d[(d.case_name==kind)&(d.design==design)].copy();g['x']=g.scenario_id.map(position);g=g.sort_values('x')
        if g.empty:continue
        ax[0,0].plot(g.x,g.river_common_grid_max_abs_m,marker=marker,color=color,lw=.8,label=label,ms=4)
        ax[0,1].plot(g.x,g.fine_minus_coarse_initial_discharge_excess_m3.abs(),marker=marker,color=color,lw=.8,ms=4)
        ax[1,0].plot(g.x,g.coarse_balance_relative,marker=marker,color=color,ls='--',lw=.8,ms=4)
        ax[1,0].plot(g.x,g.fine_balance_relative,marker=marker,color=color,ls='-',lw=1,ms=4)
        ax[1,1].plot(g.x,g.latch_record_mismatch_count,marker=marker,color=color,lw=.8,ms=4)
    ax[0,0].set_title('River trajectory: maximum common-node difference')
    ax[0,0].set_ylabel('Absolute difference (m)');ax[0,0].legend(fontsize=8)
    ax[0,1].set_title('Initial-reference outlet excess volume')
    ax[0,1].set_ylabel('|7.5 s − 15 s| (m³)')
    for a in ax[0]:a.ticklabel_format(axis='y',style='sci',scilimits=(0,0));a.set_ylim(bottom=0)
    ax[1,0].set_title('Post-integration water-balance audit')
    ax[1,0].set_ylabel('Relative residual');ax[1,0].set_yscale('log')
    ax[1,0].axhline(1e-4,color='black',ls=':',label='Unchanged tolerance 1e−4')
    ax[1,0].plot([],[],color='#555555',ls='--',label='15 s');ax[1,0].plot([],[],color='#555555',ls='-',label='7.5 s')
    ax[1,0].legend(fontsize=8,loc='lower right')
    ax[1,1].set_title('Feedback latch disagreement at shared 120 s clock')
    ax[1,1].set_ylabel('Mismatched controller records')
    if d.latch_record_mismatch_count.max()==0:ax[1,1].set_ylim(-.05,.5)
    ax[1,1].text(.03,.92,f'Max actual-opening difference: {d.gate_common_grid_max_abs_m.max():.3g} m',
                 transform=ax[1,1].transAxes,fontsize=9,va='top')
    for a in ax[1]:a.set_xticks(range(len(ids)),labels,rotation=55,ha='right');a.set_xlabel('Original scenario ID (Det. = September event)')
    fig.suptitle('15 s versus 7.5 s: finite numerical and controller comparisons',fontsize=14)
    fig.text(.5,.012,'Only complete valid comparable trajectories are plotted; cold exclusions remain in the accounting tables.\n'
             'Lines connect case summaries, not time evolution. Small numerical differences are not physical calibration or a convergence proof.',
             ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.07,1,.96]);save(fig,'04_resolution_and_controller_checks.png')


def main():
    PLOTS.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.18})
    metrics=pd.read_csv(RESULTS/'metrics.csv');pairs=pd.read_csv(RESULTS/'paired_metrics.csv');weather=pd.read_csv(RESULTS/'weather_summary.csv')
    original=pairs[(pairs.dt_s==15)&(pairs.scenario_id!=DET)]
    if original.scenario_id.duplicated().any() or len(original)!=20:raise ValueError('Expected one retained primary row per original weather ID')
    valid=original[original.paired_valid.eq(True)]
    cold=int((original.wetland_classification=='meteorological_liquid_water_applicability').sum())
    initial_diagnosis(metrics,weather);primary_ecdf(valid,cold);response_ecdf(valid,cold)
    f=RESULTS/'convergence.csv'
    if f.exists():convergence_plot(pd.read_csv(f))
    else:print('04 pending: convergence.csv does not yet exist; no placeholder figure generated')


if __name__=='__main__':main()
