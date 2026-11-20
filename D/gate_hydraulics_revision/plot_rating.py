"""Compare the original and candidate ratings; no ODE or meteorology changes."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/cee690-stage-d-revision-matplotlib')
from pathlib import Path
from math import sqrt
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from B.gate_hydraulics import GateParameters, sluice_gate_flow
from .hydraulics import Assumptions, gate_flow


def main():
    out=Path(__file__).resolve().parent/'plots';out.mkdir(exist_ok=True)
    p=GateParameters(gravity=9.80665,gate_width=1.,gate_max_opening=.3,
        gate_sill_z=100.,Cd_gate_free=.6,Cd_gate_submerged=.8,neglect_approach_velocity_head=True)
    fig,axs=plt.subplots(1,2,figsize=(12.2,4.7),layout='constrained')
    r=np.linspace(.5,4,1000);H=.12
    colors=['#668eae','#1b637d','#80a7a1']
    for cu,color in zip((2.6,2.8,3.),colors):
        s=Assumptions(Cw_US=cu)
        q=[gate_flow(100+H,98,H/z/.3,p,s).discharge/H**1.5 for z in r]
        axs[0].plot(r,q,color=color,label=f'D candidate: Cw = {cu:g} US',lw=2 if cu==2.8 else 1.3)
        heads=np.linspace(.05000001,.45,900)
        q=[gate_flow(100+h,98,.1/.3,p,s).discharge for h in heads]
        axs[1].plot(heads,q,color=color,lw=2 if cu==2.8 else 1.3)
    use=r>=1.25
    axs[0].plot(r[use],.6*sqrt(2*p.gravity)/r[use],color='#b04c36',ls='--',label='Original B (supported domain)')
    axs[0].axhline(sqrt(p.gravity)*(2/3)**1.5,color='black',lw=1,ls=':',label='Ideal broad-crest energy bound')
    hs=np.linspace(.125,.45,500)
    axs[1].plot(hs,[sluice_gate_flow(100+h,98,.1/.3,p).discharge for h in hs],color='#b04c36',ls='--')
    axs[0].axvspan(1,1.5,alpha=.09,color='#d9a327')
    axs[1].axvspan(.1,.15,alpha=.09,color='#d9a327')
    axs[1].axvspan(0,.05,color='gray',alpha=.16)
    axs[0].set(xlabel='Upstream energy head / opening, H/a',ylabel='Q / (b H$^{3/2}$)  [m$^{1/2}$/s]',
        title='Fixed head H = 0.12 m; unchanged clear width')
    axs[1].set(xlabel='Upstream energy head H [m]',ylabel='Gate discharge [m³/s]',
        title='Fixed actual opening a = 0.10 m',xlim=(0,.45),ylim=(0,None))
    axs[1].text(.024,.17,'Outside\ncrest\ndomain',rotation=90,ha='center',va='center',fontsize=9)
    axs[0].legend(fontsize=8,loc='upper right')
    for ax in axs:ax.grid(alpha=.18)
    fig.suptitle('Stage D candidate rating: consistency checked, transition not empirically validated',fontsize=12)
    fig.savefig(out/'hydraulic_rating_comparison.png',dpi=180)
    plt.close(fig)


if __name__=='__main__':main()
