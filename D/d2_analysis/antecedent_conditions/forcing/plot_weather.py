"""Plot already prepared antecedent/event inputs; no ODE or data generation."""
from pathlib import Path
import hashlib
import json
import os
import sys
sys.dont_write_bytecode = True
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/cee690-antecedent-weather-mpl')
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
from D.g5_diagnostics import runner as historical
from D.d2_analysis.antecedent_conditions.forcing import provider


def main():
    cases = [('historical_20190825T0200Z', 'August 2019'), ('historical_20240804T1500Z', 'August 2024')]
    specs = [(2,3.6e6,'Precipitation (mm/h)','#2970a1'), (3,1.,'Air temperature (°C)','#b9504b'),
             (4,1.,'Solar irradiance (W/m²)','#be8521'), (5,1.,'Wind speed (m/s)','#26856e')]
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
    fig,axes=plt.subplots(4,2,figsize=(14,10.5),sharex=True,sharey='row',layout='constrained')
    records=[]
    row_maxima=np.zeros(4)
    for column,(event_id,label) in enumerate(cases):
        record=provider.record_for(event_id,72)
        frame,_=provider.read_forcing(record['forcing_file'])
        fixed=historical.load_fixed(dt=15.,design='best')
        fixed['model_config']['numerics']['duration_s']=record['total_duration_s']
        fixed['manifest']['duration_s']=record['total_duration_s']
        validation=provider.forcing_conservation(frame,fixed)
        assert validation['pass_conservation'] and validation['event_weather_and_boundary_inputs_bit_exact']
        segments=list(provider.hourly_input_segments(frame,fixed))
        # Keep duplicate hourly-boundary coordinates: left and right input values
        # are distinct for piecewise-constant rain/solar and the documented seam.
        time=np.concatenate([t for t,_ in segments])/3600-72
        values=np.concatenate([u for _,u in segments],axis=1)
        for row,(index,scale,ylabel,color) in enumerate(specs):
            axis=axes[row,column]
            row_maxima[row]=max(row_maxima[row],float(np.max(values[index]*scale)))
            axis.axvspan(-72,0,color='#dbe7ee',alpha=.65,zorder=0)
            axis.axvline(-24,color='#6e627e',linestyle='--',lw=1.1)
            axis.axvline(0,color='0.25',lw=1.2)
            axis.plot(time,values[index]*scale,color=color,lw=1.2)
            axis.set_ylabel(ylabel)
            axis.grid(alpha=.16)
            axis.set_xlim(-72,144)
            axis.set_xticks(np.arange(-72,145,24))
            if row==0:
                axis.set_title(f"{label} | event t = 0: {record['event_start_utc']}\nPrefix rain: {record['antecedent_rain_mm']:.1f} mm",loc='left',fontsize=11)
            if row in (0,2,3):axis.set_ylim(bottom=0)
        records.append({'event_id':event_id,'forcing_sha256':record['forcing_sha256'],
                        'original_event_forcing_sha256':record['original_event_forcing_sha256'],
                        'mapping_check':validation})
    for row in (0,2,3):
        axes[row,0].set_ylim(0,row_maxima[row]*1.06)
    for axis in axes[-1]:axis.set_xlabel('Time relative to original event start (h)')
    fig.suptitle('Saved weather inputs: 72-hour antecedent period + unchanged 144-hour event\nHourly rain/solar steps and the original temperature/wind mapping are retained',fontsize=14)
    fig.legend(handles=[Patch(facecolor='#dbe7ee',alpha=.65,label='72 h antecedent period'),
                        Line2D([0],[0],color='#6e627e',ls='--',label='−24 h: shorter conditioning start'),
                        Line2D([0],[0],color='0.25',label='0 h: original event start')],
               loc='outside lower center',ncols=3,frameon=False)
    output=HERE.parent/'plots/05_weather_inputs.png';output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output);plt.close(fig)
    metadata={'figure':str(output.relative_to(ROOT)),'source':records,'time_grid_s':15,
              'model_executed':False,'input_generation':'Existing provider mapping only; no new weather values or downloads',
              'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'weather_figure_provenance.json').write_text(json.dumps(metadata,indent=2)+'\n')
    (HERE/'weather_figure_caption.md').write_text('''# Weather-input figure caption

[Figure 05](../plots/05_weather_inputs.png). The four meteorological inputs for the August 2019 and August 2024 experiments are shown over the existing 72-hour antecedent period and the unchanged 144-hour event. Time zero marks the original event start; the shaded region is the 72-hour prefix and the dashed −24-hour line marks the shorter conditioning start. The plotted rain and solar steps and temperature/wind curves are the actual separately mapped inputs on the 15-second numerical grid, not simulated weather or interpolated model outputs. Duplicate coordinates at hourly boundaries retain both one-sided values.

Rainfall is the existing overlap allocation of native station accumulations; temperature and wind are regional station estimates; solar is a POWER grid estimate. The 2019 prefix contains zero quantified rain. The 2024 72-hour prefix contains 18.5 mm, while its last 24 hours contain zero quantified rain. Trace lower bounds and the original quality flags are retained in the source files. These are regional inputs transferred to the hypothetical installation, not observations of its storages.

The original event's endpoint holds are preserved by mapping it separately from the prefix. This produces an explicit temperature/wind representation jump at time zero: −1.5 °C/−1.65 m/s for 2019 and +1.27631579 °C/0 m/s for 2024. These jumps are not claims of measured instantaneous weather changes. At the beginning of a separate 24-hour run, the original mapper also applies its own first-half-hour temperature/wind hold; the dashed line identifies that start time but the displayed 72-hour curve is not asserted identical to the 24-hour mapping during that first half-hour. Both conditioning versions preserve the subsequent original event inputs exactly.

This figure executes no ODE and makes no claim of physical-state equilibration. Source hashes and exact event-mapping checks are in `weather_figure_provenance.json`. Reproduce using `/opt/miniconda3/bin/python3 -B D/d2_analysis/antecedent_conditions/forcing/plot_weather.py`.
''')
    print(output)

if __name__=='__main__':main()
