# Wetland Restoration and Flood Mitigation

**CEE690 — Systems Modeling for Climate Response, Duke University**

**HW A.2 — application of the ODE module to the wetland restoration and flood mitigation model**

This project implements an eight-state ODE model to explore water storage and routing through two dam–gate–wetland branches and a downstream river. The workflow is **model definition → algebraic relationships → ODE implementation → baseline simulation → interpretation**. The numerical example is **illustrative, synthetic, and uncalibrated**.

## System and inputs

The states are upland surface water, upland soil water, two wetlands, two dam pools, groundwater, and a river reach, in this order:

`x = [S_u, S_m, S_w1, S_w2, S_d1, S_d2, S_g, S_r]`

Each main branch follows **upland → dam pool → gate → wetland → river**. Direct upland runoff to the river and the soil–groundwater pathway remain included. The river outlet is the downstream boundary; external inflow and signed external groundwater exchange also cross the system boundary.

- **Controls:** two prescribed gate-opening histories, each between 0 and 1.
- **Fixed wetland settings:** vegetation category and fractional cover, selected independently for each branch and held fixed during the run.
- **Disturbances:** precipitation, temperature, solar radiation, wind, external river inflow, downstream water level, external groundwater level, and outlet blockage.

Storage is in m³, time in seconds, and flow in m³/s. Internal flows are one-way; receiving-water levels can suppress discharge.

## Main files

| File | Purpose |
|---|---|
| [model_notes.md](model_notes.md) | System boundary, eight ODEs, main algebraic relationships, figure guide, and baseline interpretation. |
| [wetland_flood_model.py](wetland_flood_model.py) | Numerical constants in `class c`, model RHS, synthetic inputs, official solver integration, and baseline outputs. |
| [algebraic_plots.py](algebraic_plots.py) | Generates A–G without loading the ODE solver. |
| [requirements.txt](requirements.txt) | NumPy and Matplotlib versions used locally. |
| [.gitignore](.gitignore) | Excludes local records, results, environments, caches, and private/reference material. |
| [plots/](plots/) | Existing A–J figures for this illustrative example. |

## Install and run

Run these commands from the project root on macOS/Linux with Python 3.13 available as `python3`. The existing environment used Python **3.13.5**, NumPy **2.2.6**, and Matplotlib **3.10.6**.

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

### Algebraic plots only

```sh
python algebraic_plots.py
```

This generates A–G using NumPy and Matplotlib; it does not require the ODE solver.

### Baseline ODE simulation

The baseline uses Henri P. Gavin's official [multivarious solver](https://github.com/hpgavin/multivarious), specifically `ode4u`. The verified source revision is `06d7aed5b541ec7fe26a6dc69be6aa7fb610aa02` (package version 1.2.1). Its package initializer also requires SciPy; SciPy **1.16.2** was used locally. Install the baseline dependencies separately (the VCS install requires Git):

```sh
python -m pip install scipy==1.16.2
python -m pip install --no-deps "multivarious @ git+https://github.com/hpgavin/multivarious.git@06d7aed5b541ec7fe26a6dc69be6aa7fb610aa02"
python wetland_flood_model.py
```

If the verified official checkout already exists at `work/multivarious`, the loader uses it before the installed package, so the VCS installation can be omitted. That checkout remains an external, ignored dependency. No solver is bundled or substituted.

The script generates A–G, runs one six-day baseline with a 120 s step, and generates H–J. Figures go to `plots/`; numerical arrays and diagnostics go to the ignored `results/baseline.npz` and `results/baseline_summary.json`. Rerunning overwrites these outputs.

**Execution record:** algebraic plotting and the baseline were run in the existing environment with the verified official checkout. A fresh virtual-environment installation and the VCS installation commands above have not been tested. This documentation cleanup did not rerun or change the simulation.

## Figures and interpretation

A–G show gate release, wetland drainage, wetland overflow, upland infiltration, evaporative demand, river-outlet blockage, and groundwater exchange, respectively. H shows all eight storage histories; I shows river and wetland responses; J shows the two controls and eight disturbances. [Model notes](model_notes.md) explain each figure.

![Illustrative baseline river and wetland responses](plots/I_baseline_response.png)

The recorded baseline gives peak river outflow **1.357896 m³/s** and peak river level **97.879407 m**. Storage and release timing shape the response; this single run does not establish a wetland-restoration benefit. Fixed geometry, empirical flow laws, simplified evaporation, and synthetic inputs limit interpretation. Numerical balance checks are not physical validation, and the results are not real flood predictions.

The implementation follows the course-provided `earth_heat_2.py` template (`class c`, `system(t,x,u,c)`, and `ode4u`). Professor handouts and local development records are excluded from the submitted files. No project license has been selected.
