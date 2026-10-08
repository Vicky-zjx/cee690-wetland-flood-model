#!/usr/bin/env python3
"""Stage A: accepted eight-state model; all site values are illustrative.

Run from any directory: python3 /path/to/Modeling/wetland_flood_model.py
The official multivarious checkout is kept in work/multivarious.
Workflow: inspect seven algebraic plots, run one baseline, inspect its outputs.
"""
from pathlib import Path
import inspect
import json
import os
import sys
import numpy as np

ROOT = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "work" / "mplconfig"))
LOCAL_PACKAGE = ROOT / "work" / "multivarious"
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

STATE_NAMES = ["S_u", "S_m", "S_w1", "S_w2", "S_d1", "S_d2", "S_g", "S_r"]
INPUT_NAMES = ["g1", "g2", "P", "T_a", "R_s", "v", "Q_ext", "h_d", "h_g", "B"]
U = {name: i for i, name in enumerate(INPUT_NAMES)}
OUTPUT_NAMES = [
    "H_w1", "H_w2", "H_d1", "H_d2", "H_g", "H_r", "Q_G1", "Q_G2",
    "Q_wr1", "Q_wr2", "Q_dsp1", "Q_dsp2", "Q_wsp1", "Q_wsp2", "Q_inf",
    "Q_mg", "Q_wg1", "Q_wg2", "Q_gr", "Q_gx", "Q_out", "E_tot", "S_tot",
    "E_star", "Q_gx_in", "Q_gx_out", "Q_ud1", "Q_ud2", "Q_ur", "E_m",
    "E_w1", "E_w2", "E_d1", "E_d2", "E_r", "P_total", "external_net",
]
Y = {name: i for i, name in enumerate(OUTPUT_NAMES)}


def load_ode4u():
    """Load the official solver only when a baseline simulation is requested."""
    if LOCAL_PACKAGE.is_dir() and str(LOCAL_PACKAGE) not in sys.path:
        sys.path.insert(0, str(LOCAL_PACKAGE))
    try:
        from multivarious.ode import ode4u
    except ImportError as error:
        raise ImportError(
            "The baseline requires the official multivarious ode4u and its dependencies. "
            "Provide the project-local checkout at work/multivarious or an importable "
            "installation. Algebraic plots need only NumPy and Matplotlib. "
            f"Import detail: {error}"
        ) from error
    return ode4u


class c:
    """Parameters, designs, forcing settings, and simulation settings (SI).

    Except for v3's physical constants, these are uncalibrated assumptions.
    Branch arrays always use [branch 1, branch 2].
    """
    def __init__(self):
        # Simulation and reporting; thresholds are not hard storage bounds.
        self.hour = 3600.0
        self.day = 24 * self.hour
        self.duration = 6 * self.day
        self.time_step = 120.0
        self.H_flood = 97.90                 # m; illustrative reporting threshold
        self.H_r_valid_max = 99.0            # m; assumed geometry-use limit
        self.H_g_valid_max = 97.9            # m; lowest wetland bed
        self.soil_bottom = 100.0             # m; separate upper soil region
        self.balance_rtol = 1e-4             # illustrative numerical criterion

        # Disjoint surface rainfall footprints; aquifer is NOT added to rain area.
        self.A_u = 600000.0                  # m2
        self.A_d = np.array([10000., 8000.])
        self.A_w = np.array([12000., 10000.])
        self.A_r = 20000.0
        self.A_g = self.A_u + sum(self.A_d) + sum(self.A_w)
        self.S_y = 0.20
        self.z_d = np.array([100.0, 99.8])    # m, common vertical datum
        self.z_w = np.array([98.0, 97.9])
        self.z_r = 97.0
        self.z_g = 94.0

        # Soil and upland routing (v3 Section 5.1).
        self.C_m = 120000.0                  # m3 = A_u * 0.20 m water capacity
        self.S_fc = 60000.0                  # m3 = A_u * 0.10 m
        self.K_u = 1e-6                      # m/s
        self.k_mg = 1 / (5 * self.day)        # 1/s
        self.k_ud = 1 / (np.array([8., 10.]) * self.hour)
        self.k_ur = 1 / (20 * self.hour)
        self.Q_udcap = np.array([1.2, 1.0])   # m3/s
        self.Q_urcap = 0.7

        # Dam pools and gates (v3 Sections 5.2-5.3).
        self.z_G = np.array([100.05, 99.85])  # m
        self.R_G = np.array([1.6, 1.8])       # s/m2; empirical linear rating
        self.Q_Gcap = np.array([1.0, 0.8])    # m3/s
        self.z_dsp = np.array([100.9, 100.7]) # m
        self.c_dsp = np.array([4.0, 3.2])     # m^(3/2)/s

        # Wetlands: independent categories and cover, fixed during each run.
        self.z_wo = np.array([98.05, 97.95])
        self.z_wsp = np.array([98.6, 98.5])
        self.c_wsp = np.array([3.0, 2.5])
        self.Q_wrcap = np.array([0.65, 0.55])
        self.R_bare = np.array([0.8, 1.0])    # s/m2
        self.R_veg = {"type_a": np.array([2.8, 3.0]),
                      "type_b": np.array([4.0, 4.5])}
        self.K_bare = np.array([1.0, 1.0])    # dimensionless ET coefficients
        self.K_veg = {"type_a": np.array([1.1, 1.1]),
                      "type_b": np.array([1.25, 1.25])}
        self.tau = ["type_a", "type_a"]       # synthetic labels, not species data
        self.rho = np.array([0.65, 0.60])

        # Groundwater and final river outlet.
        self.C_wg = np.array([0.015, 0.012])  # m2/s
        self.C_gr = 0.08                     # m2/s
        self.C_gx = 0.035                    # m2/s
        self.z_ro = 97.1                     # m
        self.R_r = 0.45                      # s/m2
        self.Q_rcap = 3.0                    # m3/s

        # Half-availability volumes V* (constitutive, not numerical corrections).
        self.V_u = 600.0                     # m3 = 1 mm * A_u
        self.V_m = 6000.0                    # m3 = 10 mm * A_u
        self.V_w = 0.01 * self.A_w           # m3 = 1 cm depth
        self.V_d = 0.01 * self.A_d
        self.V_g = 0.05 * self.S_y * self.A_g # m3 = 5 cm groundwater head
        self.V_r = 0.01 * self.A_r

        # Reduced Penman-type demand from v3 Section 6 (not FAO ET0).
        self.lambda_v = 2.45e6               # J/kg, v3 physical approximation
        self.c_p = 1013.0                    # J/(kg K)
        self.epsilon = 0.622
        self.rho_w = 1000.0                  # kg/m3
        self.R_v = 461.5                     # J/(kg K)
        self.p_0 = 101325.0                  # Pa; illustrative constant pressure
        self.r_H = 0.65                      # assumed fixed relative humidity
        self.eta = 0.60
        self.k_0 = 0.001                     # m/s
        self.k_v = 0.002                     # dimensionless
        self.K_ETm = 0.85
        self.K_ETd = np.array([1.0, 1.0])
        self.K_ETr = 1.0

        # Initial conditions derived from assumed depths/wetness.
        self.u_depth_init = 0.001            # m
        self.soil_fraction_init = 0.75
        self.w_depth_init = np.array([0.25, 0.25])
        self.d_depth_init = np.array([0.30, 0.30])
        self.r_depth_init = 0.35
        self.H_g_init = 97.20                # m

        # Synthetic forcing knots in hours; each triple is start/peak/end.
        self.rain_pulses = [(18., 22., 26., 15.), (66., 71., 76., 20.)] # mm/h
        self.ext_pulses = [(20., 25., 36., 0.35), (68., 76., 90., 0.70)]
        self.Q_ext_base = 0.08               # m3/s; outside rainfall footprint
        self.weather_period_h = 24.0        # h
        self.temperature_phase_h = 9.0      # h; sine peaks at 15 h
        self.sunrise_phase_h = 6.0          # h
        self.wind_period_h = 36.0           # h
        self.groundwater_period_h = 144.0    # h
        self.T_mean = 22.0                   # degC
        self.T_amplitude = 4.0
        self.R_s_peak = 550.0                # W/m2
        self.wind_mean = 2.0                 # m/s
        self.wind_amplitude = 0.6
        self.h_d_base = 97.15                # m
        self.h_d_pulse = (66., 77., 94., 0.10)
        self.h_g_base = 97.15                # m
        self.h_g_amplitude = 0.025           # m; slow prescribed boundary signal
        self.B_base = 0.05
        self.gate1_hours = [0, 17, 18, 26, 28, 50, 65, 66, 76, 80, 110, 144]
        self.gate1_values = [.35,.35,.20,.20,.50,.35,.35,.20,.20,.50,.35,.35]
        self.gate2_hours = [0, 18, 20, 27, 30, 52, 66, 68, 78, 82, 112, 144]
        self.gate2_values = [.45,.45,.30,.30,.60,.45,.45,.30,.30,.60,.45,.45]


    def refresh_geometry_defaults(self):
        """Explicitly reapply five linked defaults; overwrites custom values of them.

        Call after geometry edits when these depth/footprint assumptions are wanted.
        Independent soil, hydraulic, ET, and other empirical parameters stay unchanged.
        """
        self.A_g = self.A_u + sum(self.A_d) + sum(self.A_w)
        self.V_w = 0.01 * self.A_w
        self.V_d = 0.01 * self.A_d
        self.V_g = 0.05 * self.S_y * self.A_g
        self.V_r = 0.01 * self.A_r


    @property
    def time(self):
        n = int(round(self.duration / self.time_step))
        if not np.isclose(n * self.time_step, self.duration):
            raise ValueError("Duration must be a whole number of timesteps.")
        return np.arange(n + 1, dtype=float) * self.time_step

    @property
    def A_rain(self):
        return self.A_u + sum(self.A_d) + sum(self.A_w) + self.A_r

    @property
    def C_d(self):
        return self.A_d * (self.z_dsp - self.z_d)

    @property
    def C_w(self):
        return self.A_w * (self.z_wsp - self.z_w)

    @property
    def x_init(self):
        return np.array([self.A_u * self.u_depth_init,
                         self.C_m * self.soil_fraction_init,
                         *(self.A_w * self.w_depth_init),
                         *(self.A_d * self.d_depth_init),
                         (self.H_g_init-self.z_g)*self.S_y*self.A_g,
                         self.A_r*self.r_depth_init])

    def design(self):
        # v3 Section 7: categorical endpoints, independently selected per branch.
        Rv = np.array([self.R_veg[self.tau[i]][i] for i in range(2)])
        Kv = np.array([self.K_veg[self.tau[i]][i] for i in range(2)])
        return ((1-self.rho)*self.R_bare + self.rho*Rv,
                (1-self.rho)*self.K_bare + self.rho*Kv)

def positive_part(z):
    return np.maximum(z, 0.0)

def availability(S, V):
    return S / (S + V)

def stage(S, A, z):
    return z + S / A

def head_flow(Ha, Hb, zo, R, Qcap):
    return np.minimum(Qcap, positive_part(Ha-np.maximum(Hb, zo))/R)

def evaporative_demand(Ta, Rs, v, c):
    # Exact reduced demand in v3 Section 6.1, all pressures in Pa.
    TK = Ta + 273.15
    es = 610.8 * np.exp(17.27*Ta / (Ta+237.3))
    Delta = 17.27*237.3*es / (Ta+237.3)**2
    D_rho = (1-c.r_H)*es / (c.R_v*TK)
    k_a = c.k_0 + c.k_v*v
    gamma = c.c_p*c.p_0 / (c.epsilon*c.lambda_v)
    return (Delta/(Delta+gamma) * c.eta*Rs/(c.rho_w*c.lambda_v)
            + gamma/(Delta+gamma) * k_a*D_rho/c.rho_w)

def system(t, x, u, c):
    """v3 Eqs. (1)-(8). ode4u supplies one 10-element input column.

    Reject an inadmissible RK stage rather than hide it by altering states.
    Controls/disturbances are interpolated by ode4u, not looked up again here.
    """
    if not np.all(np.isfinite(x)) or np.any(x < 0) or x[1] > c.C_m:
        raise FloatingPointError(f"Inadmissible state at t={t:g} s: {x}. Reduce timestep.")
    S_u, S_m, S_w1, S_w2, S_d1, S_d2, S_g, S_r = x
    g1, g2, P, T_a, R_s, v, Q_ext, h_d, h_g, B = u
    Sw, Sd = np.array([S_w1,S_w2]), np.array([S_d1,S_d2])
    Hw, Hd = stage(Sw,c.A_w,c.z_w), stage(Sd,c.A_d,c.z_d)
    Hg = stage(S_g,c.S_y*c.A_g,c.z_g)
    Hr = stage(S_r,c.A_r,c.z_r)
    Rw, KETw = c.design()

    # v3 Section 5.1: each upland link is counted once at each end.
    Qud = np.minimum(c.Q_udcap,c.k_ud*S_u)
    Qur = min(c.Q_urcap,c.k_ur*S_u)
    Qinf = c.A_u*c.K_u*availability(S_u,c.V_u)*positive_part(1-S_m/c.C_m)
    Qmg = c.k_mg*positive_part(S_m-c.S_fc)

    # v3 Sections 5.2-5.3: pool -> gate -> wetland -> river.
    QG = np.array([g1,g2])*head_flow(Hd,Hw,c.z_G,c.R_G,c.Q_Gcap)
    Qwr = head_flow(Hw,Hr,c.z_wo,Rw,c.Q_wrcap)
    Qout = (1-B)*head_flow(Hr,h_d,c.z_ro,c.R_r,c.Q_rcap)
    Qdsp = c.c_dsp*positive_part(Hd-np.maximum(c.z_dsp,Hw))**1.5
    Qwsp = c.c_wsp*positive_part(Hw-np.maximum(c.z_wsp,Hr))**1.5

    # v3 Section 5.4: only the external groundwater exchange is signed.
    Qwg = c.C_wg*availability(Sw,c.V_w)*positive_part(Hw-Hg)
    Qgr = c.C_gr*availability(S_g,c.V_g)*positive_part(Hg-Hr)
    Qgx_in = c.C_gx*positive_part(h_g-Hg)
    Qgx_out = c.C_gx*availability(S_g,c.V_g)*positive_part(Hg-h_g)
    Qgx = Qgx_in-Qgx_out

    # v3 Section 6.2: soil footprint is not a second rain-receiving area.
    E_star = evaporative_demand(T_a,R_s,v,c)
    Em = c.A_u*c.K_ETm*availability(S_m,c.V_m)*E_star
    Ew = c.A_w*KETw*availability(Sw,c.V_w)*E_star
    Ed = c.A_d*c.K_ETd*availability(Sd,c.V_d)*E_star
    Er = c.A_r*c.K_ETr*availability(S_r,c.V_r)*E_star
    Etot = Em + sum(Ew) + sum(Ed) + Er

    dSu = P*c.A_u-Qinf-Qud[0]-Qud[1]-Qur                  # Eq. (1)
    dSm = Qinf-Qmg-Em                                    # Eq. (2)
    dSw1 = P*c.A_w[0]+QG[0]+Qdsp[0]-Qwr[0]-Qwsp[0]-Qwg[0]-Ew[0] # (3)
    dSw2 = P*c.A_w[1]+QG[1]+Qdsp[1]-Qwr[1]-Qwsp[1]-Qwg[1]-Ew[1] # (4)
    dSd1 = P*c.A_d[0]+Qud[0]-QG[0]-Qdsp[0]-Ed[0]          # Eq. (5)
    dSd2 = P*c.A_d[1]+Qud[1]-QG[1]-Qdsp[1]-Ed[1]          # Eq. (6)
    dSg = Qmg+Qwg[0]+Qwg[1]+Qgx-Qgr                       # Eq. (7)
    dSr = P*c.A_r+Qur+sum(Qwr)+sum(Qwsp)+Qgr+Q_ext-Qout-Er # Eq. (8)
    dxdt = np.array([dSu,dSm,dSw1,dSw2,dSd1,dSd2,dSg,dSr])
    external_net = P*c.A_rain+Q_ext+Qgx-Qout-Etot          # v3 Section 9.3
    y = np.array([*Hw,*Hd,Hg,Hr,*QG,*Qwr,*Qdsp,*Qwsp,Qinf,Qmg,
                  *Qwg,Qgr,Qgx,Qout,Etot,sum(x),E_star,Qgx_in,Qgx_out,
                  *Qud,Qur,Em,*Ew,*Ed,Er,P*c.A_rain,external_net])
    return dxdt, y

def triangle(hours, start, peak, end, amplitude):
    """Continuous piecewise-linear pulse, sampled in SI time histories."""
    return amplitude*positive_part(np.minimum((hours-start)/(peak-start),
                                              (end-hours)/(end-peak)))

def make_inputs(c):
    hours = c.time / c.hour
    rain_mmh = sum((triangle(hours,*pulse) for pulse in c.rain_pulses), np.zeros_like(hours))
    P = rain_mmh * 1e-3 / 3600.0               # explicit mm/hour -> m/second
    Ta = c.T_mean + c.T_amplitude*np.sin(2*np.pi*(hours-c.temperature_phase_h)/c.weather_period_h)
    Rs = c.R_s_peak*positive_part(np.sin(2*np.pi*(hours-c.sunrise_phase_h)/c.weather_period_h))
    v = c.wind_mean+c.wind_amplitude*np.sin(2*np.pi*hours/c.wind_period_h)
    Qext = c.Q_ext_base+sum((triangle(hours,*p) for p in c.ext_pulses), np.zeros_like(hours))
    hd = c.h_d_base+triangle(hours,*c.h_d_pulse)
    hg = c.h_g_base+c.h_g_amplitude*np.sin(2*np.pi*hours/c.groundwater_period_h)
    B = np.full_like(hours,c.B_base)
    g1 = np.interp(hours,c.gate1_hours,c.gate1_values)
    g2 = np.interp(hours,c.gate2_hours,c.gate2_values)
    return np.vstack([g1,g2,P,Ta,Rs,v,Qext,hd,hg,B])


def check_groundwater_levels(levels, params, context):
    """Check this configuration's storage-separation limits without modifying water."""
    if not np.all(np.isfinite(params.z_w)):
        raise ValueError("Groundwater geometry limit z_w must be finite.")
    if not np.all(np.isfinite(levels)):
        raise ValueError(f"{context}: groundwater levels must be finite.")
    limits = [("H_g_valid_max", params.H_g_valid_max),
              ("soil_bottom", params.soil_bottom), ("min(z_w)", np.min(params.z_w))]
    for name, limit in limits:
        if not np.isfinite(limit):
            raise ValueError(f"Groundwater geometry limit {name} must be finite.")
        if np.any(np.asarray(levels) >= limit):
            raise ValueError(
                f"{context}: maximum H_g={np.max(levels):g} m violates "
                f"H_g < {name}={limit:g} m. Review the geometry/storage interpretation."
            )


def validate(c, u):
    """Check domains and full input length before ode4u can zero-pad anything."""
    if u.shape != (10,len(c.time)) or not np.all(np.isfinite(u)):
        raise ValueError("Expected finite u with shape (10, len(c.time)).")
    if c.time_step <= 0 or c.duration <= 0 or not np.all(np.diff(c.time)>0):
        raise ValueError("Time must be strictly increasing.")
    for name in ["g1","g2","B"]:
        if np.any((u[U[name]]<0)|(u[U[name]]>1)):
            raise ValueError(f"{name} must lie in [0,1].")
    for name in ["P","R_s","v","Q_ext"]:
        if np.any(u[U[name]]<0):
            raise ValueError(f"{name} must be nonnegative.")
    if np.any(u[U["T_a"]] <= 0):
        raise ValueError("This liquid-water demonstration requires T_a > 0 degC.")
    positive = ["A_u","A_d","A_w","A_r","A_g","S_y","C_m","V_u","V_m",
                "V_w","V_d","V_g","V_r","R_G","Q_Gcap","Q_udcap","Q_urcap",
                "k_ud","k_ur","k_mg","R_bare","Q_wrcap","R_r","Q_rcap",
                "c_dsp","c_wsp","lambda_v","c_p","epsilon","rho_w","R_v","p_0"]
    nonnegative = ["K_u","C_wg","C_gr","C_gx","K_bare","K_ETm","K_ETd",
                   "K_ETr","k_0","k_v"]
    for name in positive+nonnegative:
        a = np.asarray(getattr(c,name))
        if not np.all(np.isfinite(a)) or np.any(a <= 0 if name in positive else a < 0):
            raise ValueError(f"Invalid parameter: {name}")
    for name in ["z_d","z_w","z_r","z_g","z_G","z_wo","z_dsp","z_wsp","z_ro",
                 "H_g_valid_max","soil_bottom","H_r_valid_max","H_g_init"]:
        if not np.all(np.isfinite(getattr(c,name))):
            raise ValueError(f"Nonfinite elevation: {name}")
    if not (0 <= c.S_fc <= c.C_m and 0 < c.S_y <= 1 and
            0 <= c.r_H <= 1 and 0 <= c.eta <= 1 and
            np.all((c.rho >= 0)&(c.rho <= 1))):
        raise ValueError("Invalid soil, yield, humidity, radiation, or design domain.")
    for category in set(c.tau):
        if category not in c.R_veg or category not in c.K_veg:
            raise ValueError("Unknown categorical vegetation type.")
        if (not np.all(np.isfinite(c.R_veg[category])) or
            not np.all(np.isfinite(c.K_veg[category])) or
            np.any(c.R_veg[category]<=0) or np.any(c.K_veg[category]<0)):
            raise ValueError("Invalid design endpoints.")
    # Ordering assumptions of this project configuration, not universal laws.
    for i in range(2):
        if c.z_dsp[i] <= c.z_G[i]:
            raise ValueError(f"Branch {i+1}: z_dsp[{i}] must be > z_G[{i}].")
        if c.z_wsp[i] <= c.z_wo[i]:
            raise ValueError(f"Branch {i+1}: z_wsp[{i}] must be > z_wo[{i}].")
    if (np.any(c.z_G<c.z_d) or np.any(c.z_wo<c.z_w) or c.z_ro<c.z_r or
        np.any(c.z_dsp<=c.z_d) or np.any(c.z_wsp<=c.z_w)):
        raise ValueError("Outlet/spill threshold lies below required donor bed.")
    check_groundwater_levels(c.H_g_init, c, "Initial condition")
    system(0,c.x_init,u[:,0],c)

def cumulative_trapezoid(values, time):
    return np.r_[0.,np.cumsum(0.5*(values[:-1]+values[1:])*np.diff(time))]

def flood_duration(time, H, threshold):
    """Time above threshold with linear crossings between sampled stages."""
    a, b = H[:-1]-threshold, H[1:]-threshold
    fraction = ((a>0)&(b>0)).astype(float)
    up = (a<=0)&(b>0)
    down = (a>0)&(b<=0)
    fraction[up] = b[up]/(b[up]-a[up])
    fraction[down] = a[down]/(a[down]-b[down])
    return float(np.sum(fraction*np.diff(time)))


def baseline_diagnostics(params, time, states, derivatives, inputs, outputs, solver):
    """A compact execution check; this is not physical model validation."""
    if not all(np.all(np.isfinite(a)) for a in [states, derivatives, outputs]):
        raise FloatingPointError("The solver returned incomplete or nonfinite results.")
    if np.any(states < 0) or np.any(states[1] > params.C_m):
        raise FloatingPointError("Storage bounds were violated; inspect the timestep.")
    check_groundwater_levels(outputs[Y["H_g"]], params, "Computed trajectory")
    if not np.isfinite(params.H_r_valid_max):
        raise ValueError("River geometry limit H_r_valid_max must be finite.")
    if np.any(outputs[Y["H_r"]] >= params.H_r_valid_max):
        raise ValueError("Computed river level violates H_r < H_r_valid_max; review geometry.")
    net_flux = outputs[Y["external_net"]]
    balance = (outputs[Y["S_tot"]] - outputs[Y["S_tot"], 0]
               - cumulative_trapezoid(net_flux, time))
    event_inflow = np.trapezoid(outputs[Y["P_total"]] + inputs[U["Q_ext"]]
                               + outputs[Y["Q_gx_in"]], time)
    scale = event_inflow if event_inflow > 0 else outputs[Y["S_tot"], 0]
    absolute_error = float(np.max(np.abs(balance)))
    relative_error = absolute_error / scale if scale > 0 else (0. if absolute_error == 0 else float("inf"))
    summary = {
        "example": "Illustrative, synthetic, uncalibrated six-day baseline",
        "solver_file": inspect.getfile(solver), "dt_s": params.time_step,
        "Q_peak_m3s": float(np.max(outputs[Y["Q_out"]])),
        "H_peak_m": float(np.max(outputs[Y["H_r"]])),
        "V_out_m3": float(np.trapezoid(outputs[Y["Q_out"]], time)),
        "D_flood_h": flood_duration(time, outputs[Y["H_r"]], params.H_flood) / params.hour,
        "balance_max_abs_m3": absolute_error,
        "balance_final_m3": float(balance[-1]),
        "balance_max_relative": float(relative_error),
        "event_inflow_m3": float(event_inflow),
        "derivative_balance_max_m3s": float(np.max(np.abs(derivatives.sum(axis=0) - net_flux))),
        "minimum_storage_m3": {name: float(np.min(states[i])) for i, name in enumerate(STATE_NAMES)},
        "max_soil_fraction": float(np.max(states[1]) / params.C_m),
        "max_groundwater_level_m": float(np.max(outputs[Y["H_g"]])),
    }
    if relative_error >= params.balance_rtol:
        raise ValueError("Water-balance error exceeds the illustrative numerical tolerance.")
    return summary, balance


def save_baseline_figure(fig, filename, title):
    fig.suptitle(title + " — illustrative, uncalibrated", fontsize=13)
    for ax in fig.axes:
        ax.grid(alpha=0.25)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(ROOT / "plots" / filename, dpi=160)
    plt.close(fig)


def plot_baseline(params, time, states, inputs, outputs):
    """Three figures: states, selected responses, and the ten input channels."""
    hours = time / params.hour
    fig, axes = plt.subplots(4, 2, figsize=(11, 10), sharex=True)
    for i, ax in enumerate(axes.flat):
        ax.plot(hours, states[i], label=STATE_NAMES[i])
        ax.set(title=STATE_NAMES[i], ylabel="Storage (m³)", xlabel="Elapsed time (h)")
    save_baseline_figure(fig, "H_baseline_storages.png", "Eight storage states")

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    groups = [("H_r",), ("Q_out",), ("Q_G1", "Q_G2"),
              ("Q_wr1", "Q_wr2", "Q_wsp1", "Q_wsp2")]
    titles = ["River level", "Final river outflow", "Gate releases", "Wetland drainage and overflow"]
    for i, (ax, names) in enumerate(zip(axes.flat, groups)):
        for name in names:
            ax.plot(hours, outputs[Y[name]], label=name)
        ax.set(title=titles[i], xlabel="Elapsed time (h)",
               ylabel="Level (m, common datum)" if i == 0 else "Flow (m³/s)")
        if i == 0:
            ax.axhline(params.H_flood, color="gray", ls="--", label="Illustrative flood threshold")
        ax.legend(fontsize=8)
    save_baseline_figure(fig, "I_baseline_response.png", "Interpreting the baseline response")

    fig, axes = plt.subplots(5, 2, figsize=(11, 13), sharex=True)
    units = ["Opening (0–1)", "Opening (0–1)", "Rain (mm/h)", "Temperature (°C)",
             "Radiation (W/m²)", "Wind (m/s)", "Inflow (m³/s)", "Level (m)", "Level (m)", "Blockage (0–1)"]
    for i, ax in enumerate(axes.flat):
        values = inputs[i] * 3.6e6 if INPUT_NAMES[i] == "P" else inputs[i]
        ax.plot(hours, values)
        ax.set(title=INPUT_NAMES[i], ylabel=units[i], xlabel="Elapsed time (h)")
    save_baseline_figure(fig, "J_baseline_inputs.png", "Two controls and eight disturbances")


if __name__ == "__main__":
    # 1. Define model parameters, fixed designs, time, and initial conditions.
    params = c()
    (ROOT / "plots").mkdir(exist_ok=True)
    (ROOT / "results").mkdir(exist_ok=True)
    x_init = params.x_init

    # 2. Inspect the algebraic relationships before the whole-system simulation.
    from A.algebraic_plots import make_algebraic_plots
    make_algebraic_plots(params)

    # 3. Build controls/disturbances and integrate with the professor's solver.
    u = make_inputs(params)
    validate(params, u)
    try:
        ode4u = load_ode4u()
    except ImportError as error:
        raise SystemExit(str(error)) from error
    time, x, x_drv, y = ode4u(system, params.time, x_init, u, params)

    # 4. Check the baseline, save a small result set, and inspect state/output plots.
    summary, balance = baseline_diagnostics(params, time, x, x_drv, u, y, ode4u)
    np.savez_compressed(ROOT / "results" / "baseline.npz", time=time, x=x,
                        x_drv=x_drv, u=u, y=y, balance=balance,
                        state_names=STATE_NAMES, input_names=INPUT_NAMES, output_names=OUTPUT_NAMES)
    (ROOT / "results" / "baseline_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    plot_baseline(params, time, x, u, y)
    print("Stage A numerical workflow completed; see the remaining-work checklist.")
    print(f"Q_peak = {summary['Q_peak_m3s']:.6f} m³/s; H_peak = {summary['H_peak_m']:.6f} m")
    print(f"Maximum water-balance residual = {summary['balance_max_abs_m3']:.6f} m³")
    print("All site values are illustrative and uncalibrated. See model_notes.md.")
