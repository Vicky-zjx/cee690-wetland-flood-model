#!/usr/bin/env python3
"""Seven educational views of the accepted v3 algebraic relationships.

Run: python3 -m A.algebraic_plots
These are illustrative parameter comparisons, not calibrated site ranges.
The model file sets the Matplotlib cache/backend before pyplot is imported.
"""
from A.wetland_flood_model import (
    c, positive_part, availability, head_flow, evaporative_demand, ROOT,
)
import numpy as np
import matplotlib.pyplot as plt


def make_algebraic_plots(params=None):
    """Save A-G in plots/ and return their paths; no ODE simulation is run."""
    params = c() if params is None else params
    plot_dir = ROOT / "plots"
    plot_dir.mkdir(exist_ok=True)
    paths = []

    def finish(fig, filename, note):
        for ax in fig.axes:
            ax.grid(alpha=0.25)
            ax.legend(fontsize=9)
        fig.text(0.5, 0.02, note + "\nIllustrative comparisons; not calibrated ranges.",
                 ha="center", va="bottom", fontsize=8)
        fig.tight_layout(rect=(0, 0.11, 1, 1))
        path = plot_dir / filename
        fig.savefig(path, dpi=160)
        plt.close(fig)
        paths.append(path)

    # A. Q_G = g * F: change gate opening at fixed resistance and rated capacity.
    # Zero/negative effective head gives no release; opening scales slope AND cap.
    # This matches one-way routing: out of dam pool (5/6), into wetland (3/4).
    head = np.linspace(-0.2, 3.0, 401)
    resistance, capacity = params.R_G[0], params.Q_Gcap[0]
    fig, ax = plt.subplots(figsize=(8, 5))
    for opening in (0.25, 0.50, 0.75, 1.00):
        discharge = opening * head_flow(head, 0, 0, resistance, capacity)
        ax.plot(head, discharge,
                label=f"g={opening:.2f}; cap={opening*capacity:.2f} m³/s")
    ax.axvline(0, color="0.5", linestyle=":", linewidth=1)
    ax.axvline(resistance*capacity, color="0.5", linestyle="--", linewidth=1)
    ax.set(xlabel=r"Effective gate head $H_d-\max(H_w,z_G)$ (m)",
           ylabel=r"Gate release $Q_G$ (m³/s)",
           title="A. Gate release: opening changes slope and effective capacity")
    finish(fig, "A_gate_release.png",
           f"Branch 1: R_G={resistance:g} s/m²; rated capacity={capacity:g} m³/s. "
           f"Dotted: zero head; dashed: cap begins at {resistance*capacity:g} m.")

    # B. Q_wr = F: compare bare and two vegetation designs for Wetland 1.
    # Higher effective resistance lowers the slope and delays the same flow cap.
    # This is a design assumption, not species evidence; outflow (3), inflow (8).
    head = np.linspace(-0.2, 3.0, 401)
    cover = 0.65
    designs = [("bare, cover=0", params.R_bare[0])]
    for category in ("type_a", "type_b"):
        resistance = ((1-cover)*params.R_bare[0]
                      + cover*params.R_veg[category][0])
        designs.append((f"{category}, cover={cover:g}", resistance))
    fig, ax = plt.subplots(figsize=(8, 5))
    for name, resistance in designs:
        discharge = head_flow(head, 0, 0, resistance, params.Q_wrcap[0])
        ax.plot(head, discharge, label=f"{name}; R={resistance:.2f} s/m²")
    ax.axvline(0, color="0.5", linestyle=":", linewidth=1)
    ax.set(xlabel=r"Effective wetland head $H_w-\max(H_r,z_{w,o})$ (m)",
           ylabel=r"Normal drainage $Q_{wr}$ (m³/s)",
           title="B. Wetland drainage: resistance delays the shared discharge cap")
    finish(fig, "B_wetland_drainage.png",
           "At the same effective head, higher resistance reduces the uncapped discharge "
           "and delays reaching the shared discharge cap.\n"
           f"Wetland 1; all curves retain the same rated cap "
           f"{params.Q_wrcap[0]:g} m³/s. Vegetation categories are synthetic labels.")

    # C. Q_wsp = c_wsp * [effective overflow head]_+^(3/2).
    # Vary the overflow coefficient; zero below onset and convex growth above it
    # are consistent with the accepted spill rating. Leaves (3), enters (8).
    head = np.linspace(-0.2, 1.0, 401)
    fig, ax = plt.subplots(figsize=(8, 5))
    for multiplier in (0.5, 1.0, 1.5):
        coefficient = multiplier*params.c_wsp[0]
        discharge = coefficient*positive_part(head)**1.5
        ax.plot(head, discharge, label=f"c_wsp={coefficient:g} m^(3/2)/s")
    ax.axvline(0, color="0.5", linestyle=":", linewidth=1)
    ax.set(xlabel=r"Overflow head $H_w-\max(z_{wsp},H_r)$ (m)",
           ylabel=r"Wetland overflow $Q_{wsp}$ (m³/s)",
           title="C. Wetland overflow: onset at zero, then 3/2-power growth")
    finish(fig, "C_wetland_overflow.png",
           "Zero marks the effective threshold, including river tailwater. "
           "This approximate overflow rating has no imposed flow cap.")

    # D. Q_inf = A_u*K_u*phi_u(S_u)*[1-S_m/C_m]_+; vary infiltration rate K_u.
    # More surface water supports infiltration, but full soil saturation stops it.
    # Both supply limits make physical sense; water leaves (1) and enters (2).
    saturation = np.linspace(0, 1, 301)
    surface_storage = np.linspace(0, 10*params.V_u, 301)
    fixed_storage, fixed_saturation = params.V_u, 0.50
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for multiplier in (0.5, 1.0, 2.0):
        rate = multiplier*params.K_u
        label = f"K_u={rate:.1e} m/s"
        flow = (params.A_u*rate*availability(fixed_storage, params.V_u)
                * positive_part(1-saturation))
        axes[0].plot(saturation, flow, label=label)
        flow = (params.A_u*rate*availability(surface_storage, params.V_u)
                * (1-fixed_saturation))
        axes[1].plot(surface_storage, flow, label=label)
    axes[0].axvline(1, color="0.5", linestyle=":", linewidth=1)
    axes[0].set(xlabel=r"Soil saturation $S_m/C_m$ (dimensionless)",
                ylabel=r"Infiltration $Q_{inf}$ (m³/s)",
                title=f"D1. Fixed upland storage S_u={fixed_storage:g} m³")
    axes[1].set(xlabel=r"Upland surface storage $S_u$ (m³)",
                ylabel=r"Infiltration $Q_{inf}$ (m³/s)",
                title=f"D2. Fixed soil saturation={fixed_saturation:g}")
    finish(fig, "D_upland_infiltration.png",
           f"A_u={params.A_u:g} m²; V_u*={params.V_u:g} m³; "
           f"soil capacity C_m={params.C_m:g} m³. Zero supply or full saturation stops infiltration.")

    # E. Accepted reduced demand E_* versus temperature; vary solar radiation.
    # Radiation adds energy; temperature changes vapor deficit and Penman weights.
    # This is a simplified demand, not standard ET0 or actual water withdrawal;
    # availability/area/coefficients turn it into losses in (2)-(6) and (8).
    temperature = np.linspace(5, 35, 301)
    wind = 2.0
    fig, ax = plt.subplots(figsize=(8, 5))
    for radiation in (100, 300, 600):
        demand = evaporative_demand(temperature, radiation, wind, params)
        ax.plot(temperature, demand*1000*3600,
                label=f"R_s={radiation:g} W/m²")
    ax.set(xlabel=r"Air temperature $T_a$ (°C)",
           ylabel=r"Reduced evaporative demand $E_*$ (mm/h)",
           title="E. Evaporative demand under three radiation conditions")
    finish(fig, "E_evaporative_demand.png",
           f"Wind={wind:g} m/s; assumed relative humidity={params.r_H:g}; "
           f"energy fraction η={params.eta:g}. Only display units change from m/s to mm/h.")

    # F. Q_out = (1-B)*F: vary outlet blockage at fixed resistance and capacity.
    # Blockage scales both the slope and cap; no discharge for adverse head.
    # This affects only the final outlet, a boundary loss from river ODE (8).
    head = np.linspace(-0.2, 2.0, 401)
    fig, ax = plt.subplots(figsize=(8, 5))
    for blockage in (0, 0.25, 0.50, 0.75):
        discharge = (1-blockage)*head_flow(head, 0, 0, params.R_r, params.Q_rcap)
        ax.plot(head, discharge,
                label=f"B={blockage:.2f}; cap={(1-blockage)*params.Q_rcap:.2f} m³/s")
    ax.axvline(0, color="0.5", linestyle=":", linewidth=1)
    ax.axvline(params.R_r*params.Q_rcap, color="0.5", linestyle="--", linewidth=1)
    ax.set(xlabel=r"Effective river head $H_r-\max(h_d,z_{r,o})$ (m)",
           ylabel=r"Final outlet discharge $Q_{out}$ (m³/s)",
           title="F. River outlet: blockage lowers effective capacity")
    finish(fig, "F_river_outlet_blockage.png",
           f"R_r={params.R_r:g} s/m²; unblocked cap={params.Q_rcap:g} m³/s. "
           "Dotted: zero head; dashed: capacity onset. B=1 would give zero everywhere.")

    # G. Q_gx = C_gx*[h_g-H_g]_+ - C_gx*phi_g(S_g)*[H_g-h_g]_+.
    # Vary conductance at fixed S_g: a higher outside level imports water; a lower
    # one exports available groundwater. The signed boundary term enters (7).
    head_difference = np.linspace(-1, 1, 401)
    groundwater_storage = params.x_init[6]
    available = availability(groundwater_storage, params.V_g)
    fig, ax = plt.subplots(figsize=(8, 5))
    for multiplier in (0.5, 1.0, 2.0):
        conductance = multiplier*params.C_gx
        exchange = conductance*(positive_part(head_difference)
                               - available*positive_part(-head_difference))
        ax.plot(head_difference, exchange, label=f"C_gx={conductance:g} m²/s")
    ax.axhline(0, color="0.5", linewidth=0.8)
    ax.axvline(0, color="0.5", linestyle=":", linewidth=1)
    ax.text(0.98, 0.08, "Negative: outward\nPositive: inward", transform=ax.transAxes,
            ha="right", fontsize=9)
    ax.set(xlabel=r"External minus internal groundwater level $h_g-H_g$ (m)",
           ylabel=r"Signed groundwater exchange $Q_{gx}$ (m³/s)",
           title="G. Groundwater boundary: direction reverses at equal levels")
    finish(fig, "G_groundwater_exchange.png",
           f"Fixed S_g={groundwater_storage:g} m³; outward availability φ_g={available:.4f}. "
           "Only outward exchange is reduced by the availability factor.")
    return paths


if __name__ == "__main__":
    for path in make_algebraic_plots():
        print(path)
