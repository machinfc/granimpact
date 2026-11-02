"""Result plots (matplotlib optional).

Lazy import on purpose: reading and summarizing the CSV works without matplotlib,
so that the simulator can be used on a server with no graphical environment.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Sequence

if TYPE_CHECKING:  # pragma: no cover
    from .results import Case


def _pyplot():
    try:
        import matplotlib

        matplotlib.use("Agg")  # headless: save to PNG
        import matplotlib.pyplot as plt

        return plt
    except ImportError as exc:  # pragma: no cover
        raise SystemExit(
            "matplotlib is not installed; tables and summaries still work.\n"
            "Install it with: pip install matplotlib"
        ) from exc


def plot_case(case: "Case", out: str | Path) -> Path:
    """4-panel figure: morphometry, energy, contacts and time fit."""
    plt = _pyplot()
    out = Path(out)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    fig.suptitle(f"GranImpact · {case.name}")

    # 1) Crater morphometry
    ax = axes[0][0]
    for column, label in (("D", "D [mm]"), ("d_exc", "d_exc [mm]"), ("h_rim", "h_rim [mm]")):
        t = [x for x, y in zip(case.data["t"], case.data[column]) if y is not None]
        y = [v * 1000.0 for v in case.finite(column)]
        if y:
            ax.plot(t, y, label=label)
    ax.set_xlabel("t [s]")
    ax.set_ylabel("length [mm]")
    ax.set_title("Morphometry")
    ax.legend(fontsize=8)

    # 2) Energy
    ax = axes[0][1]
    ax.plot(case.time(), case.finite("kinetic_energy"), label="kinetic E [J]")
    pot = case.finite("potential_energy")
    if pot:
        ax.plot(case.time(), pot, label="potential E [J]", alpha=0.7)
    ax.set_xlabel("t [s]")
    ax.set_ylabel("energy [J]")
    ax.set_title("Energy")
    ax.legend(fontsize=8)

    # 3) Contacts
    ax = axes[1][0]
    ax.plot(case.time(), case.finite("contacts"), color="tab:orange")
    ax.set_xlabel("t [s]")
    ax.set_ylabel("contacts")
    ax.set_title("Contact network")

    # 4) Maximum speed
    ax = axes[1][1]
    ax.plot(case.time(), case.finite("max_speed"), color="tab:green")
    ax.set_xlabel("t [s]")
    ax.set_ylabel("|v| max [m/s]")
    ax.set_title("Maximum speed")

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


def plot_sweep(cases: Sequence["Case"], out: str | Path) -> Path:
    """D vs impact energy, in log-log, with the 1/4 reference slope."""
    plt = _pyplot()
    from .models import compare_models

    result = compare_models(cases)
    points = result["points"]
    energies = [e for e in points["E_J"] if e is not None and e > 0]
    if points["n"] < 2 or len(energies) < points["n"]:
        raise SystemExit(
            "sweep plot needs 2+ cases with a positive launch energy "
            f"(got {points['n']}); run scripts/energy_sweep.py over the speeds first")
    out = Path(out)

    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    ax.plot(points["E_J"], points["D_mm"], "o-", label="measured")

    fit = result["power"]
    if "a" in fit:
        es = points["E_J"]
        ax.plot(es, [fit["a"] * e ** fit["b"] for e in es], "--",
                label=f"fit D=aE^b  (b={fit['b']:.3f}, R2={fit['R2']:.4f})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("impact E [J]")
    ax.set_ylabel("D [mm]")
    ax.set_title("Energy sweep")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


def plot_subsurface_chains(
    case: "Case",
    subsurface_csv: str | Path,
    out: str | Path,
) -> Path:
    """Visualization 1: Subsurface compaction saturation & 2D Force Chains (Heckel model).

    Panels:
      (a) 2D (r, z) Force Chains / Love-Weber Virial Stress sigma(r, z) [kPa] at peak stress
      (b) 2D (r, z) Subsurface Packing Fraction phi(r, z) vs Random Close Packing (0.64)
      (c) Time evolution of subsurface packing phi_subsurface(t) & Heckel exponential curve
      (d) Virial stress evolution vs projectile yield strength sigma_y & bond comminution
    """
    import math
    from .results import load_subsurface_chains

    plt = _pyplot()
    snapshots = load_subsurface_chains(subsurface_csv)
    if not snapshots:
        raise SystemExit(f"No subsurface data in {subsurface_csv}")

    # Pick the snapshot with the highest maximum virial stress (or final snapshot if all 0)
    peak_t = max(
        snapshots.keys(),
        key=lambda t: (max(c["stress_kPa"] for c in snapshots[t]), t),
    )
    final_t = max(snapshots.keys())

    def _to_grid(cells: list[dict[str, float]], key: str):
        rs = sorted({c["r_m"] * 1000.0 for c in cells})
        zs = sorted({c["z_m"] * 1000.0 for c in cells})
        r_idx = {r: i for i, r in enumerate(rs)}
        z_idx = {z: i for i, z in enumerate(zs)}
        grid = [[0.0 for _ in rs] for _ in zs]
        for c in cells:
            iz = z_idx[c["z_m"] * 1000.0]
            ir = r_idx[c["r_m"] * 1000.0]
            grid[iz][ir] = c[key]
        return rs, zs, grid

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    fig.suptitle(
        f"GranImpact · Subsurface Force Chains & Uniaxial Heckel Saturation ({case.name})",
        fontsize=12,
        fontweight="bold",
    )

    # Panel (a): 2D Force Chain Virial Stress sigma(r, z) [kPa]
    ax = axes[0][0]
    rs, zs, stress_grid = _to_grid(snapshots[peak_t], "stress_kPa")
    im1 = ax.imshow(
        stress_grid,
        origin="lower",
        aspect="auto",
        extent=[min(rs), max(rs), min(zs), max(zs)],
        cmap="inferno",
    )
    fig.colorbar(im1, ax=ax, label=r"Virial Stress $\sigma(r,z)$ [kPa]")
    ax.set_xlabel("Radial distance $r$ [mm]")
    ax.set_ylabel("Height $z$ [mm]")
    ax.set_title(f"(a) 2D Force Chains beneath Impact ($t = {peak_t:.4f}$ s)")

    # Panel (b): 2D Local Packing Fraction phi(r, z)
    ax = axes[0][1]
    rs_f, zs_f, phi_grid = _to_grid(snapshots[final_t], "phi_local")
    im2 = ax.imshow(
        phi_grid,
        origin="lower",
        aspect="auto",
        extent=[min(rs_f), max(rs_f), min(zs_f), max(zs_f)],
        cmap="viridis",
        vmin=0.40,
        vmax=0.64,
    )
    fig.colorbar(im2, ax=ax, label=r"Local Packing $\phi(r,z)$ [-]")
    ax.set_xlabel("Radial distance $r$ [mm]")
    ax.set_ylabel("Height $z$ [mm]")
    ax.set_title(f"(b) Subsurface Compaction Field ($t = {final_t:.4f}$ s)")

    # Panel (c): Subsurface packing evolution + Heckel analytical law
    ax = axes[1][0]
    if "phi_subsurface" in case.data and case.finite("phi_subsurface"):
        ts = [t for t, v in zip(case.data["t"], case.data["phi_subsurface"]) if v is not None]
        phis = case.finite("phi_subsurface")
        ax.plot(ts, phis, "o-", color="tab:blue", label=r"Measured $\phi_{\mathrm{subsurface}}(t)$")
    phi0 = float(case.meta.get("settled_phi") or case.meta.get("phi") or 0.524)
    ax.axhline(phi0, linestyle="--", color="gray", label=f"Initial bed $\\phi_0 = {phi0:.3f}$")
    ax.axhline(0.64, linestyle=":", color="tab:red", label=r"RCP saturation $\phi_{\max} = 0.64$")
    ax.set_xlabel("t [s]")
    ax.set_ylabel(r"Packing fraction $\phi$ [-]")
    ax.set_title("(c) Uniaxial Compaction Saturation vs RCP Limit")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel (d): Virial stress & Cohesive Bond Comminution (Regime 1 vs Regime 2)
    ax = axes[1][1]
    if "max_virial_stress_kPa" in case.data and case.finite("max_virial_stress_kPa"):
        ts = [t for t, v in zip(case.data["t"], case.data["max_virial_stress_kPa"]) if v is not None]
        s_max = case.finite("max_virial_stress_kPa")
        s_mean = case.finite("mean_virial_stress_kPa")
        ax.plot(ts, s_max, "-", color="tab:red", label=r"Peak stress $\sigma_{\max}$ [kPa]")
        if len(s_mean) == len(ts):
            ax.plot(ts, s_mean, "--", color="tab:orange", label=r"Mean column stress $\langle\sigma\rangle$ [kPa]")
    ax.axhline(17.15, linestyle=":", color="black", label=r"Projectile yield $\sigma_y = 17.15$ kPa")
    ax.set_xlabel("t [s]")
    ax.set_ylabel("Stress [kPa]")
    ax.set_title(r"(d) Subsurface Stress vs Projectile Yield $\sigma_y$")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.3)

    out_path = Path(out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def plot_lateral_saturation(case: "Case", out: str | Path) -> Path:
    """Visualization 2: Lateral opening & Rim height saturation mechanism.

    Panels:
      (a) Vertical-to-horizontal momentum transfer kappa_vh(t) & Diameter D(t) (Granular Bulldozing)
      (b) Rim height h_rim(t) & Inner wall slope tan(alpha)(t) vs Angle of Repose limit
      (c) Volume non-conservation: V_in(t) vs V_ex(t) vs Excess Compaction Volume V_exc(t)
      (d) Energy partitioning: Translational KE, Rotational KE, and Cumulative Frictional Work W_fric
    """
    import math

    plt = _pyplot()
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    fig.suptitle(
        f"GranImpact · Lateral Opening & Rim Saturation Theory ({case.name})",
        fontsize=12,
        fontweight="bold",
    )

    def _paired(col: str, scale: float = 1.0):
        if col not in case.data:
            return [], []
        ts, vs = [], []
        for t, v in zip(case.data["t"], case.data[col]):
            if t is not None and v is not None and math.isfinite(v):
                ts.append(t)
                vs.append(v * scale)
        return ts, vs

    # Panel (a): Momentum transfer kappa_vh & Lateral opening D(t)
    ax = axes[0][0]
    t_d, d_mm = _paired("D", 1000.0)
    t_k, kappa = _paired("momentum_transfer_vh", 1.0)
    if d_mm:
        ax.plot(t_d, d_mm, "o-", color="tab:blue", label="Crater diameter $D(t)$ [mm]")
    ax.set_xlabel("t [s]")
    ax.set_ylabel("$D(t)$ [mm]", color="tab:blue")
    ax.grid(True, alpha=0.3)
    if kappa:
        ax2 = ax.twinx()
        ax2.plot(t_k, kappa, "--", color="tab:red", label=r"Momentum transfer $\kappa_{vh}(t)$")
        ax2.set_ylabel(r"$\kappa_{vh} = P_{\mathrm{horiz}} / P_{0,z}$ [-]", color="tab:red")
    ax.set_title("(a) Lateral Opening & Horizontal Momentum Transfer")

    # Panel (b): Rim height saturation & Inner wall slope vs Angle of Repose
    ax = axes[0][1]
    t_r, h_rim_mm = _paired("h_rim", 1000.0)
    t_s, slope = _paired("rim_slope", 1.0)
    if h_rim_mm:
        ax.plot(t_r, h_rim_mm, "s-", color="tab:green", label=r"Rim height $h_{\mathrm{rim}}(t)$ [mm]")
    h_model = case.meta.get("h_rim_model_m")
    if h_model is not None:
        ax.axhline(
            float(h_model) * 1000.0,
            linestyle=":",
            color="darkgreen",
            label=f"Log-Log model $h_{{\\mathrm{{rim}}}} = {float(h_model)*1000.0:.2f}$ mm",
        )
    ax.set_xlabel("t [s]")
    ax.set_ylabel(r"$h_{\mathrm{rim}}$ [mm]", color="tab:green")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.3)
    if slope:
        ax2 = ax.twinx()
        ax2.plot(t_s, slope, "--", color="tab:purple", label=r"Inner wall slope $\tan\theta(t)$")
        ax2.axhline(math.tan(math.radians(40.0)), linestyle=":", color="tab:purple", alpha=0.6)
        ax2.set_ylabel(r"Inner wall slope $\tan\theta$ [-]", color="tab:purple")
    ax.set_title(r"(b) Rim Height Saturation & Angle-of-Repose Limit")

    # Panel (c): Volume balance V_in vs V_ex vs V_exc (subsurface compression)
    ax = axes[1][0]
    t_vin, v_in = _paired("v_in", 1e6)
    t_vex, v_ex = _paired("v_ex", 1e6)
    t_vexc, v_exc = _paired("v_exc", 1e6)
    if v_in:
        ax.plot(t_vin, v_in, "-", color="tab:blue", label=r"Cavity volume $V_{\mathrm{in}}$ [cm$^3$]")
    if v_ex:
        ax.plot(t_vex, v_ex, "--", color="tab:orange", label=r"Rim volume $V_{\mathrm{ex}}$ [cm$^3$]")
    if v_exc:
        ax.plot(t_vexc, v_exc, "-.", color="tab:red", label=r"Compacted volume $V_{\mathrm{exc}}$ [cm$^3$]")
    ax.set_xlabel("t [s]")
    ax.set_ylabel(r"Volume [cm$^3$]")
    ax.set_title(r"(c) Volume Balance ($V_{\mathrm{exc}} = V_{\mathrm{in}} - V_{\mathrm{ex}}$)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel (d): Energy Partitioning & Frictional Dissipation W_fric
    ax = axes[1][1]
    t_ke, ke = _paired("kinetic_energy", 1.0)
    t_rot, rke = _paired("rotational_energy", 1.0)
    t_dis, w_fric = _paired("dissipated_energy", 1.0)
    if ke:
        ax.plot(t_ke, ke, "-", color="tab:blue", label="Translational $E_k(t)$ [J]")
    if rke:
        ax.plot(t_rot, rke, "--", color="tab:cyan", label="Rotational $E_{\\mathrm{rot}}(t)$ [J]")
    if w_fric:
        ax.plot(t_dis, w_fric, "-", color="tab:red", label="Frictional dissipation $W_{\\mathrm{fric}}(t)$ [J]")
    ax.set_xlabel("t [s]")
    ax.set_ylabel("Energy [J]")
    ax.set_title("(d) Frictional & Rotational Energy Dissipation")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    out_path = Path(out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def plot_crater_formation(
    case: "Case",
    radial_csv: str | Path,
    depth_map_csv: str | Path,
    out: str | Path,
) -> Path:
    """Visualization 3: Time-evolution of crater formation, central peak & 2D ToF depth map.

    Panels:
      (a) Radial profile time-evolution z_bed(r, t) during impact
      (b) Final radial cross-section: Pure Bed Cavity z_bed(r) vs Total Surface z_total(r) (Porous Piston)
      (c) 2D Topographic Depth Map delta_z(x, y) [mm] (Time-of-Flight camera view)
    """
    from .results import load_depth_map, load_radial_evolution

    plt = _pyplot()
    rad_snaps = load_radial_evolution(radial_csv)
    dm = load_depth_map(depth_map_csv)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    fig.suptitle(
        f"GranImpact · Crater Formation Process & Porous Piston ({case.name})",
        fontsize=12,
        fontweight="bold",
    )

    # Panel (a): Radial profile evolution z_bed(r, t)
    ax = axes[0]
    times = sorted(rad_snaps.keys())
    cmap = plt.get_cmap("plasma")
    for idx, t in enumerate(times):
        color = cmap(idx / max(1, len(times) - 1))
        rs_mm = [b["r_m"] * 1000.0 for b in rad_snaps[t]]
        zs_mm = [b["z_bed_mm"] for b in rad_snaps[t]]
        ax.plot(rs_mm, zs_mm, color=color, label=f"t = {t:.3f} s")
    ax.axhline(0.0, linestyle="--", color="gray", alpha=0.6)
    ax.set_xlabel("Radial distance $r$ [mm]")
    ax.set_ylabel(r"Bed elevation $\Delta z_{\mathrm{bed}}(r)$ [mm]")
    ax.set_title("(a) Radial Excavation Evolution")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3)

    # Panel (b): Final Bed Cavity vs Total Surface (Central Peak / Porous Piston)
    ax = axes[1]
    if times:
        last_bins = rad_snaps[times[-1]]
        rs_mm = [b["r_m"] * 1000.0 for b in last_bins]
        z_bed = [b["z_bed_mm"] for b in last_bins]
        z_tot = [b["z_total_mm"] for b in last_bins]
        ax.plot(rs_mm, z_bed, "--", color="tab:blue", linewidth=2.0, label="Bed cavity ($z_{\\mathrm{bed}}$)")
        ax.plot(rs_mm, z_tot, "-", color="tab:red", linewidth=2.0, label="Total surface ($z_{\\mathrm{total}}$)")
        ax.fill_between(
            rs_mm,
            z_bed,
            z_tot,
            where=[zt >= zb for zb, zt in zip(z_bed, z_tot)],
            color="tab:orange",
            alpha=0.35,
            label="Porous piston / Central peak",
        )
    ax.axhline(0.0, linestyle=":", color="black", alpha=0.6)
    ax.set_xlabel("Radial distance $r$ [mm]")
    ax.set_ylabel(r"Elevation $\Delta z(r)$ [mm]")
    ax.set_title("(b) Final Cavity & Central Peak Cross-Section")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel (c): 2D Synthetic ToF Depth Map
    ax = axes[2]
    xs = dm["x_mm"]
    ys = dm["y_mm"]
    z2d = dm["z_total_mm"]
    im = ax.imshow(
        z2d,
        origin="lower",
        extent=[min(xs), max(xs), min(ys), max(ys)],
        cmap="RdBu_r",
    )
    fig.colorbar(im, ax=ax, label=r"Elevation $\Delta z(x,y)$ [mm]")
    ax.set_xlabel("$x$ [mm]")
    ax.set_ylabel("$y$ [mm]")
    ax.set_title("(c) 2D ToF Depth Map (craterslab compatible)")

    out_path = Path(out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path

