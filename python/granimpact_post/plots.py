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
