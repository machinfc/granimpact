"""Reading and summarizing the time series written by the simulator.

Expected format (written by io::CsvWriter):

    # GranImpact default | dt=5e-06 t_settle=0.15 t_end=0.08 phi=0.524 v=2 angle=90 gz=-9.81
    t,particles,contacts,kinetic_energy,potential_energy,max_speed,d_exc,D,epsilon,...
    0.0,21904,0,0.0,0.0,0.0,,,...

Design key: the CSV is self-contained (the header carries the run parameters
case), so a results folder can be analyzed without access to the configs.
"""

from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, Sequence

HEADER_PREFIX = "# GranImpact"
FOOTER_PREFIX = "# result:"
_META_RE = re.compile(r"(\w+)=([^\s]+)")

# Columns the analyzer leaves empty until there is a crater.
CRATER_COLUMNS = (
    "d_exc",
    "D",
    "epsilon",
    "h_rim",
    "aspect_zd",
    "v_in",
    "ejecta_mass",
)


def _to_number(text: str) -> float | None:
    text = (text or "").strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


@dataclass
class Case:
    """A simulated case: header parameters + time series."""

    name: str
    path: Path | None = None
    meta: dict[str, str] = field(default_factory=dict)
    columns: list[str] = field(default_factory=list)
    data: dict[str, list[float | None]] = field(default_factory=dict)
    labels: dict[str, list[str]] = field(default_factory=dict)  # text columns
    rows: int = 0

    # ------------------------------------------------------------------ series --
    def series(self, column: str) -> list[float | None]:
        if column not in self.data:
            raise KeyError(f"column '{column}' does not exist; available: {', '.join(self.columns)}")
        return self.data[column]

    def finite(self, column: str) -> list[float]:
        """Series without the empty stretches (e.g. before the impact)."""
        return [v for v in self.data[column] if v is not None and math.isfinite(v)]

    def time(self) -> list[float]:
        return self.finite("t")

    # ------------------------------------------------------------- key moment -----
    def impact_row(self) -> int | None:
        """First row where the projectile already touches something, i.e. the
        start of the impact. If no contacts are recorded, the first row with
        morphometry is computed. `None` if the projectile has not reached the bed."""
        for i, value in enumerate(self.data.get("contacts", [])):
            if value:
                return i
        for i, value in enumerate(self.data.get("D", [])):
            if value is not None:
                return i
        return None

    def param(self, key: str, cast=float):
        """Header parameter converted to the requested type."""
        if key not in self.meta:
            raise KeyError(f"the header has no '{key}'; it has: {', '.join(self.meta)}")
        return cast(self.meta[key])

    # ------------------------------------------------------------------ summary --
    def summary(self) -> dict[str, float | str | None]:
        """Final case result (last row) + maxima of interest."""
        impact = self.impact_row()

        def last(column: str):
            values = self.finite(column)
            return values[-1] if values else None

        def peak(column: str):
            values = self.finite(column)
            return max(values) if values else None

        def opt_last(column: str):
            if column not in self.data:
                return None
            values = self.finite(column)
            return values[-1] if values else None

        def opt_peak(column: str):
            if column not in self.data:
                return None
            values = self.finite(column)
            return max(values) if values else None

        classification = self.labels.get("classification", [])
        return {
            "case": self.name,
            "t_final_s": last("t"),
            "impact_speed_m_s": self.meta.get("impact_speed_m_s"),
            "D_mm": _mm(last("D")),
            "d_exc_mm": _mm(last("d_exc")),
            "h_rim_mm": _mm(last("h_rim")),
            "h_cp_mm": _mm(opt_last("h_cp")),
            "rim_slope": opt_last("rim_slope"),
            "epsilon": last("epsilon"),
            "Z/D": last("aspect_zd"),
            "V_in_cm3": _cm3(last("v_in")),
            "V_exc_cm3": _cm3(opt_last("v_exc")),
            "phi_subsurface": opt_last("phi_subsurface"),
            "max_stress_kPa": opt_peak("max_virial_stress_kPa"),
            "class": classification[-1] if classification else None,
            "KE_max_J": peak("kinetic_energy"),
            "max_contacts": peak("contacts"),
            "rows": self.rows,
            "impact_row": impact,
        }

    def as_text(self) -> str:
        s = self.summary()
        width = max(len(k) for k in s)
        return "\n".join(f"{k.ljust(width)} : {v}" for k, v in s.items())


def _mm(value: float | None) -> float | None:
    return None if value is None else value * 1000.0


def _cm3(value: float | None) -> float | None:
    return None if value is None else value * 1e6


def load_case(path: str | Path) -> Case:
    """Reads a GranImpact CSV."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)

    case = Case(name=path.stem, path=path)
    with path.open(encoding="utf-8", newline="") as handle:
        header = handle.readline().strip()
        if header.startswith(HEADER_PREFIX):
            body = header[len(HEADER_PREFIX):].strip()
            name, _, rest = body.partition("|")
            case.name = name.strip() or case.name
            case.meta = {k: v for k, v in _META_RE.findall(rest)}
            column_line = handle.readline()
        else:  # CSV without its own header: the first line is assumed to be columns
            column_line = header + "\n"

        reader = csv.reader(_chain_first(column_line, handle))
        header_row = next(reader, [])
        if not header_row:
            raise ValueError(f"{path} has no column header")
        case.columns = [c.strip() for c in header_row]
        case.data = {c: [] for c in case.columns}
        case.labels = {c: [] for c in case.columns}

        for row in reader:
            if row and row[0].startswith(FOOTER_PREFIX.strip()):
                # Final block written when the run ends: measured metrics.
                case.meta.update(_META_RE.findall(row[0]))
                continue
            if not row or not any(row):
                continue
            for column, raw in zip(case.columns, row):
                number = _to_number(raw)
                if number is None:
                    case.labels[column].append(raw.strip())
                    case.data[column].append(None)
                else:
                    case.data[column].append(number)
                    case.labels[column].append("")
            case.rows += 1
    return case


def _chain_first(first: str, rest: Iterator[str]) -> Iterator[str]:
    yield first
    yield from rest


def load_many(paths: Iterable[str | Path]) -> list[Case]:
    return [load_case(p) for p in sorted(Path(p) for p in paths)]


def sweep_table(cases: Sequence[Case]) -> list[dict[str, float | str | None]]:
    """Comparison table ready to print or plot (energy sweep)."""
    rows = []
    for case in cases:
        s = case.summary()
        energy = _launch_energy(case)
        rows.append(
            {
                "case": case.name,
                "v_m_s": case.meta.get("v"),
                "E_J": energy,
                "D_mm": s["D_mm"],
                "d_exc_mm": s["d_exc_mm"],
                "h_rim_mm": s["h_rim_mm"],
                "Z/D": s["Z/D"],
                "eps": s["epsilon"],
                "class": s["class"],
                "path": str(case.path) if case.path else "",
            }
        )
    return rows


def _launch_energy(case: Case) -> float | None:
    """Launch kinetic energy of the projectile, in joules."""
    from .mass import projectile_mass_kg  # local import: avoids circular dependencies

    footer = case.meta.get("launch_energy_J")
    if footer is not None:
        try:
            return float(footer)
        except ValueError:
            pass
    mass = projectile_mass_kg(case)
    speed = case.meta.get("v")
    if mass is None or speed is None:
        return None
    return 0.5 * mass * float(speed) ** 2


def resolve_sibling_csv(case_or_path: Case | str | Path, suffix: str) -> Path:
    """Finds `<name>_<suffix>.csv` next to a `Case` or `<name>_crater_evolution.csv`."""
    if isinstance(case_or_path, Case):
        if case_or_path.path is None:
            raise FileNotFoundError(f"Case '{case_or_path.name}' has no file path")
        base_dir = case_or_path.path.parent
        name = case_or_path.name
    else:
        p = Path(case_or_path)
        if p.is_file() and p.name.endswith(f"_{suffix}.csv"):
            return p
        base_dir = p.parent
        name = p.stem.replace("_crater_evolution", "")
    candidate = base_dir / f"{name}_{suffix}.csv"
    if not candidate.is_file():
        raise FileNotFoundError(f"Companion file not found: {candidate}")
    return candidate


def load_depth_map(path: str | Path) -> dict:
    """Loads `<name>_depth_map.csv` into a structured dict with 2D arrays (in mm)."""
    path = Path(path)
    meta: dict[str, float] = {}
    rows: list[tuple[int, int, float, float, float, float]] = []
    max_ix = 0
    max_iy = 0
    with path.open(encoding="utf-8", newline="") as handle:
        first = handle.readline().strip()
        if first.startswith("#"):
            for k, v in _META_RE.findall(first):
                num = _to_number(v)
                if num is not None:
                    meta[k] = num
            header_line = handle.readline()
        else:
            header_line = first + "\n"
        reader = csv.DictReader(_chain_first(header_line, handle))
        for row in reader:
            ix = int(row["ix"])
            iy = int(row["iy"])
            x_mm = float(row["x_mm"])
            y_mm = float(row["y_mm"])
            zb = float(row["z_bed_mm"])
            zt = float(row["z_total_mm"])
            max_ix = max(max_ix, ix)
            max_iy = max(max_iy, iy)
            rows.append((ix, iy, x_mm, y_mm, zb, zt))

    nx = max_ix + 1
    ny = max_iy + 1
    z_bed = [[0.0 for _ in range(nx)] for _ in range(ny)]
    z_total = [[0.0 for _ in range(nx)] for _ in range(ny)]
    x_coords = [0.0 for _ in range(nx)]
    y_coords = [0.0 for _ in range(ny)]
    for ix, iy, x_mm, y_mm, zb, zt in rows:
        z_bed[iy][ix] = zb
        z_total[iy][ix] = zt
        x_coords[ix] = x_mm
        y_coords[iy] = y_mm

    return {
        "meta": meta,
        "nx": nx,
        "ny": ny,
        "x_mm": x_coords,
        "y_mm": y_coords,
        "z_bed_mm": z_bed,
        "z_total_mm": z_total,
    }


def export_craterslab_npz(
    depth_map_csv: str | Path,
    out_npz: str | Path | None = None,
    use_total_surface: bool = True,
) -> Path:
    """Exports `<name>_depth_map.csv` to a `.npz` file compatible with `craterslab`.

    `craterslab.sensors.DepthMap` loads `.npz` arrays where depth/elevation is in mm
    along with `SensorResolution(x_res, y_res, z_res, 'mm')`.
    """
    try:
        import numpy as np
    except ImportError as exc:  # pragma: no cover
        raise SystemExit("numpy is required to export .npz files: pip install numpy") from exc

    dm = load_depth_map(depth_map_csv)
    src_path = Path(depth_map_csv)
    if out_npz is None:
        stem = src_path.stem.replace("_depth_map", "")
        out_npz = src_path.parent / f"{stem}_craterslab.npz"
    out_path = Path(out_npz)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    surface_key = "z_total_mm" if use_total_surface else "z_bed_mm"
    dmap = np.asarray(dm[surface_key], dtype=np.float32)
    z_bed = np.asarray(dm["z_bed_mm"], dtype=np.float32)
    z_tot = np.asarray(dm["z_total_mm"], dtype=np.float32)
    dx_mm = float(dm["meta"].get("dx_mm", 2.8025))
    dy_mm = float(dm["meta"].get("dy_mm", 2.8025))

    np.savez_compressed(
        out_path,
        dmap=dmap,
        z_bed_mm=z_bed,
        z_total_mm=z_tot,
        x_mm=np.asarray(dm["x_mm"], dtype=np.float32),
        y_mm=np.asarray(dm["y_mm"], dtype=np.float32),
        resolution_mm=np.asarray([dx_mm, dy_mm, 1.0], dtype=np.float32),
    )
    return out_path


def load_subsurface_chains(path: str | Path) -> dict[float, list[dict[str, float]]]:
    """Loads `<name>_subsurface_chains.csv` grouped by snapshot time `t`."""
    snapshots: dict[float, list[dict[str, float]]] = {}
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            t = round(float(row["t"]), 8)
            snapshots.setdefault(t, []).append(
                {
                    "r_m": float(row["r_m"]),
                    "z_m": float(row["z_m"]),
                    "phi_local": float(row["phi_local"]),
                    "stress_kPa": float(row["stress_kPa"]),
                    "coordination": float(row["coordination"]),
                    "speed_m_s": float(row["speed_m_s"]),
                }
            )
    return snapshots


def load_radial_evolution(path: str | Path) -> dict[float, list[dict[str, float]]]:
    """Loads `<name>_radial_evolution.csv` grouped by snapshot time `t`."""
    snapshots: dict[float, list[dict[str, float]]] = {}
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            t = round(float(row["t"]), 8)
            snapshots.setdefault(t, []).append(
                {
                    "r_m": float(row["r_m"]),
                    "z_bed_mm": float(row["z_bed_mm"]),
                    "z_total_mm": float(row["z_total_mm"]),
                }
            )
    return snapshots

