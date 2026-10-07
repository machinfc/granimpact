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

    # ------------------------------------------------------------- momento clave --
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

        classification = self.labels.get("classification", [])
        return {
            "case": self.name,
            "t_final_s": last("t"),
            "impact_speed_m_s": self.meta.get("impact_speed_m_s"),
            "D_mm": _mm(last("D")),
            "d_exc_mm": _mm(last("d_exc")),
            "h_rim_mm": _mm(last("h_rim")),
            "epsilon": last("epsilon"),
            "Z/D": last("aspect_zd"),
            "V_in_cm3": _cm3(last("v_in")),
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
    """Launch kinetic energy = 1/2 m v^2, with the mass measured in the run.

    Computed from the CSV (not from nominal values): the projectile is a
    set of grains generated with a seed, so its mass is a result.
    """
    from .mass import projectile_mass_kg  # local import: avoids circular dependencies

    mass = projectile_mass_kg(case)
    speed = case.meta.get("v")
    if mass is None or speed is None:
        return None
    return 0.5 * mass * float(speed) ** 2
