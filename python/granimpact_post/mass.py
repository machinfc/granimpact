"""Projectile mass from the data left in the output.

The projectile has no tabulated mass: it is a sphere of grains generated with a
seed. The exact mass is recovered from the initial kinetic energy of the series:
in the first row the bed is at rest and only the projectile moves, so
    m = 2 * KE / v^2 (v = launch speed from the header).
If the first row already shows bed motion, the estimate is discarded.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from .results import Case


def projectile_mass_kg(case: "Case") -> float | None:
    speed = case.meta.get("v")
    if speed is None or float(speed) <= 0:
        return None
    kinetic = case.data.get("kinetic_energy", [])
    contacts = case.data.get("contacts", [])
    for i, ke in enumerate(kinetic):
        if ke is None:
            continue
        # Only valid while the bed touches nothing (contacts = 0).
        if i < len(contacts) and contacts[i]:
            continue
        return 2.0 * ke / float(speed) ** 2
    return None
