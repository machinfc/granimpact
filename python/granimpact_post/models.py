"""Fit of D(E) to the two models compared in this work.

  * power law D = a * E^b (Uehara et al. propose b ~ 1/4)
  * logarithmic law D = a + b * ln E (the thesis model, for slow impacts)

Here only what is in the CSV is FITTED: this module ships no coefficients
preset nor expected results. The reference exponents are printed as
hypothesis to test, and the verdict comes from the R2 of the user's data.

Practical advantage: the energy sweep is six long runs, so the fit
is done on the measurements, not on numbers from the literature.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Iterable, Sequence

if TYPE_CHECKING:  # pragma: no cover
    from .results import Case


def _points(cases: Iterable["Case"]) -> tuple[list[float], list[float], list[str]]:
    from .results import sweep_table

    energies, diameters, names = [], [], []
    for row in sweep_table(cases):
        if row["E_J"] is None or row["D_mm"] is None:
            continue
        energies.append(float(row["E_J"]))
        diameters.append(float(row["D_mm"]))
        names.append(str(row["case"]))
    return energies, diameters, names


def _r_squared(observed: Sequence[float], predicted: Sequence[float]) -> float:
    mean = sum(observed) / len(observed)
    ss_tot = sum((o - mean) ** 2 for o in observed)
    ss_res = sum((o - p) ** 2 for o, p in zip(observed, predicted))
    return float("nan") if ss_tot == 0 else 1.0 - ss_res / ss_tot


def fit_power_law(energies: Sequence[float], diameters: Sequence[float]) -> dict:
    """D = a*E^b by least squares in log-log."""
    if len(energies) < 2:
        return {"model": "power", "n": len(energies), "note": "2+ points needed"}
    xs = [math.log(e) for e in energies]
    ys = [math.log(d) for d in diameters]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    b = sxy / sxx if sxx else float("nan")
    a = math.exp(my - b * mx)
    predicted = [a * e**b for e in energies]
    return {
        "model": "D = a*E^b",
        "n": n,
        "a": a,
        "b": b,
        "b_esperado_uehara": 0.25,
        "R2": _r_squared(diameters, predicted),
    }


def fit_log(energies: Sequence[float], diameters: Sequence[float]) -> dict:
    """D = a + b*ln E by least squares."""
    if len(energies) < 2:
        return {"model": "logarithm", "n": len(energies), "note": "2+ points needed"}
    xs = [math.log(e) for e in energies]
    ys = list(diameters)
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    b = sxy / sxx if sxx else float("nan")
    a = my - b * mx
    predicted = [a + b * x for x in xs]
    return {"model": "D = a + b*ln E", "n": n, "a": a, "b": b, "R2": _r_squared(ys, predicted)}


def compare_models(cases: Iterable["Case"]) -> dict[str, dict]:
    """Fits both models and also returns which one wins on R2."""
    energies, diameters, names = _points(cases)
    result = {
        "points": {"n": len(energies), "cases": ", ".join(names) or "-",
                   "E_J": energies, "D_mm": diameters},
        "power": fit_power_law(energies, diameters),
        "logarithm": fit_log(energies, diameters),
    }
    r2_pow = result["power"].get("R2")
    r2_log = result["logarithm"].get("R2")
    if r2_pow is not None and r2_log is not None and not math.isnan(r2_pow) and not math.isnan(r2_log):
        result["best_by_R2"] = {"model": "power" if r2_pow >= r2_log else "logarithm",
                                  "R2_power": r2_pow, "R2_logarithm": r2_log,
                                  "note": "comparison against the few sweep points: "
                                          "indicative, not conclusive"}
    return result
