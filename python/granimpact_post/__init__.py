"""GranImpact post-processing: read the CSV files, summarize and plot them.

Only the standard library is needed to read and summarize; matplotlib is optional
and it is imported only when plots are requested.

Typical use::

    from granimpact_post import load_case, sweep_table
    case = load_case("data/outputs/default/default_crater_evolution.csv")
    print(case.summary())
"""

from .results import (  # noqa: F401
    Case,
    export_craterslab_npz,
    load_case,
    load_depth_map,
    load_many,
    load_radial_evolution,
    load_subsurface_chains,
    sweep_table,
)

__all__ = [
    "Case",
    "export_craterslab_npz",
    "load_case",
    "load_depth_map",
    "load_many",
    "load_radial_evolution",
    "load_subsurface_chains",
    "sweep_table",
]
__version__ = "0.2.0"
