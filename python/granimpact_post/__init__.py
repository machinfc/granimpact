"""GranImpact post-processing: read the CSV files, summarize and plot them.

Only the standard library is needed to read and summarize; matplotlib is optional
and it is imported only when plots are requested.

Typical use::

    from granimpact_post import load_case, sweep_table
    case = load_case("data/outputs/default/default_crater_evolution.csv")
    print(case.summary())
"""

from .results import Case, load_case, load_many, sweep_table  # noqa: F401

__all__ = ["Case", "load_case", "load_many", "sweep_table"]
__version__ = "0.1.0"
