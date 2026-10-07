# granimpact_post

Post-processing of GranImpact results. It only uses the standard library to read
and summarize; `matplotlib` is optional and only imported when plots are requested.

```bash
python3 -m granimpact_post summary ../data/outputs/default/default_crater_evolution.csv
python3 -m granimpact_post list      "../data/outputs/*/*.csv"
python3 -m granimpact_post sweep     "../data/outputs/energy_sweep_*/*.csv" -o table.csv --plot D_vs_E.png
python3 -m granimpact_post fit       "../data/outputs/energy_sweep_*/*.csv"
python3 -m granimpact_post plot      ../data/outputs/default/default_crater_evolution.csv -o default.png
```

Install (optional; it can also be used with `PYTHONPATH=python`):

```bash
pip install --user .          # o:  pip install --user ".[plots]"
```

## Modules

| Module | Content |
|---|---|
| `results.py` | `Case`, `load_case`, `load_many`, `sweep_table`: reading the self-contained CSV |
| `mass.py` | projectile mass inferred from the energy series (not from nominal values) |
| `models.py` | D = a·E^b and D = a + b·ln E fits on the measured data |
| `plots.py` | plots (requires matplotlib) |
| `cli.py` | command-line interface |

Method note: `models.py` ships no literature coefficients; it fits what is in the
CSVs. The reference exponents (Uehara's b ≈ 1/4, the thesis logarithmic law) are
printed as hypotheses to be tested, and the verdict comes from the R² of the data.
