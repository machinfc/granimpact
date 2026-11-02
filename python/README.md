# granimpact_post — Scientific Post-Processing & `craterslab` Exporter

`granimpact_post` is the Python post-processing package for **GranImpact**. It parses the
self-contained simulation CSVs, generates the 3 publication-grade validation figures
(subsurface Love–Weber force chains & Heckel saturation, lateral momentum redirection & rim
height saturation, and radial crater evolution with porous-piston anatomy), fits analytical
scaling laws, and exports 2D elevation maps to `.npz` files directly compatible with
`craterslab`.

---

## 1. Installation & Usage

Install in editable mode from the repository root (or run directly with `PYTHONPATH=python`):

```bash
pip install -e ./python
```

### Generate All 3 Scientific Visualizations + `craterslab` `.npz` (`visualize-all`)

```bash
# Loose bed cases
python3 -m granimpact_post visualize-all \
    data/outputs/loose/default/default_crater_evolution.csv

python3 -m granimpact_post visualize-all \
    data/outputs/loose/loose_vertical/loose_vertical_crater_evolution.csv

# Compacted bed cases
python3 -m granimpact_post visualize-all \
    data/outputs/compact/compact_vertical/compact_vertical_crater_evolution.csv
```

### Individual CLI Subcommands

```bash
# 1. Visualization 1: 2D (r,z) subsurface force chains & Heckel RCP saturation
python3 -m granimpact_post chains \
    data/outputs/loose/loose_vertical/loose_vertical_crater_evolution.csv

# 2. Visualization 2: Lateral opening, momentum transfer & rim height saturation
python3 -m granimpact_post rim-theory \
    data/outputs/loose/loose_vertical/loose_vertical_crater_evolution.csv

# 3. Visualization 3: Radial crater formation, porous piston & 2D ToF depth map
python3 -m granimpact_post formation \
    data/outputs/loose/loose_vertical/loose_vertical_crater_evolution.csv

# 4. Export 2D ToF elevation map directly to craterslab-compatible .npz
python3 -m granimpact_post export-npz \
    data/outputs/loose/loose_vertical/loose_vertical_depth_map.csv \
    -o data/outputs/loose/loose_vertical/loose_vertical_craterslab.npz

# 5. Inspect case summaries, energy sweep tables, and D(E) scaling fits
python3 -m granimpact_post summary "data/outputs/*/*/*_crater_evolution.csv"
python3 -m granimpact_post list    "data/outputs/*/*/*_crater_evolution.csv"
python3 -m granimpact_post sweep   "data/outputs/loose/energy_sweep_*/*.csv" \
    -o data/outputs/loose/sweep_table.csv --plot data/outputs/loose/D_vs_E.png
python3 -m granimpact_post fit     "data/outputs/loose/energy_sweep_*/*.csv"
```

---

## 2. Package Modules & `craterslab` `.npz` Schema

| Module | Responsibility |
| :--- | :--- |
| `results.py` | `Case`, `load_case`, `load_many`, `sweep_table`, `export_craterslab_npz`: self-contained CSV & depth-map parser |
| `plots.py` | Generates `_force_chains.png`, `_lateral_saturation.png`, `_crater_formation.png`, `_overview.png`, and energy sweep plots |
| `models.py` | Least-squares fits of power-law (`D = a E^b`) and logarithmic (`D = a + b ln E`) models with `R²` comparison |
| `mass.py` | Infers effective projectile mass from the initial kinetic energy series |
| `cli.py` | Command-line interface exposing all 10 subcommands |

### `craterslab` `.npz` Array Keys

Every `<case>_craterslab.npz` file exported by `export-npz` or `visualize-all` contains:
- `xp`, `yp`: 1D coordinate axes in millimeters (`mm`), centered on the container (`[-225, +225] mm`).
- `M`: 2D elevation matrix `Δz(y, x)` in millimeters (`mm`, `+` = rim / central peak, `-` = excavated cavity).
- `x_mm`, `y_mm`, `z_total_mm`, `z_bed_mm`: Full 2D meshes separating pure bed excavation (`z_bed_mm`) from the total topography including the comminuted projectile (`z_total_mm`).
