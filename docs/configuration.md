# Configuration, cases and outputs

## 1. JSON schema

`SimConfig::fromFile` + `validate()`. SI units throughout (metres, seconds, kilograms,
pascals). Unknown keys are ignored; any value out of range aborts at startup with a
concrete message.

| Block | Key | Unit | Code default | Meaning |
|---|---|---|---|---|
| `simulation` | `name` | — | `default` | prefix of the output files |
| | `dt` | s | 1e-6 | integration step (see the Rayleigh limit) |
| | `t_settle` | s | 0.02 | maximum settling budget |
| | `t_end` | s | 0.05 | duration of the impact phase |
| | `output_every` | steps | 500 | write a CSV row every N steps |
| | `seed` | — | 20240517 | seed (all randomness goes through it) |
| `domain` | `x`,`y`,`z` | m | 0.45, 0.45, 0.25 | experimental box 45×45×25 cm |
| `particles` | `r_min`,`r_max` | m | 0.4, 0.6 mm | sand radius range |
| | `rho_grain` | kg/m³ | 2650 | grain density (quartz) |
| | `packing_fraction` | — | 0.524 | density **target**; the real one is measured |
| | `bed_height_fraction` | — | 0.60 | fraction of the box height where the bed is seeded |
| | `young` | Pa | 70e9 | grain Young's modulus (the configs use 1e8, see below) |
| | `poisson` | — | 0.25 | Poisson's ratio |
| | `restitution` | — | 0.55 | grain–grain restitution coefficient |
| | `friction` | — | 0.40 | Coulomb friction |
| `projectile` | `radius` | m | 0.0354 | projectile radius (7.07 cm diameter) |
| | `density` | kg/m³ | 1310 | effective density of the granular projectile |
| | `phi` | — | 0.50 | internal packing fraction |
| | `young`,`poisson` | Pa, — | 5e6, 0.20 | projectile properties |
| | `yield_kPa` | kPa | 17.15 | yield stress measured in the thesis |
| | `speed` | m/s | 2.0 | initial speed |
| | `angle_deg` | degrees | 90 | 90 = vertical; 45 = oblique |
| | `start_z` | m | 0.20 | projectile centre height |
| `contact` | `tangential` | — | true | enables Mindlin + Coulomb |
| | `rolling_friction` | — | false | reserved (extended stage) |
| | `tangential_stiffness_ratio` | — | 1.0 | tangential stiffness scale |
| `settle` | `gravity_boost` | — | 1.0 | multiplies g during settling only (compaction) |
| | `ke_tolerance` | J/particle | 1e-8 | stop criterion on kinetic energy |
| | `min_steps` | steps | 2000 | minimum settling steps |
| | `max_steps` | steps | 0 = no cap | hard cap on settling steps |
| | `velocity_damping` | 1/s | 0 | decay `v *= exp(-λ·dt)`; 0 = pure physics |
| `gravity` | `gz` | m/s² | −9.81 | −1.62 for the lunar analogue |
| `boundaries` | `periodic_xy` | — | false | periodic side walls |
| | `wall_restitution` | — | 0.30 | restitution on walls and floor |
| | `fixed_bottom` | — | true | rigid floor at z = 0 |
| `output` | `dir` | — | `data/outputs` | output folder (created if missing) |
| | `vtk` | — | true | write `.vtu` + `.pvd` |
| | `csv` | — | true | write the series CSV |
| | `vtk_every` | steps | 2000 | write a `.vtu` every N steps |

### Two method notes on the values

- **`young = 1e8 Pa` in all six cases.** That is a softened stiffness compared with real
  quartz (≈ 70 GPa), the standard DEM practice to keep the time step tractable. At
  70 GPa the Rayleigh time drops by a factor ≈ 26, and with it the admissible `dt`.
- **`packing_fraction` is a target, not data.** The generator seeds the nested lattice
  and expands it to land on that density; the real bed value is measured and reported
  (`seeded phi`, `settled phi` in the log; `seeded_phi` and `settled_phi` at the end of
  the CSV).

## 2. The six cases

| Config | v (m/s) | angle | target φ | g (m/s²) | Bed grains | Purpose |
|---|---|---|---|---|---|---|
| `default.json` | 2.0 | 90° | 0.524 | 9.81 | 21 904 | quick case: validates the full pipeline |
| `loose_vertical.json` | 2.0 | 90° | 0.524 | 9.81 | ≈ 50k | thesis replica, loose bed |
| `compact_vertical.json` | 2.0 | 90° | 0.574 | 9.81 | ≈ 60k | compacted bed (`gravity_boost = 3`) |
| `loose_oblique45.json` | 2.8284 | 45° | 0.524 | 9.81 | ≈ 50k | circularity paradox (same `v_n`) |
| `lunar_impact.json` | 5.0 | 75° | 0.524 | 1.62 | ≈ 50k | lunar analogue |
| `energy_sweep.json` | 1–6 m/s | 90° | 0.524 | 9.81 | ≈ 50k | D(E) law (driven by `scripts/energy_sweep.py`) |

Grain size is what sets the particle count: all cases use the same 45 × 45 × 25 cm
sandbox with the bed at 60 % height. `default` uses 5.0–5.5 mm radii (21 904 grains,
minutes per case) and the others 3.7–4.2 mm (≈ 50k grains), the step that motivates the
OpenMP stage. The energy sweep has an extra block that only the script reads:

```json
"sweep": { "parameter": "projectile.speed", "values": [1.0, 1.5, 2.0, 2.5, 3.0, 4.0] }
```

## 3. Command line

```
granimpact <config.json> [options]

  --backend <serial|openmp|cuda>   contact solver (default: serial)
  --out <dir>                      output folder
  --name <name>                    case name (file prefix)
  --dt <s>                         integration step
  --t-end <s>                      impact duration
  --t-settle <s>                   settling duration
  --seed <int>                     seed
  --speed <m/s>                    projectile speed (sweeps)
  --no-vtk                         no VTK files (much faster)
  --backends                       list the backends available in this build
  -h, --help
```

Examples:

```bash
granimpact configs/default.json --out data/outputs      # quick case
granimpact configs/lunar_impact.json --out data/outputs --no-vtk
granimpact configs/energy_sweep.json --speed 4.0 --name sweep_v4
```

## 4. Outputs

| File | Content |
|---|---|
| `<case>_crater_evolution.csv` | time series + final block with the measured metrics |
| `<case>_summary.txt` | human-readable summary: parameters, result, timings |
| `<case>_<step>.vtu`, `<case>.pvd` | per-particle fields; open the `.pvd` in ParaView |

CSV columns (all SI unless noted):

| Column | Meaning |
|---|---|
| `t` | simulated time |
| `particles`, `contacts` | number of particles and of active contacts |
| `kinetic_energy`, `potential_energy` | kinetic and gravitational potential energy [J] |
| `max_speed` | maximum speed [m/s] |
| `d_exc` | maximum crater excavation [m] |
| `D` | major diameter of the cavity (ellipse by moments) [m] |
| `epsilon` | ellipse eccentricity (0 = circular) |
| `h_rim` | rim height [m] |
| `aspect_zd` | depth-to-diameter ratio |
| `v_in` | cavity volume [m³] |
| `ejecta_mass` | mass flagged as ejecta [kg] |
| `classification` | morphological class (`Simple`, `Simple (deep)`, `Sand Mound / no crater`, …) |

The **last line** of the CSV is a `# result:` block with the impact speed measured at
first contact, the settling steps, the backend, the measured densities and the contrast
model values. It is what lets a results folder be analyzed without the configs
(`granimpact_post` uses it).

## 5. Adding a case

1. Copy an existing config and change as little as possible (one case = one physics question).
2. `granimpact <config>.json --no-vtk --t-end 0.005` to check that it starts.
3. Add it to `scripts/run_all.sh` if it belongs to the full batch.
4. Document it in the table above: which question it answers and what makes it different.
