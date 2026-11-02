# Configuration, Reference Cases & Output Specification

## 1. JSON Configuration Schema (`SimConfig`)

`SimConfig::fromFile()` + `SimConfig::validate()` parse and validate every parameter in SI
units (meters, seconds, kilograms, pascals) unless explicitly suffixed (`_kPa`, `_deg`).
Any out-of-range value aborts immediately at startup with a descriptive diagnostic.

| JSON Block | Key | Unit | Default | Physical / Numerical Meaning |
| :--- | :--- | :---: | :---: | :--- |
| `simulation` | `name` | — | `default` | Case identifier and output file prefix |
| | `dt` | `s` | `1e-6` | Integration timestep (checked against Rayleigh limit) |
| | `t_settle` | `s` | `0.02` | Maximum bed settling & compaction duration |
| | `t_end` | `s` | `0.05` | Duration of the impact phase after projectile release |
| | `output_every` | `steps` | `500` | Interval for CSV time-series and radial/subsurface snapshots |
| | `seed` | — | `20240517` | Master `std::mt19937_64` seed (strict determinism) |
| `domain` | `x`, `y`, `z` | `m` | `0.45, 0.45, 0.25` | Experimental sandbox dimensions (`450 × 450 × 250 mm`) |
| `particles` | `r_min`, `r_max` | `m` | `4e-4, 6e-4` | Polydisperse bed grain radius bounds (`r_max/r_min = 1.35`) |
| | `rho_grain` | `kg/m³` | `2650.0` | Intrinsic quartz grain material density |
| | `packing_fraction` | — | `0.524` | Target bed packing fraction (`0.524` loose, `0.574` compact) |
| | `bed_height_fraction` | — | `0.60` | Fraction of `domain.z` filled during initial seeding |
| | `young` | `Pa` | `70e9` | Grain Young's modulus (`1e8 Pa` softened DEM in configs) |
| | `poisson` | — | `0.25` | Grain Poisson's ratio |
| | `restitution` | — | `0.55` | Normal grain–grain coefficient of restitution `e` |
| | `friction` | — | `0.40` | Sliding Coulomb friction coefficient `μ` (`0.50` compact) |
| `projectile` | `radius` | `m` | `0.0354` | Spherical projectile radius (`D_p = 70.7 mm`) |
| | `density` | `kg/m³` | `1310.0` | Bulk density of consolidated projectile (`M_p = 243.46 g`) |
| | `phi` | — | `0.50` | Projectile internal solid packing fraction `φ_p` |
| | `young`, `poisson` | `Pa`, — | `5e6, 0.20` | Projectile elastic modulus and Poisson's ratio |
| | `yield_kPa` | `kPa` | `17.15` | Uniaxial compressive yield strength `σᵧ` (Thesis / PRE 2024) |
| | `speed` | `m/s` | `2.0` | Initial release velocity magnitude `v₀` |
| | `angle_deg` | `deg` | `90.0` | Impact angle above horizontal (`90°` vertical, `45°` oblique) |
| | `start_z` | `m` | `0.20` | Initial release height of projectile center |
| `contact` | `tangential` | — | `true` | Enables Mindlin–Deresiewicz history + Coulomb cap |
| | `rolling_friction` | — | `false` | Enables Ai et al. (2011) rotational rolling resistance |
| | `rolling_friction_coeff` | — | `0.15` | Dimensionless rolling resistance coefficient `μ_r` |
| | `tangential_stiffness_ratio` | — | `1.0` | Multiplier on Mindlin tangential stiffness `k_t` |
| | `cohesive_projectile` | — | `true` | Enables breakable cement bonds inside projectile |
| | `bed_cohesion_kPa` | `kPa` | `0.0` | Capillary / moisture / van der Waals cohesion in bed |
| | `magnetic_susceptibility` | — | `0.0` | Grain magnetic susceptibility `χ_m` for magnetic sand |
| | `magnetic_field_T` | `T` | `0.0` | Vertical magnetic field `B₀` for dipole–dipole forces |
| `settle` | `gravity_boost` | — | `1.0` | Settling gravity multiplier (`1.0` loose, `3.0` compact) |
| | `compaction_pressure_kPa` | `kPa` | `0.0` | Uniaxial top-lid press during settling (`4.46 kPa` compact) |
| | `ke_tolerance` | `J/grain` | `1e-8` | Kinetic energy per particle threshold for settling convergence |
| | `min_steps` | `steps` | `2000` | Minimum settling steps before checking convergence |
| | `max_steps` | `steps` | `0` | Hard cap on settling steps (`0` = use `t_settle / dt`) |
| | `velocity_damping` | `1/s` | `0.0` | Settling-only viscous drag (`0` during impact phase) |
| `gravity` | `gz` | `m/s²` | `-9.81` | Vertical gravity (`-9.81` Earth, `-1.62` Moon) |
| `boundaries` | `periodic_xy` | — | `false` | Periodic lateral boundaries in `X` and `Y` |
| | `wall_restitution` | — | `0.30` | Restitution coefficient against container walls/floor |
| | `wall_friction` | — | `0.40` | Coulomb friction against container walls/floor |
| | `wall_damping` | — | `0.15` | Viscoelastic damping factor on wall collisions |
| | `fixed_bottom` | — | `true` | Rigid floor at `z = 0` |
| `output` | `dir` | — | `data/outputs/loose/default` | Structured output path (`data/outputs/{loose,compact}/<case>`) |
| | `vtk` | — | `true` | Write `.vtu` snapshots and `.pvd` ParaView collection |
| | `csv` | — | `true` | Write time-series, depth map, radial, and subsurface CSVs |
| | `vtk_every` | `steps` | `2000` | Timestep interval between `.vtu` snapshots |

---

## 2. Canonical Reference Cases (`configs/*.json`)

All configurations model the `450 × 450 × 250 mm` sandbox with a `D_p = 70.7 mm` (`M_p = 243.46 g`,
`σᵧ = 17.15 kPa`) granular projectile. Outputs are automatically organized into
`data/outputs/loose/<case>/` or `data/outputs/compact/<case>/`:

| Config File | Radii (mm) | Bed Grains | v₀ (m/s) | Angle | Target φ | Lid Press | Gravity (m/s²) | Output Directory |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| `default.json` | `4.1–5.5` | `25,857` | `2.0000` | `90°` | `0.524` | `0.00 kPa` | `-9.81` | `data/outputs/loose/default` |
| `loose_vertical.json` | `3.1–4.2` | `58,926` | `2.0000` | `90°` | `0.524` | `0.00 kPa` | `-9.81` | `data/outputs/loose/loose_vertical` |
| `compact_vertical.json` | `3.1–4.2` | `66,735` | `2.0000` | `90°` | `0.574` | `4.46 kPa` | `-9.81` | `data/outputs/compact/compact_vertical` |
| `loose_oblique45.json` | `3.1–4.2` | `58,926` | `2.8284` | `45°` | `0.524` | `0.00 kPa` | `-9.81` | `data/outputs/loose/loose_oblique45` |
| `compact_oblique45.json` | `3.1–4.2` | `66,735` | `2.8284` | `45°` | `0.574` | `4.46 kPa` | `-9.81` | `data/outputs/compact/compact_oblique45` |
| `lunar_impact.json` | `3.1–4.2` | `58,926` | `5.0000` | `75°` | `0.524` | `0.00 kPa` | `-1.62` | `data/outputs/loose/lunar_impact` |
| `lunar_oblique45.json` | `3.1–4.2` | `58,926` | `2.8284` | `45°` | `0.524` | `0.00 kPa` | `-1.62` | `data/outputs/loose/lunar_oblique45` |
| `energy_sweep.json` | `3.1–4.2` | `58,926` | `1.0–20.0` | `90°` | `0.524` | `0.00 kPa` | `-9.81` | `data/outputs/loose/energy_sweep` |

---

## 3. Command-Line Interface (`granimpact_serial`)

Every physical, geometric, and numerical parameter can be overridden from the command line
without editing JSON files:

```text
Usage:
  ./build/release/bin/granimpact_serial [config.json] [options]

Core & Execution Options:
  -c, --config <path>      path to JSON configuration file
  --backend <name>         serial | openmp | cuda  (default: serial)
  -o, --out, --output <d>  output folder (overrides config.output.dir)
  -n, --name <name>        case name and file prefix
  --dt <seconds>           integration timestep [s]
  --t-settle <seconds>     bed settling budget [s]
  --t-end <seconds>        impact duration [s]
  --seed <uint64>          random generator seed
  --no-vtk                 disable VTK .vtu/.pvd output for faster batch sweeps
  --backends               list compiled contact solvers and exit
  -h, --help               show help message

Kinematics, Bed & Projectile Overrides:
  -v, --speed <m/s>        initial projectile speed v_0 [m/s]
  -a, --angle <deg>        impact angle in (0, 90] deg (90 = vertical)
  -p, --phi <0..0.64>      target bed packing fraction phi_0
  -g, --gravity <m/s2>     vertical gravity gz (-9.81 Earth, -1.62 Moon)
  --box-x <m>              container width X [m]
  --box-y <m>              container length Y [m]
  --box-z <m>              container height Z [m]
  --bed-height-frac <0..1> fraction of Z filled by the initial bed
  --r-min <m>              minimum grain radius [m] (scales particle count N)
  --r-max <m>              maximum grain radius [m] (scales particle count N)
  --proj-radius <m>        projectile sphere radius R_p [m]
  --proj-phi <0..0.64>     projectile internal packing fraction phi_p
  --proj-density <kg/m3>   projectile bulk density rho_p [kg/m3]
  --proj-yield <kPa>       projectile cohesive yield strength sigma_y [kPa]
  --compaction-kpa <kPa>   uniaxial top-lid compaction pressure [kPa]
  --gravity-boost <factor> settling gravity multiplier
  --bed-cohesion <kPa>     capillary / moisture bed cohesion [kPa]
  --mag-chi <val>          magnetic susceptibility chi_m for magnetic sand
  --mag-b <Tesla>          ambient vertical magnetic field B_0 [Tesla]
```

### Practical CLI Examples

```bash
# 1. Run loose vertical impact and generate all 3 plots + craterslab .npz
./build/release/bin/granimpact_serial configs/loose_vertical.json
python3 -m granimpact_post visualize-all \
    data/outputs/loose/loose_vertical/loose_vertical_crater_evolution.csv

# 2. Grazing 20-deg tangential lunar impact via CLI override
./build/release/bin/granimpact_serial configs/lunar_oblique45.json \
    --angle 20 --name lunar_tangential20 \
    --out data/outputs/loose/lunar_tangential20

# 3. Cohesive wet sand or magnetic regolith experiment
./build/release/bin/granimpact_serial configs/default.json \
    --bed-cohesion 2.5 --mag-chi 0.5 --mag-b 0.05 \
    --name cohesive_magnetic_test --out data/outputs/loose/cohesive_magnetic_test
```

---

## 4. Output Files, `craterslab` `.npz` Contract & ParaView Guide

Each simulation run populates `data/outputs/{loose,compact}/<case>/` with the following
self-contained files:

| File Pattern | Generated By | Description |
| :--- | :--- | :--- |
| `<case>_crater_evolution.csv` | C++ (`CsvWriter`) | Time series (`26` columns) + `# result:` metadata footer |
| `<case>_depth_map.csv` | C++ (`writeDepthMapCsv`) | `96 × 96` Virtual ToF grid (`x_mm`, `y_mm`, `z_bed_mm`, `z_total_mm`) |
| `<case>_subsurface_chains.csv` | C++ (`appendSubsurfaceFieldCsv`) | Cylindrical `(r, z)` snapshots of `phi_local`, `stress_kPa`, `coordination` |
| `<case>_radial_evolution.csv` | C++ (`appendRadialProfileCsv`) | Azimuthally averaged radial profiles `z_bed_mm(r, t)` & `z_total_mm(r, t)` |
| `<case>_summary.txt` | C++ (`writeCaseSummary`) | Human-readable summary of parameters, observables, models, and wall time |
| `<case>_<step>.vtu`, `<case>.pvd` | C++ (`VtkWriter`) | 3D unstructured grid snapshots + ParaView time collection |
| `<case>_craterslab.npz` | Python (`granimpact_post`) | Direct `craterslab` elevation tensor (`xp`, `yp`, `M`, `z_bed_mm`, `z_total_mm`) |
| `<case>_force_chains.png` | Python (`granimpact_post`) | Visualization 1: 2D `(r, z)` virial stress chains & Heckel RCP saturation |
| `<case>_lateral_saturation.png` | Python (`granimpact_post`) | Visualization 2: Momentum redirection `κ_vh`, rim slope `tan θ` & `h_rim` saturation |
| `<case>_crater_formation.png` | Python (`granimpact_post`) | Visualization 3: Radial crater evolution, porous piston & 2D ToF depth map |

### 4.1. Columns of `<case>_crater_evolution.csv`

| Column(s) | Units | Physical Meaning |
| :--- | :---: | :--- |
| `t`, `particles`, `contacts` | `s`, —, — | Simulation time, active grain count, and active contact count |
| `kinetic_energy`, `potential_energy`, `rotational_energy` | `J` | Translational, gravitational, and 6-DOF rotational energies |
| `dissipated_energy` | `J` | Cumulative frictional + viscoelastic work dissipated during impact (`W_fric`) |
| `max_speed` | `m/s` | Maximum grain speed in the container |
| `d_exc`, `D`, `epsilon` | `m`, `m`, — | Excavated crater depth, major diameter, and cavity eccentricity |
| `h_rim`, `rim_slope`, `aspect_zd` | `m`, —, — | Raised rim height, inner crater wall slope `tan θ`, and aspect ratio `d_exc / D` |
| `v_in`, `v_ex`, `v_exc` | `m³` | Excavated cavity volume, raised rim volume, and compressed excess volume |
| `h_cp` | `m` | Central peak / porous piston height above cavity floor |
| `phi_subsurface` | — | Packing fraction measured in the sub-crater compression column (`φ_sub`) |
| `mean_virial_stress_kPa`, `max_virial_stress_kPa` | `kPa` | Mean and peak Love–Weber virial compressive stress beneath the crater |
| `momentum_transfer_vh` | — | Horizontal-to-vertical momentum redirection ratio `κ_vh` |
| `active_bonds`, `broken_bonds` | — | Intact and fractured cohesive cement bonds in the projectile |
| `ejecta_mass`, `classification` | `kg`, — | Escaped ejecta mass and automated class (`Simple`, `Complex (central peak)`, etc.) |

### 4.2. Visualizing 3D `.pvd` / `.vtu` Files in ParaView

1. **Install ParaView on Ubuntu/Debian**:
   ```bash
   sudo apt update && sudo apt install -y paraview
   ```
2. **Open the simulation collection**:
   ```bash
   paraview data/outputs/loose/default/default.pvd &
   ```
   Click **Apply** in the left *Properties* panel.
3. **Render true-scale spheres (`Point Gaussian`)**:
   - Change representation from *Surface* to **Point Gaussian**.
   - Enable **Use Scale Array**, set **Gaussian Scale Array** to `radius`, and set
     **Gaussian Radius** to `1.0`.
4. **Color by physical per-grain fields**:
   - `type`: Distinguishes bed sand (`0`), projectile (`1`), and escaped ejecta (`2`).
   - `virial_stress_kPa`: Visualizes 3D Love–Weber compressive force chains `[kPa]`.
   - `contact_force`: Total normal contact force `Σ|F_n|` acting on each grain `[N]`.
   - `coordination`: Active contact count `Z_i` per grain.
   - `speed` and `omega`: Translational speed `[m/s]` and angular speed `[rad/s]`.
5. **Cross-sectional cut (`Clip`)**:
   - Add a **Clip** filter with **Plane Normal** `(0, 1, 0)` through the box center
     `y = 0.225 m` to inspect the subsurface compacted plug and comminuted projectile.
