# Physical Model & Governing Equations

Everything described in this document is implemented in `src/physics/`, `src/generation/`,
`src/sim/`, and `src/analysis/`. The mathematical formulation reproduces the low-velocity
granular-granular impact experiments of Corrales-Machin et al., *Physical Review E* **109**,
054907 (2024), *Computational Particle Mechanics* (2023), and the PhD Thesis (DICIM 2024).

---

## 1. 6-DOF Equations of Motion (`Integrator.cpp`)

Each spherical grain `i` has radius `r_i`, mass `m_i = (4/3) π ρ_i r_i³`, solid-sphere
moment of inertia `I_i = (2/5) m_i r_i²`, position `x_i`, translational velocity `v_i`,
and angular velocity `ω_i`. The system is integrated via the symplectic **6-DOF Velocity
Verlet** (*kick–drift–kick*) scheme:

```text
1. First half-kick (stepFirstHalf):
   v_i(t + Δt/2) = v_i(t) + (F_i(t) / m_i) · (Δt/2)
   ω_i(t + Δt/2) = ω_i(t) + (τ_i(t) / I_i) · (Δt/2)
   x_i(t + Δt)   = x_i(t) + v_i(t + Δt/2) · Δt

2. Force & torque evaluation (ContactModel::computeForces):
   Recompute F_i(t + Δt) and τ_i(t + Δt) from all active contacts and external fields.

3. Second half-kick (stepSecondHalf):
   v_i(t + Δt)   = v_i(t + Δt/2) + (F_i(t + Δt) / m_i) · (Δt/2)
   ω_i(t + Δt)   = ω_i(t + Δt/2) + (τ_i(t + Δt) / I_i) · (Δt/2)
```

Symplectic energy conservation and exact quadratic trajectory integration under constant
acceleration are verified in `tests/test_integrator.cpp`.

---

## 2. 6-DOF Grain–Grain Contact Mechanics (`ContactModel.cpp`)

For two touching grains `i` and `j` with overlap `δ = (r_i + r_j) - ||x_j - x_i|| > 0`
and unit contact normal `n = (x_j - x_i) / ||x_j - x_i||`, the effective pair parameters
are:

```text
1/E* = (1 - ν_i²)/E_i + (1 - ν_j²)/E_j,     1/G* = (2 - ν_i)/G_i + (2 - ν_j)/G_j
R*   = (r_i r_j)/(r_i + r_j),               m_eff = (m_i m_j)/(m_i + m_j)
```

### 2.1. Hertz–Tsuji Normal Viscoelastic Force
```text
k_n = (4/3) E* √R*,     F_n,el = k_n δ^(3/2)
γ_n = -2 ln(e) √(m_eff · 2 E* √(R* δ)) / √(π² + ln²(e))
F_n = max(0, F_n,el - γ_n v_n)
```
where `v_n = (v_j - v_i) · n` is the normal relative velocity and `e` is the restitution
coefficient.

### 2.2. Mindlin–Deresiewicz Tangential History & 6-DOF Slip Velocity
Unlike 3-DOF point-mass approximations, the relative velocity at the contact point
includes surface rotation of both spheres:

```text
v_c = (v_j - v_i) - v_n n - (r_i ω_i + r_j ω_j) × n
```

The elastic tangential displacement `δ_t` is stored in the flat 16-slot per-particle SoA
table (`ContactModel::history_slots_`), rotated into the current tangent plane at every
timestep, and incremented by `v_c Δt`:

```text
k_t = 8 G* √(R* δ) · ratio_t
F_t = -k_t δ_t,     capped by Coulomb sliding: ||F_t|| ≤ μ |F_n|
```
When Coulomb slip occurs (`||F_t|| > μ |F_n|`), `δ_t` is rescaled to the slip cone and
the incremental frictional work `dW_fric = ||F_t|| · ||v_c|| Δt` is accumulated into
`dissipated_energy_`.

### 2.3. Rotational Torque & Ai et al. (2011) Rolling Resistance
Tangential contact forces exert opposite lever-arm torques on both grains, complemented
by directional rolling friction (`rolling_friction_coeff` `μ_r = 0.15–0.18`) that models
grain angularity and sustains experimental angles of repose (`40°–44°`):

```text
ω_rel = ω_i - ω_j
τ_roll = -μ_r R* |F_n| (ω_rel / ||ω_rel||)
τ_i   -= r_i (n × F_t) - τ_roll
τ_j   -= r_j (n × F_t) + τ_roll
```

### 2.4. Weakly Consolidated Projectile Bonds & Comminution (`σᵧ = 17.15 kPa`)
In the experiments (Corrales-Machin et al., *PRE* 2024), the projectile (`D_p = 70.7 mm`,
`M_p = 243.46 g`, `ρ_p = 1310 kg/m³`, `φ_p = 0.50`) is cemented by a weak gypsum/Portland
matrix with uniaxial compressive yield strength `σᵧ = 17.15 kPa`. At projectile creation,
`ContactModel::registerCohesiveBonds()` registers breakable cement bridges across a neck
area `A_bond = π (min(r_i, r_j))²`:
- While intact (`cohesive = 1`), the bond resists tension up to `F_pull = σᵧ A_bond` and
  shear up to `F_shear = c₀ A_bond + μ |F_n|`.
- When impact chains exceed `σᵧ` (`P₀ = ½ ρ_p v₀² ≥ σᵧ`, **Regime 2**), the bond fractures
  irreversibly (`cohesive = 0`, incrementing `broken_bonds_`), transforming the projectile
  into a comminuted **porous piston** that deposits a central peak inside the crater.

### 2.5. Capillary/Moisture Bed Cohesion & Magnetic Dipole Interactions
For cohesive planetary regoliths or wet/magnetic sands, `ContactModel` supports two
additional physical interactions controlled via JSON or CLI flags:
1. **Capillary / van der Waals Bed Cohesion (`bed_cohesion_kPa`)**: Adds a short-range
   tensile pull-off resistance `F_coh = σ_coh · π (R*)²` between contacting bed grains,
   shifting the Coulomb yield envelope to `||F_t|| ≤ μ (|F_n| + F_coh)`.
2. **Magnetic Sand Dipole–Dipole Attraction (`magnetic_susceptibility`, `magnetic_field_T`)**:
   Under an external vertical magnetic field `B = (0, 0, B₀)`, grains with susceptibility
   `χ_m` acquire an induced magnetic moment `m_i = (4/3) π r_i³ (χ_m / μ₀) B` and exert
   anisotropic dipole–dipole forces (`F_mag ∝ (1 - 3 cos²θ_B) / d⁴`), enhancing vertical
   force chains and column stiffness.

---

## 3. Neighbor Search & Flat SoA Contact History (`CellList.cpp`)

`CellList` partitions the container into a Compressed Sparse Row (CSR) grid with cell side
`h_cell = 2 · r_max · 1.05`:
- Internal vectors (`head_`, `next_`, `cell_of_`) are reused across steps with **zero
  dynamic heap allocations**.
- Each cell visits only `14` of the `27` neighbor offsets (self + 13 forward cells) with
  canonical ordering `i < j`, visiting each candidate pair exactly once in `O(N)` time.
- Tangential history and cohesive flags are stored on the lower index `i < j` inside a
  contiguous vector `history_slots_` of size `N × 16` (`kMaxContactsPerParticle = 16`).
  Stale contacts are recycled automatically via `step_stamp_` without hash-table overhead.

---

## 4. Boundaries, Settling & Uniaxial Top-Lid Compaction

### 4.1. Container Boundaries (`Integrator.cpp`)
| Boundary | Physical Treatment |
| :--- | :--- |
| **Floor (`z = 0`)** | Rigid frictional boundary with restitution `wall_restitution`, Coulomb friction `wall_friction`, and normal damping `wall_damping` |
| **Side walls (`x, y`)** | Rigid frictional walls (`450 × 450 mm` container) or periodic if `periodic_xy = true` |
| **Ceiling (`z > L_z`)** | Grains escaping above the container are flagged as `ParticleKind::Ejecta` so they cannot fall back into the cavity |

### 4.2. Loose vs. Pre-Compacted Bed Preparation (`Simulation.cpp`)
To replicate the two experimental bed regimes without altering impact-phase physics:
1. **Loose Bed (`packing_fraction = 0.524`, `compaction_pressure_kPa = 0.0`)**: Grains
   settle under normal gravity (`gz = -9.81 m/s²`), converging to `φ₀ = 0.527` (`H ≈ 109 mm`).
2. **Pre-Compacted Bed (`packing_fraction = 0.574`, `compaction_pressure_kPa = 4.46`)**:
   Replicates the experimental `92 kg` metallic press (`P_lid = 92 × 9.81 / 0.45² = 4.46 kPa`)
   applied to the top grain layer during the first `75%` of settling, followed by elastic
   relaxation during the final `25%` before recording the reference surface.

---

## 5. Polydisperse Generation, Ballistic Lead & Particle-Count Scaling

### 5.1. Polydisperse Expanded Nested Lattice (`Generators.cpp`)
Grains are seeded on a non-overlapping staggered lattice (`s = 2.01 r_max`, `dz = 1.42 r_max`)
with radii drawn uniformly in `[r_min, r_max]` (`r_max / r_min = 1.35` breaks hexagonal
crystallization), isotropically expanded to match the target packing fraction, and
perturbed with an `8%` random jitter verified against a spatial hash grid.

### 5.2. Exact Projectile Mass Calibration & Oblique Ballistic Lead
- **Mass calibration**: Because a discrete sphere of radius `R_p` holds `N_p` grains of
  varying radii, `Generators::generateProjectile()` calls `scaleMassOf()` on the projectile
  grains so the total projectile mass equals the exact experimental value
  `M_p = ρ_p (4/3) π R_p³ = 243.46 g`.
- **Oblique & Lunar Ballistic Lead Targeting**: For oblique impacts (`angle_deg < 90°`)
  and low-gravity lunar runs (`gz = -1.62 m/s²`), `generateProjectile()` solves the exact
  ballistic quadratic equation from release height `z₀` to the settled surface `z_bed`:
  ```text
  t_hit = (v_z0 + √(v_z0² + 2 |g_z| Δz)) / |g_z|,     x_0 = x_target - v_x0 · t_hit
  ```
  This guarantees that the projectile always strikes the exact center of the container
  `(L_x/2, L_y/2, z_bed)` regardless of impact angle or planetary gravity.

### 5.3. Scaling Particle Count (`N`) vs. Packing Fraction (`φ`)
The solid volume of a polydisperse grain with uniform radius distribution in `[r_min, r_max]`
has mean volume:
```text
<V_g> = (4/3) π <r³>,     where <r³> = (r_max⁴ - r_min⁴) / (4 (r_max - r_min))
```
Therefore, for a bed of volume `V_bed = L_x L_y H_bed` and a spherical projectile of
volume `V_p = (4/3) π R_p³`:
```text
N_bed ≈ (φ_0 · L_x L_y H_bed) / <V_g>,     N_p ≈ φ_p · (R_p³ / <r³>)
```
- **To keep packing fraction `φ` constant while changing dimensions (`R_p` or `L_x, L_y, H_bed`)**:
  Keep `r_min, r_max` fixed (`N` scales proportionally to container/projectile volume), or
  scale grain radii proportionally (`r' = α r`) if you want to keep `N` constant.
- **To vary packing fraction `φ` at fixed dimensions**: Change `--phi` (or `packing_fraction`)
  and `--compaction-kpa` for the bed, or `--proj-phi` and `--proj-density` for the projectile.

---

## 6. Rayleigh Time-Step Stability

The critical Rayleigh wave time step for Hertzian spheres is computed at startup:
```text
t_R = π r_min √(ρ_grain / G) / (0.1631 ν + 0.8766),     Δt ≤ t_R / 20
```
With softened DEM Young's modulus `E = 10⁸ Pa` (`ν = 0.25`, `ρ = 2650 kg/m³`), `t_R ≈ 98 μs`
for `r_min = 3.1 mm`, so `Δt = 2.0 μs` (`< t_R / 45`) satisfies Rayleigh stability with
ample margin even at `v₀ = 20 m/s`.

---

## 7. Virtual ToF Camera, Observables & Archimedes Subsurface Fields

### 7.1. 3-Pass Virtual Time-of-Flight (ToF) Surface Reconstruction (`CraterAnalyzer.cpp`)
To emulate the experimental overhead optical ToF sensor without inter-grain ray leakage:
1. **Pass 1 (Center Projection)**: Records the top surface `z_i + r_i` of the highest grain
   in each cell of the `96 × 96` grid (`Δx ≈ 4.69 mm`).
2. **Pass 2 (Sphere-Cap Optical Occlusion)**: Projects each upper-layer sphere onto neighboring
   cells within `1.15 r_i` using the upper-hemisphere chord height, eliminating artificial
   inter-grain gaps.
3. **Pass 3 (`3 × 3` Pinhole Pore Filter)**: Fills isolated sub-grain pores using the local
   75th percentile of valid neighbor cells and masks a `3`-cell perimeter border to
   eliminate rigid-wall meniscus artifacts.

Combined with pure bed-grain solid volume accounting (`solidsVolumeOf(ParticleKind::Bed)`),
the settled loose bed measures `φ₀ = 0.527` (matching the `0.524` target within `0.6%`).

### 7.2. Exact Archimedes Spherical-Slice Subsurface Fields (`subsurfaceField`)
Point-center binning of coarse grains creates severe horizontal stripes in `(r, z)` packing
maps. `CraterAnalyzer::subsurfaceField()` eliminates discretization artifacts by integrating
the **exact analytical volume of intersection** of every sphere `[z_i - r_i, z_i + r_i]`
with each vertical slab `[z_low, z_high]` via Archimedes' spherical slice formula:
```text
V_slice(a, b) = π [ r_i² (b - a) - (b³ - a³)/3 ]
```
and evaluating the Love–Weber virial compressive stress per grain:
```text
σ_virial,i = (r_i · Σ_j |F_n,ij|) / (3 · V_grain,i)
```

---

## 8. Analytical Contrast Models & Rim Height Saturation Theory

`CraterAnalyzer::models` evaluates five physical reference models alongside every run:

| Physical Model | Mathematical Expression | Physical Regime & Role |
| :--- | :--- | :--- |
| **Uehara (2003)** | `D = 0.92 (ρ_bed / (ρ_grain μ²))^(1/4) D_p^(3/4) H^(1/4)` | Solid-sphere excavation scaling in cohesionless sand |
| **Logarithmic (Thesis Eq. 3.10)** | `D = 0.089 · log₁₀(E)` | Low-velocity granular-granular diameter scaling (`E > 1 J`) |
| **Uniaxial Heckel Packing (PRE 2024)** | `φ(P) = φ_max - (φ_max - φ₀) exp(-K max(0, P - P₀))` | Subsurface column compaction saturating at `φ_max = 0.64` (RCP) |
| **Heckel Excavation Depth** | `d_Heckel = h_col (1 - φ₀ / φ(P))` | Permanent crater depth component due to uniaxial bed compression |
| **Double-Log Rim Saturation** | `h_rim(E) = h_scale · ln(1 + ln(1 + E / E₀))` | Asymptotic saturation of rim height (`h_rim ≈ 18–20 mm`) |

### Why Crater Rim Height (`h_rim`) Saturates at High Impact Energy
Simulation diagnostics (`<case>_lateral_saturation.png` and `<case>_force_chains.png`)
prove that rim saturation arises from three coupled granular mechanisms:
1. **Subsurface RCP Jamming (`φ_sub → 0.64`)**: Under high stagnation pressure (`P₀ ≥ σᵧ`),
   the column beneath the projectile rapidly compresses to the Random Close Packing limit
   (`φ_RCP = 0.64`), acting as a rigid anvil that halts vertical penetration.
2. **Vertical-to-Horizontal Momentum Redirection (`κ_vh`)**: Once the subsurface column
   jams, projectile comminution redirects vertical momentum radially outward (`κ_vh ≈ 1.5`),
   widening the crater diameter `D` and rim volume `V_ex` horizontally rather than raising `h_rim`.
3. **Angle-of-Repose Avalanche & Frictional Sink (`W_fric`)**: As the rim piles up, its
   inner wall slope `tan θ` exceeds the dynamic stability limit fixed by Coulomb (`μ`) and
   rolling friction (`μ_r`), triggering surface avalanches back toward the angle of repose
   while rotational and sliding friction dissipate the remaining kinetic energy (`W_fric`).
