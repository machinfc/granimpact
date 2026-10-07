# Physical model

Everything below is implemented in `src/physics/`, `src/generation/` and
`src/analysis/`. The references are the project's; every number that appears is either
an input parameter of the cases or a measurement from a real run.

## 1. Equations of motion

**Velocity Verlet** integration in *kick–drift–kick* form (`Integrator.cpp`):

1. `v += a·dt/2`
2. `x += v·dt`
3. recompute contact forces
4. `v += a·dt/2`

Advantages: second order in position, symplectic (energy does not drift without
dissipation) and compatible with forces that depend on position and on contact history
(Mindlin tangential). Displacement under constant forces is exact to order dt², which
`tests/test_integrator.cpp` verifies against the analytical free-fall and projectile
solutions.

## 2. Grain–grain contact law (Hertz–Mindlin)

`ContactModel.cpp`. Let δ be the penetration, R* the equivalent radius, E* and G* the
equivalent moduli:

```
1/E* = (1-ν₁²)/E₁ + (1-ν₂²)/E₂            1/G* = (2-ν₁)/G₁ + (2-ν₂)/G₂
R*   = r₁r₂/(r₁+r₂)                        G    = E/(2(1+ν))

F_n  = (4/3) E* √R* δ^{3/2}   +  γ · v_n            (Hertz + damping)
k_s  = 8 G* √(R* δ)                                  (Mindlin tangential stiffness)
F_t  = -k_s δ_t ,   with  |F_t| ≤ μ |F_n|            (per-pair history)
```

Normal damping uses the **Di Renzo–Di Maio** form with γ = -2 ln e · √(m_eff k) /
√(π² + ln² e), i.e. γ consistent with the restitution coefficient `e` in the config.
This is checked in `tests/test_contact_model.cpp`: the restitution measured in a binary
collision stays within ±35 % of the requested value (the formula is a linear-spring
approximation applied to Hertz).

The force is **repulsive**: the pair receives `F` and `-F` with the normal defined from
body i to body j (`force[i] -= n·f_n + f_t`, `force[j] += ...`). That sign is what makes
two overlapping grains separate; a test pins it down.

## 3. Neighbor search

`CellList.cpp`: CSR grid with cells of side `2·r_max·1.05` (the classic maximum-diameter
rule, with margin). Every pair is enumerated **once** using half the offset list (14 of
the 27 neighbours), with canonical pair `(min, max)`. This halves the checks and avoids
counting the dissipation twice. Cost per step is O(N + number of candidate pairs).

## 4. Boundaries

`Integrator.cpp::applyBoundaries`, controlled by the `boundaries` block:

| Boundary | Behaviour |
|---|---|
| Rigid floor (z = 0) | position corrected; bounce with `wall_restitution` and tangential loss ∝ μ |
| Side walls | rigid with `wall_restitution`, or periodic if `periodic_xy = true` |
| Ceiling | the particle is flagged `Ejecta` (leaves the system, accounted separately) |

Ejecta flagging is what allows measuring ejected mass and mean speed without letting
the grains fall back into the crater.

## 5. Bed and projectile generation

`Generators.cpp`. The bed is **not** seeded by random deposition (RSA): with coarse
particles, RSA leaves a bed so loose that it collapses as soon as the simulation starts,
and the subsidence would be read as a false crater. The method is:

1. **Dense nested lattice**: horizontal spacing `s = 2.01·r_max`, vertical spacing
   `dz = 1.42·r_max`, alternating layers shifted by half a cell. The minimum distance
   between a particle and its nearest neighbours is `√(s²/2 + dz²) ≈ 1.00·s = 2.01·r_max`,
   above `2·r_max`: the lattice does not overlap for any radius in the range, although
   the margin is thin (0.5 %), which is why the jitter is checked (point 3).
2. **Isotropic expansion** by `cbrt(φ_dense / φ_target)`, which brings the block density
   exactly to the requested value. If the lattice is already looser than the target
   (very wide radius range) it is not compressed: it is seeded as is and the resulting
   density is reported — compressing would create overlaps.
3. **Small jitter** (8 % of the spacing) with a hash-grid neighbour check, which reduces
   the radius if needed and discards the site if it falls below `r_min`. This preserves
   the radius distribution `[r_min, r_max]`.
4. **Settling**: gravity (+ `settle.gravity_boost` if the case compacts the bed) until
   the kinetic energy per particle drops below `settle.ke_tolerance` **and** the mean
   bed height stabilizes (change < 10⁻³·r_mean), with a minimum of `min_steps` steps.
   Without that second criterion, the energy can drop while the bed is still sinking,
   and the subsidence contaminates the crater reading.

The bed density is therefore a **measured result**: `seeded phi` (lattice volume) and
`settled phi` (real occupied volume) are both reported. The `packing_fraction` in the
config is the target that sets the block size.

The projectile is a sphere of radius R filled with grains on the same nested lattice,
with initial velocity `(v·cosθ, 0, −v·sinθ)`.

## 6. Time stability

Rayleigh criterion for Hertz spheres:

```
t_R = π · r · √(ρ/G) / (0.1631·ν + 0.8766)          dt ≤ t_R / 20
```

The binary computes it at startup and warns if the config's `dt` violates it. With
E = 1·10⁸ Pa (the value in the configs) and r = 5.5 mm, `dt_max ≈ 7.7·10⁻⁶ s`, so the
`dt = 5·10⁻⁶ s` of the default case fits with margin. Method note: 10⁸ Pa is a
**softened** stiffness compared with real quartz (≈ 70 GPa); it is standard DEM practice
to keep the time step tractable (at 70 GPa, `t_R` drops by a factor ≈ 26 and the cost
per step grows in the same proportion).

## 7. Crater observables

`CraterAnalyzer.cpp`. The crater is measured on a per-cell map of z maxima (grid
96 × 96, 4.7 mm per cell in the default), smoothed over a neighbourhood of ≈ 2 cells to
remove granularity noise:

1. **Drift correction**: the bed keeps settling by a few tenths of a grain during the
   impact phase. The mean deformation (reference − current) is measured on the outer
   ring — cells more than 4·R_projectile from the impact — and **the level** is
   corrected, not the shape: fitting a plane to the ring and evaluating it at the centre
   would extrapolate (the settled surface is not flat) and introduced an error of
   several millimetres.
2. **Excavation**: `delta = (z_ref − z) − drift`. `d_exc` is the maximum over the cells
   connected to the deepest point (8-neighbour connected component, threshold
   `margin_factor · r_max` with `margin_factor = 1.0` by default), which prevents a
   stray grain outside the crater from inflating the reading. That threshold is also the
   method's minimum sensitivity: anything deeper than one grain radius is detected.
3. **Ellipse by moments**: variances and covariance are computed over the cavity cells;
   from them come the diameter `D = 4√λ₁` (two standard deviations per side), the minor
   axis, the ellipticity ε = √(1 − λ₂/λ₁) and the orientation angle.
4. **Rim**: `h_rim` positive envelope around the cavity, bounded by radius and clipped
   to 2·r_max.
5. **Volume**: `V_in = π/3 · a² · d_exc` (cavity cone), in cm³.
6. **Classification**: morphological class from Z/D and ε, plus the *deep* flag when the
   excavation exceeds 0.20·D.

Ejected mass and mean speed, and the measured bed density, are reported too.

## 8. Analytical contrast models

`CraterAnalyzer.cpp::models`, with the formulas as used in the analysis (they are
references, **not** fits to data):

| Model | Expression | Declared validity |
|---|---|---|
| Uehara | `D = 1.84 (ρ_bed/ρ_grain)^{1/4} D_proj^{3/4} E^{1/4}` | excavation scaling in sands |
| Logarithmic (thesis) | `D = 0.089 · log₁₀(E)` | slow impacts, E > 1 J |
| Heckel | `d = h (1 − φ₀/φ_pressure)` | compaction by projectile pressure |

The first two can also be fitted to the measured data with
`python3 -m granimpact_post fit` (power law and logarithmic law) and compared by R²:
the verdict comes from the measurements, not from literature coefficients.

## 9. Known limits of the model

### 9.1 The bed compacts too much (main limitation, measured)

The `default` case leaves the bed at **φ ≈ 0.82 measured** (79.8 mm of height out of the
150 mm seeded) when the config target is 0.524. It is not a reading error: the nested
lattice is expanded by 7 % to land on the target, and with that margin the grains no
longer touch, so under gravity they fall into the gaps of the layer below and the bed
compacts to a nearly crystalline packing (the measured density also includes the elastic
compression of the contacts: with E = 1e8 Pa and an 80 mm column, the overlap at the
bottom is of the order of 0.04 mm).

Direct consequence: such a bed is much harder than the loose sand of the experiment and
the crater comes out small — 26.5 mm of measured diameter against 151 mm of Uehara
reference for that same impact (E = 0.156 J).

What has to change (first task of the next iteration): prepare the bed by **deposition**
(a rain of grains with real friction and dissipation, which is how a loose bed is
prepared in the experiment) instead of relaxing a lattice, and validate that the final
density lands within ±5 % of the target before using the results to compare with the
literature.

- **Spherical grains** and no breakage: angularity and fragmentation are not modelled.
- **No cohesion** (0 % humidity): valid for dry sand, not for real regolith.
- **Rolling friction off** by default (`contact.rolling_friction = false`); reserved for
  the extended stage.
- **Tangential history without pruning**: `ContactModel` stores the tangential
  displacement per pair and does not delete it when the contact opens. In long runs with
  many ephemeral pairs that grows in memory; it has to be pruned before production use
  (noted in `docs/roadmap.md`).
- **Softened stiffness** (see §6): it changes the duration of the contact, not the
  excavation scale in the quasi-static regime; it should be checked with a sensitivity
  case before publishing.
