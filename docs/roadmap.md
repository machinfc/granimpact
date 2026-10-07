# Roadmap by stages

The project advances in stages: each one adds **a compute backend** and raises the
particle count. The skeleton is already prepared for all of them: what remains is adding
code at the marked points, not reorganizing anything (see `architecture.md`).

## 1. Particle ladder

The sandbox is the same at every stage (0.45 × 0.45 × 0.25 m, bed at 60 % height, target
φ 0.524). What changes is the grain size, which is what sets the particle count:

| Stage | Radii | Grains | Lattice (mm) | Backend |
|---|---|---|---|---|
| 3 (current) | 5.0–5.5 mm | **21 904** | 11.80 × 8.33 | serial |
| 3-4 (replica) | 3.7–4.2 mm | **52 822** | 8.88 × 6.27 | serial → OpenMP |
| 3-4 (compact) | 3.7–4.2 mm, φ 0.574 | **59 823** | 8.61 × 6.09 | serial → OpenMP |
| 4 | 1.75–2.00 mm | **540 225** | 4.22 × 2.98 | OpenMP |
| 5 | 1.35–1.55 mm | **1 182 447** | 3.26 × 2.30 | CUDA |
| future | 0.4–0.6 mm (thesis sand) | **≈ 24 078 816** | 1.21 × 0.85 | GPU/MPI, beyond current scope |

The numbers come from the nested-lattice geometry (`docs/physics.md` §5): horizontal
spacing `2.01·r_max` expanded to land on the target density, layers every `1.42·r_max`.
The `default` case is verified against the real run (21 904 ✓). With the thesis sand
(0.4–0.6 mm) the radius range is so wide (1.5:1) that the dense lattice falls below the
target: the generator does not compress (it would create overlaps) and reports the
seeded density, which is the honest thing to do.

## 2. Status per stage

### Stage 3 — serial core ✅ done

Delivered: validated `SimConfig`, SoA `ParticleSystem`, CSR cell lists, Hertz–Mindlin
contact with tangential history, Velocity Verlet, boundaries, bed and projectile
generation, morphological analysis, CSV + VTK + summary, `install.sh`, 7 test suites
(6 C++ + 1 Python), GitHub Actions CI.

Acceptance criterion met: a full case runs end to end and produces a measurable crater,
with the pipeline documented and reproducible by seed.

### Stage 4 — OpenMP ⏳ next

What is missing:

1. `src/backends/OpenMpSolver.cpp` with `makeOpenMpContactSolver`. Starting point: the
   pair loop of `physics::ContactModel`, with per-thread force accumulation (a 3 × N
   array per thread, not `#pragma omp atomic`, which serializes).
2. Build with `-DGRANIMPACT_ENABLE_OPENMP=ON` (preset `omp`).
3. Add the equivalence test: same case with `--backend serial` and `--backend openmp`,
   compare the CSVs and require a relative difference < 1e-6 in a short case. That test
   is what guarantees the parallelization did not change the physics.
4. Measure scaling (1/2/4/8 threads) on `loose_vertical` and document the real measured
   speedup, with no estimates.

Acceptance criterion: `loose_vertical` (52 822 grains) complete in under an hour on an
8-core machine, with numerical equivalence against serial.

### Stage 5 — CUDA ⏳

What is missing: `src/backends/CudaSolver.cu` with the physics already written in the
core (the same `ContactModel` translated to a kernel), GPU memory and body exchange over
PCIe only; build with `-DGRANIMPACT_ENABLE_CUDA=ON`; repeat the equivalence test and
measure the full `energy_sweep`. Acceptance criterion: 1 M particles and the energy
sweep finished, with times measured and published exactly as they come out.

## 3. Computational cost

Measured on this machine (container, g++ 14, Release, 1 thread, `default` with 22 003
particles and 35 251 + 16 000 steps):

| Variant | Wall time | Cost per step |
|---|---|---|
| `--no-vtk` | 491 s (8 min 11 s) | 9.6 ms |
| with VTK every 2000 steps | 578 s (9 min 38 s) | 11.3 ms (≈15 % overhead) |

Cost per step grows roughly linearly with N (cell lists), and the number of required
steps grows as the grain size shrinks (`dt ∝ r` through the Rayleigh criterion). Practical
consequence: stage 3-4 with ~53k grains is **hours in serial** — exactly why stage 4
exists. Any speedup or per-case figure for stages 4-5 will be published when measured;
nothing is estimated here.

## 4. Pending work beyond the stages

| Task | Why | Where |
|---|---|---|
| **Bed preparation by deposition** (priority 1) | the settled nested lattice over-compacts: in the `default` case the bed ends at φ ≈ 0.82 (measured) instead of the 0.524 target, and with that the crater comes out far smaller than the Uehara reference (26.5 mm measured vs 151 mm reference). A real loose sand bed is prepared by depositing and letting it jam, not by relaxing a lattice | `generation/Generators.cpp` + settling phase in `sim/Simulation.cpp` |
| Prune the tangential history | pairs that separate leave entries that are never used again; in long runs that grows in memory | `physics/ContactModel.cpp` |
| Stiffness sensitivity study | the cases use 1e8 Pa (softened); the morphology must be shown not to depend on that choice | new `configs/stiffness_*.json` |
| Binary VTK or `.pvtu` | ASCII gets heavy above ~10⁵ particles | `io/Writers.cpp` |
| Calibration against the experiment | compare measured D(E) from the camera with the simulated crater; measure and publish the difference | `python/granimpact_post` + experimental data |
| Rolling friction | off today; it affects the angle of repose and the rim | `ContactConfig` + `ContactModel` |

## 5. Correspondence with the original plan

| Plan month | Goal | Actual status |
|---|---|---|
| 1–2 | serial + VTK output + CI | done (ahead of plan: includes analysis, tests, installer, docs) |
| 3 | full Hertz–Mindlin + morphology | done in code; the ~50k grain case still has to be run and validated against the literature |
| 4 | OpenMP, 500k particles | interface ready, solver pending |
| 5–6 | CUDA + benchmarks | interface ready, solver pending |

Note on the sandbox: the original plan says "45×45×15 cm". The 15 cm is the **bed
height** (60 % of the box's 25 cm), which is how it is configured (`domain.z = 0.25`,
`bed_height_fraction = 0.60`).
