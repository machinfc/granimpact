# Architecture

The design goal is that **adding a stage means adding code, not rewriting it**: the
simulator grows in stages (serial → OpenMP → CUDA; thousands → millions of particles)
and the only thing that changes between stages is the *contact solver*. Everything else
—configuration, bed generation, integration, analysis, outputs— is shared.

## 1. Layers

```
apps/ ──► sim/ ──► analysis/ io/      (orchestration, observables, files)
            │
            ├──► backends/ ──► physics/ ──► core/
            └──► generation/ ──────────────┘
```

| Layer | Folder | Responsibility | Depends on |
|---|---|---|---|
| Core | `core/` | types (`Real`, `Vec3`, `Box`), `ParticleSystem` (SoA), `SimConfig` (JSON) | — |
| Physics | `physics/` | Hertz–Mindlin contact law, cell lists, Velocity Verlet, boundaries | core |
| Generation | `generation/` | expanded nested-lattice bed, granular projectile | core |
| Solvers | `backends/` | `ContactSolver` interface + serial implementation (OpenMP and CUDA: stages 4-5) | physics |
| Simulation | `sim/` | `Simulation`: phases, stop criteria, results | all |
| Analysis | `analysis/` | `CraterAnalyzer`: surface, morphometry, classification, models | core |
| Outputs | `io/` | CSV series, VTK `.vtu`/`.pvd`, case summary | core, analysis |
| Executables | `apps/` | `main_serial.cpp` today; `main_openmp.cpp` / `main_cuda.cu` arrive with their backends in stages 4-5 | sim |

Rule: **the arrows never cross**. `core` includes nothing from the project, `physics`
knows nothing about the simulation and `analysis` knows nothing about the backend.
That is why the same analysis runs on a serial and on a GPU run.

## 2. Flow of a run

```
config.json
   │  SimConfig::fromFile + validate
   ▼
generateBed()            dense nested lattice -> isotropic expansion -> jitter + check
   ▼
settling                 gravity (+boost), fixed dt, double criterion: Ek per particle
   │                     and stable mean height  (no artificial damping by default)
   ▼
recordReferenceSurface() surface map of the settled bed (crater reference)
   ▼
generateProjectile()     sphere of grains with velocity (v, angle)
   ▼
impact loop              [ CellList -> ContactSolver -> integrate -> boundaries -> measure ]
   │                     every `output_every` steps: CSV row; every `vtk_every`: .vtu
   ▼
analyze()                bed drift -> deformation -> connected component -> ellipse
   │                     -> classification -> contrast models
   ▼
outputs                  <case>_crater_evolution.csv
                         <case>_summary.txt
                         <case>_<step>.vtu + <case>.pvd
```

Extension points and where to touch them:

| I want to… | Touch | Do not touch |
|---|---|---|
| Add the OpenMP stage | `src/backends/OpenMpSolver.cpp` | anything else |
| Add the CUDA stage | `src/backends/CudaSolver.cu` | anything else |
| Another crater observable | `analysis/CraterAnalyzer.{hpp,cpp}` + a column in `io/Writers.cpp` | the physics |
| Another output format | `io/Writers.{hpp,cpp}` | the simulation |
| Another analytical model | `analysis/CraterAnalyzer.cpp` (`models` namespace) | the rest |
| Another post-processing metric | `python/granimpact_post/` | C++ |

## 3. Key interfaces

```cpp
// Contact solver: the only piece that gets accelerated.
class ContactSolver {
    virtual std::string name() const = 0;
    virtual physics::ContactStats computeForces(core::ParticleSystem&, const physics::CellList&, core::Real dt) = 0;
};

// Particles in structure-of-arrays (SoA): contiguous per field, SIMD-friendly and
// copyable as-is to GPU memory.
class ParticleSystem {
    std::vector<Vec3> position_, velocity_, force_, torque_;
    std::vector<Real> radius_, mass_, inv_mass_, inertia_, inv_inertia_;
    std::vector<Material> material_;
    std::vector<ParticleKind> kind_;      // Bed | Projectile | Ejecta
};
```

`ParticleSystem` uses **SoA, not AoS**, precisely for stages 4-5: the same array the
serial loop walks with `__restrict` is uploaded to the GPU with no repacking. And
`ContactSolver` receives the already built `CellList`, so neighbor search can be
parallelized too without changing the physics.

## 4. Design decisions and why

- **One library (`granimpact_core`) + thin executables.** Tests and the analysis link
  against the same library as the executable: what is tested is what runs.
- **JSON config with its own validation** (`SimConfig::validate`): physical ranges,
  radius coherence, `t_settle > 0`, etc. An invalid config fails at startup, not in the
  middle of a 6-hour run.
- **Cell size = 2·r_max·1.05 and half the offset list**: two fewer pair visits per
  contact, and no double counting of dissipation.
- **VTK XML in ASCII**: readable, diffable and with no dependency on the VTK library.
  For large meshes, binary or `.pvtu` can be added later without touching the analysis
  (which reads from `ParticleSystem`, not from the files).
- **The CSV is self-contained**: the header carries the parameters and the end of the
  file the measured result (real impact speed, settling steps, measured φ). A results
  folder can be analyzed without the configs or the binary.
- **Flush after every row**: an interrupted run leaves usable results.
- **No invented numbers in the documentation**: everything in `docs/` comes from real
  runs or is an input parameter, and the two are clearly distinguished.

## 5. Tests

```
tests/ (ctest)                          python/tests/ (unittest)
  test_integrator        analytical free fall, energy conservation
  test_contact_model     Hertz, restitution, Coulomb, repulsion
  test_cell_list         unique and complete pairs vs brute force
  test_generation        no overlaps, density, material, seed
  test_crater_analysis   known depression, ellipse, classification
  test_simulation_smoke  full run + CSV + VTK
  test_postprocess       post-processing library (python)
```

- Tests live in `tests/` and are registered in `tests/CMakeLists.txt`.
- `ctest` runs all 7; the Python tests also run on their own with
  `python3 -m unittest discover -s python/tests` (and in CI).
- `-DGRANIMPACT_WARNINGS_AS_ERRORS=ON` turns any warning into an error: recommended in CI.
- Tests are **deterministic**: all randomness goes through a seeded `std::mt19937_64`.

## 6. Build

| Option | Effect |
|---|---|
| `GRANIMPACT_BUILD_TESTS` | build and register the tests (ON by default) |
| `GRANIMPACT_ENABLE_OPENMP` | adds `src/backends/OpenMpSolver.cpp` and links OpenMP |
| `GRANIMPACT_ENABLE_CUDA` | enables the CUDA language and `CudaSolver.cu` |
| `GRANIMPACT_ENABLE_SANITIZERS` | ASan + UBSan for development |
| `GRANIMPACT_WARNINGS_AS_ERRORS` | warnings as errors (CI) |
| `GRANIMPACT_USE_SYSTEM_DEPS` | use system nlohmann/json and gtest (no network) |

Presets: `cmake --preset debug | release | asan | omp | cuda`. Dependencies
(nlohmann/json 3.11.3 and GoogleTest 1.14.0) are downloaded with `FetchContent` the
first time; `install.sh` warns if there is no network and offers the system alternative.

Installation (`cmake --install`): `bin/granimpact_serial`, `bin/granimpact` (launcher
that picks the best available backend), `lib/libgranimpact_core.a`,
`include/granimpact/**`, `share/granimpact/{configs,docs,python}`.

## 7. Code style

- C++17, no mandatory dependencies beyond nlohmann/json.
- `namespace granimpact::{core,physics,generation,analysis,io,sim,backends}`.
- Warnings enabled: `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion
  -Wdouble-promotion …`; the project builds **warning-free** (verified).
- Comments sit at the point of the decision and explain the *why*. The two "whys" that
  cost the most time: the sign of the repulsive force and the bed-drift correction in
  the analysis (both have their own comment at the exact line, plus a test).
