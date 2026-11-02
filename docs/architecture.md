# Architecture

The design goal of **GranImpact** is that **adding a parallel execution stage means adding
a backend solver, never rewriting the physics or analysis pipeline**. The simulator
evolves across four engineering phases:
1. **Phase 1 (Implemented)**: Serial 6-DOF DEM reference core & `craterslab` validation.
2. **Phase 2 (Next)**: OpenMP shared-memory multi-core parallelism (`src/backends/OpenMpSolver.cpp`).
3. **Phase 3 (Planned)**: NVIDIA CUDA GPU acceleration (`src/backends/CudaSolver.cu`).
4. **Phase 4 (Planned)**: Cloud HPC / AWS Batch multi-GPU scaling & full millimetric resolution.

Between phases, only the `ContactSolver` implementation changes; configuration parsing,
polydisperse bed/projectile generation, 6-DOF integration, Virtual ToF morphometry,
Archimedes subsurface field integration, and `craterslab` `.npz` exports are shared.

---

## 1. Layered Module Hierarchy

```text
apps/ ──► sim/ ──► analysis/ & io/    (orchestration, Virtual ToF, CSV/VTK/.npz)
            │
            ├──► backends/ ──► physics/ ──► core/
            └──► generation/ ─────────────────┘
```

| Layer | Folder | Responsibility | Depends on |
| :--- | :--- | :--- | :--- |
| **Core** | `core/` | Types (`Real`, `Vec3`, `Box`, `ParticleKind`), 6-DOF SoA `ParticleSystem`, `SimConfig` (JSON + CLI validation) | — |
| **Physics** | `physics/` | 6-DOF Hertz–Mindlin–Coulomb contact law, Ai et al. (2011) rolling resistance, cohesive projectile bonds, capillary/magnetic forces, flat 16-slot SoA contact history, zero-allocation CSR `CellList`, 6-DOF Velocity Verlet & boundaries | `core` |
| **Generation** | `generation/` | Expanded nested-lattice polydisperse bed, cohesive granular projectile, oblique ballistic lead-target placement | `core` |
| **Backends** | `backends/` | `ContactSolver` interface + `CpuContactSolver` (Phase 1); `OpenMpSolver` (Phase 2) and `CudaSolver` (Phase 3) | `physics` |
| **Simulation** | `sim/` | `Simulation` lifecycle: bed seeding, uniaxial top-lid compaction (`4.46 kPa`), settling convergence, projectile release, impact loop | `core`, `physics`, `generation`, `backends`, `analysis`, `io` |
| **Analysis** | `analysis/` | `CraterAnalyzer`: 3-pass Virtual ToF camera with sphere-cap occlusion, Archimedes spherical-slice `(r, z)` packing/stress fields, porous-piston & central-peak detection, thesis contrast models | `core` |
| **Outputs** | `io/` | Time-series CSV (`26` observables), 2D ToF depth map CSV, subsurface force-chain CSV, radial profile CSV, VTK `.vtu`/`.pvd` with per-grain virial stress | `core`, `analysis` |
| **Post-Processing** | `python/granimpact_post/` | 3 publication-grade scientific plots (`chains`, `rim-theory`, `formation`), `craterslab`-compatible `.npz` exporter, energy sweep fitting | CSV outputs |

**Layering Invariant**: Dependency arrows never point upward. `core` includes nothing from
the rest of the project, `physics` is agnostic of crater analysis, and `analysis` operates
solely on `ParticleSystem` arrays regardless of whether forces were computed on a single
CPU core, OpenMP threads, or a CUDA GPU.

---

## 2. End-to-End Execution Flow

```text
config.json + CLI overrides
   │  SimConfig::fromFile() + SimConfig::validate() + SimConfig::bedRegimeFolder()
   ▼
generateBed()             Dense nested lattice -> isotropic expansion to target φ
   │                      -> polydisperse radius sampling [r_min, r_max] + jitter
   ▼
settling & compaction     Gravity (+ gravity_boost) + velocity_damping
   │                      Optional uniaxial top-lid press (compaction_pressure_kPa = 4.46)
   │                      Dual stop criterion: E_k/particle < ke_tolerance AND stable <z>
   ▼
recordReferenceSurface()  3-pass Virtual ToF elevation map z_ref(x, y) of settled bed
   ▼
generateProjectile()      Spherical granular cluster (R_p, φ_p, ρ_p) with mass calibration
   │                      + ballistic lead-target offset (x_0, z_0) for oblique impacts
   │                      + ContactModel::registerCohesiveBonds() (σ_y = 17.15 kPa)
   ▼
impact loop               1. Integrator::stepFirstHalf()   (6-DOF kick + drift)
   │                      2. CellList::build()             (zero-alloc O(N) CSR grid)
   │                      3. ContactSolver::computeForces() (Hertz-Mindlin + rolling + bonds)
   │                      4. Integrator::stepSecondHalf()  (6-DOF second kick)
   │                      5. Integrator::applyBoundaries() (floor, walls, ejecta ceiling)
   │                      Every output_every steps: CraterAnalyzer::analyze() -> CSV row
   │                      Every vtk_every steps:    VtkWriter::write()        -> .vtu
   ▼
final diagnostics         <case>_crater_evolution.csv  (time series + # result: footer)
                          <case>_depth_map.csv         (96x96 ToF map for craterslab .npz)
                          <case>_subsurface_chains.csv (2D cylindrical (r,z) φ & σ_virial)
                          <case>_radial_evolution.csv  (azimuthal profiles z_bed & z_total)
                          <case>_summary.txt           (human-readable physical report)
                          <case>_<step>.vtu + <case>.pvd (3D ParaView animation)
```

### Extension Points

| Engineering Task | Files to Modify | Untouched Modules |
| :--- | :--- | :--- |
| Add Phase 2 OpenMP solver | `src/backends/OpenMpSolver.cpp`, `src/backends/Backend.cpp` | `core`, `generation`, `analysis`, `io` |
| Add Phase 3 CUDA solver | `src/backends/CudaSolver.cu`, `src/backends/Backend.cpp` | `core`, `generation`, `analysis`, `io` |
| Add a new contact force term | `include/granimpact/physics/ContactModel.hpp`, `src/physics/ContactModel.cpp` | `CellList`, `CraterAnalyzer`, `Writers` |
| Add a new crater observable | `include/granimpact/analysis/CraterAnalyzer.hpp`, `src/analysis/CraterAnalyzer.cpp`, `src/io/Writers.cpp` | `physics`, `generation`, `backends` |
| Add a new post-processing plot | `python/granimpact_post/plots.py`, `python/granimpact_post/cli.py` | Entire C++ codebase |

---

## 3. Key HPC Data Structures & Interfaces

```cpp
// Contact solver interface: the single hot-spot kernel swapped across Phases 1-3.
class ContactSolver {
public:
    virtual ~ContactSolver() = default;
    virtual std::string name() const = 0;
    virtual physics::ContactStats computeForces(
        core::ParticleSystem& particles,
        const physics::CellList& cells,
        core::Real dt) = 0;
};

// 6-DOF Structure-of-Arrays (SoA) particle storage: contiguous per field,
// SIMD-friendly, and directly copyable to GPU device buffers.
class ParticleSystem {
    std::vector<Vec3> position_, velocity_, omega_, force_, torque_;
    std::vector<Real> radius_, mass_, inv_mass_, inertia_, inv_inertia_;
    std::vector<ParticleKind> kind_;             // Bed | Projectile | Ejecta
    std::vector<Real> young_, poisson_, friction_, restitution_;
    std::vector<Real> normal_force_sum_;         // Σ|F_n| per grain for virial stress
    std::vector<Index> coordination_;            // active contact count Z_i per grain
};
```

### Why These Architectural Choices Matter for HPC

1. **Pure 6-DOF Structure-of-Arrays (`ParticleSystem`)**: Rather than storing an Array of
   Structures (`std::vector<Particle>`), every physical field (`position_`, `velocity_`,
   `omega_`, `force_`, `torque_`, `normal_force_sum_`, `coordination_`) is stored in its
   own contiguous vector. This maximizes CPU cache-line utilization in serial/OpenMP loops
   and enables coalesced global memory transactions in CUDA without data repacking.
2. **Flat 16-Slot Per-Particle Contact History (`ContactModel::history_slots_`)**:
   Traditional DEM codes store Mindlin tangential history and cohesive bonds in a hash
   table (`std::unordered_map<uint64_t, ContactState>`), which causes pointer chasing,
   heap allocations every timestep, and thread contention. `ContactModel` replaces this
   with a flat contiguous array of size `N × 16` (`kMaxContactsPerParticle = 16`, stamped
   by `step_stamp_`), delivering a `1.89×` serial speedup and mapping directly to lock-free
   OpenMP and CUDA kernels.
3. **Zero-Allocation CSR Neighbor Search (`CellList`)**: Cell head/next buffers are
   reused across timesteps without heap reallocation. Pairs are visited canonically
   (`i < j`) across the `14` forward neighbor offsets (out of `27`), halving pair checks.
4. **Structured Regime Output Isolation (`SimConfig::bedRegimeFolder()`)**: Outputs are
   automatically routed into `data/outputs/loose/<case>/` or `data/outputs/compact/<case>/`
   based on `compaction_pressure_kPa` and `packing_fraction`, keeping loose and compacted
   experimental campaigns cleanly separated.

---

## 4. Verification & Test Suite Architecture

All `7` automated test suites link directly against `granimpact_core` and run via `ctest`:

```text
tests/ (C++17 GoogleTest)                python/tests/ (Python unittest)
  1. test_integrator                     Analytical free fall, symplectic conservation
  2. test_contact_model                  Hertz, restitution, Mindlin-Coulomb, rolling,
                                         cohesive bond fracture, capillary & magnetic
  3. test_cell_list                      Exact pair equivalence vs O(N^2) brute force
  4. test_generation                     Zero overlaps, polydispersity, oblique lead
  5. test_crater_analysis                3-pass Virtual ToF, Archimedes subsurface,
                                         Heckel & double-log rim saturation models
  6. test_simulation_smoke               End-to-end loose & compact runs + CSV/VTK
  7. test_postprocess                    CSV parser, craterslab .npz export, CLI
```

---

## 5. Build System & CMake Options

| CMake Option | Default | Effect |
| :--- | :---: | :--- |
| `GRANIMPACT_BUILD_TESTS` | `ON` | Builds and registers all 7 C++ & Python test suites in CTest |
| `GRANIMPACT_ENABLE_OPENMP` | `OFF` | Enables OpenMP compilation and `OpenMpSolver` registration |
| `GRANIMPACT_ENABLE_CUDA` | `OFF` | Enables CUDA language support and `CudaSolver` registration |
| `GRANIMPACT_ENABLE_SANITIZERS` | `OFF` | Enables AddressSanitizer + UndefinedBehaviorSanitizer (`ASan`/`UBSan`) |
| `GRANIMPACT_WARNINGS_AS_ERRORS` | `OFF` | Passes `-Werror` on top of `-Wall -Wextra -Wpedantic -Wconversion` |
| `GRANIMPACT_USE_SYSTEM_DEPS` | `OFF` | Uses system `nlohmann_json` and `GTest` instead of `FetchContent` |
