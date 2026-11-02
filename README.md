# GranImpact — 3D 6-DOF Discrete Element Method for Granular Impact Cratering

[![CI](https://github.com/machinfc/granimpact/actions/workflows/ci.yml/badge.svg)](https://github.com/machinfc/granimpact/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/machinfc/granimpact/blob/main/LICENSE)

## 1. Overview & Purpose

**GranImpact** is a 3D 6-DOF Discrete Element Method (DEM) simulation engine in C++17 and
Python post-processing suite for low-velocity granular-granular impact cratering on loose,
pre-compacted, oblique, and lunar-gravity beds (`craterslab` `.npz` compatible).

| Capability / Physical Regime | Implementation & Validated Result |
| :--- | :--- |
| **6-DOF Contact Mechanics** | Hertz–Tsuji normal viscoelasticity + Mindlin–Deresiewicz tangential history + Coulomb cap + Ai et al. (2011) rolling resistance (`μᵣ = 0.15–0.18`). |
| **Consolidated Projectile Comminution** | Breakable Portland cement bridges (`σᵧ = 17.15 kPa`, `φₚ = 0.50`, `Mₚ = 243.46 g`) capturing Regime 1 (`P₀ < σᵧ`, Simple) → Regime 2 (`P₀ ≥ σᵧ`, Porous Piston & Central Peak). |
| **Uniaxial Heckel Bed Compaction** | Loose bed (`φ₀ = 0.527`, `H = 109.0 mm`) vs. uniaxial pre-compacted bed (`4.46 kPa` top press, `H = 85.4 mm`) saturating at Random Close Packing (`φ_RCP = 0.64`). |
| **Virtual ToF Camera & `craterslab`** | 3-pass optical sphere-cap occlusion + pore filter (`96 × 96` grid), Archimedes spherical-slice `(r, z)` stress/packing fields, and direct `.npz` export for `craterslab`. |
| **HPC-Ready SoA Architecture** | Contiguous Structure-of-Arrays (`ParticleSystem`), zero-allocation `O(N)` `CellList`, and flat 16-slot per-particle contact table (`ContactModel`) ready for OpenMP & CUDA. |

```text
granimpact/
├── CMakeLists.txt                   # C++17 build, strict warnings, GoogleTest
├── configs/                         # Experimental & planetary JSON scenarios
│   ├── default.json                 # Fast pipeline validation (loose/default)
│   ├── loose_vertical.json          # Loose bed (phi_0 = 0.524), 90 deg vertical
│   ├── compact_vertical.json        # Compacted bed (phi_0 = 0.574, 4.46 kPa), 90 deg
│   ├── loose_oblique45.json         # Loose bed, 45 deg oblique impact
│   ├── compact_oblique45.json       # Compacted bed, 45 deg oblique impact
│   ├── lunar_impact.json            # Lunar gravity (g = -1.62 m/s^2), 75 deg
│   ├── lunar_oblique45.json         # Lunar gravity (g = -1.62 m/s^2), 45 deg oblique
│   └── energy_sweep.json            # Parametric kinetic energy sweep template
├── include/granimpact/              # Public headers (core, physics, analysis, io, sim)
├── src/                             # C++17 engine & CLI implementation
├── python/granimpact_post/          # CLI post-processing, 3 scientific plots & .npz
├── tests/                           # 7 unit & integration test suites (C++ & Python)
└── data/outputs/{loose,compact}/    # Structured per-regime simulation outputs
```

---

## 2. Build & Run

### 2.1. Prerequisites

| Dependency | Version | Purpose |
| :--- | :--- | :--- |
| **C++17 Compiler** | `g++ ≥ 9` or `clang++ ≥ 10` | Compiles `granimpact_core` and `granimpact_serial` with `-O3`. |
| **CMake & Ninja** | `cmake ≥ 3.16`, `ninja` | Build system and `FetchContent` for `nlohmann_json` & `googletest`. |
| **Python 3** | `python ≥ 3.10`, `numpy`, `matplotlib` | Runs `granimpact_post` plots and `craterslab` `.npz` export. |
| **ParaView** *(optional)* | `paraview ≥ 5.10` | Interactive 3D visualization of `.pvd` / `.vtu` particle & stress fields. |

### 2.2. Configure, Build & Install

```bash
cmake -B build/release -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DGRANIMPACT_WARNINGS_AS_ERRORS=ON
cmake --build build/release -j$(nproc)
pip install -e ./python
```

### 2.3. Run Canonical Scenarios & Generate Visualizations / `.npz`

```bash
# 1. Quick validation run (~25.8k bed grains, loose bed)
./build/release/bin/granimpact_serial configs/default.json

# 2. Loose vs. Compacted vertical & 45-deg oblique impacts
./build/release/bin/granimpact_serial configs/loose_vertical.json
./build/release/bin/granimpact_serial configs/compact_vertical.json
./build/release/bin/granimpact_serial configs/loose_oblique45.json
./build/release/bin/granimpact_serial configs/compact_oblique45.json

# 3. Lunar gravity oblique impacts (g = -1.62 m/s^2)
./build/release/bin/granimpact_serial configs/lunar_impact.json
./build/release/bin/granimpact_serial configs/lunar_oblique45.json

# 4. Generate all 3 scientific figures + craterslab .npz
python3 -m granimpact_post visualize-all \
    data/outputs/loose/default/default_crater_evolution.csv
```

---

## 3. Verification, Testing & Benchmarks

### 3.1. Automated Verification Suite (`CTest` & `unittest`)

```bash
ctest --test-dir build/release --output-on-failure
```

All **7 test suites** (`test_integrator`, `test_contact_model`, `test_cell_list`,
`test_generation`, `test_crater_analysis`, `test_simulation_smoke`, `test_postprocess`)
verify symplectic energy conservation, Hertz–Mindlin–Coulomb forces, rolling resistance,
breakable cohesive bonds, capillary/magnetic forces, `O(N)` neighbor search equivalence
against `O(N²)` brute force, Archimedes subsurface integration, and `.npz` export.

### 3.2. Physical Accuracy Validation (`450 × 450 × 250 mm` Sandbox, `N_bed = 25,857`)

| Observable / Metric | Loose Regime 1 (`v₀ = 2.54 m/s`) | Loose Regime 2 (`v₀ = 7.02 m/s`) | Compact Regime 2 (`v₀ = 6.97 m/s`) | Physical Reference (PRE 2024 & Thesis) |
| :--- | :---: | :---: | :---: | :--- |
| **Initial Packing (φ₀)** | `0.527` | `0.527` | `0.640` *(4.46 kPa lid)* | Loose (`φ_L = 0.524`) vs. Compacted (`φ_C = 0.574`). |
| **Impact Pressure (P₀)** | `4.09 kPa < σᵧ` | `32.10 kPa ≥ σᵧ` | `32.10 kPa ≥ σᵧ` | Transition at projectile yield `σᵧ = 17.15 kPa`. |
| **Peak Virial Stress (σ_max)** | `54.8 kPa` | `160.2 kPa` | `151.4 kPa` | Love–Weber force chains exceed `σᵧ` in Regime 2. |
| **Peak Subsurface Packing (φ_sub)** | `0.587` *(partial)* | `0.640` *(RCP jam)* | `0.640` *(RCP jam)* | Uniaxial Heckel saturation at `φ_RCP = 0.64`. |
| **Crater Diameter (D)** | `122.4 mm` | `99.0 mm` | `174.0 mm` | Compacted plug redirects vertical momentum laterally. |
| **Excavated Depth (d_exc)** | `13.48 mm` | `22.55 mm` | `34.80 mm` | Uniaxial floor subsidence + lateral excavation. |
| **Rim Height (h_rim)** | `15.21 mm` | `20.33 mm` | `18.42 mm` | Double-log saturation (`h_rim ∝ log(log E)`). |
| **Central Peak Height (h_cp)** | `0.00 mm` | `13.98 mm` | `13.20 mm` | Porous piston forms central peak only for `P₀ ≥ σᵧ`. |
| **Morphological Class** | `Simple` | `Complex (central peak)` | `Complex (central peak)` | Automated morphometric transition detection. |

### 3.3. End-to-End Serial Optimization & Scaling Roadmap

| Optimization / Scaling Stage | Contact History & Spatial Search | Particle Count (N) | Step Cost / Status |
| :--- | :--- | :---: | :--- |
| **Baseline Serial (AoS + Hash Map)** | `std::unordered_map` per-pair tangential state | `22,003` | `11.3 ms/step` (`577.7 s` total) |
| **Optimized Serial (Flat SoA + O(1) Radius)** | Pre-allocated 16-slot SoA table, zero-alloc `CellList` | `25,960` | `8.9 ms/step` (`304.9 s` total, `1.89×` speedup) |
| **Strong Scaling (OpenMP — Phase 2)** | Thread-private force reduction / cell coloring | `58,926 – 520,000` | Interface ready (`--backend openmp`) |
| **Weak Scaling (CUDA — Phase 3)** | Coalesced SoA device arrays & cell-sorted kernels | `1,750,000+` | Interface ready (`--backend cuda`) |
| **Historical Reference Comparison** | Melosh (1989) hydrocodes & Uehara (2003) solid-sphere scaling (`Z ∝ H¹/³`) omit porous target compaction and projectile comminution; GranImpact reproduces Heckel uniaxial saturation and porous-piston central peaks. | — | Verified against *Phys. Rev. E* **109**, 054907 (2024) |

---

## 4. Engineering Roadmap Status & Technical Limitations

### ✅ Implemented & Verified (Phase 1 — Serial Reference Core)
- **Full 6-DOF DEM Physics**: Hertz–Tsuji normal viscoelasticity, Mindlin–Deresiewicz
  tangential history, Coulomb sliding, Ai et al. (2011) rolling resistance, and breakable
  cement bonds (`σᵧ = 17.15 kPa`) for weakly consolidated projectiles.
- **Uniaxial Bed Compaction & Extensions**: Top-lid uniaxial press (`4.46 kPa`) during
  settling, oblique ballistic lead targeting (`45°` & `75°`), capillary/moisture bed
  cohesion (`bed_cohesion_kPa`), and dipole–dipole magnetic sand forces.
- **Virtual ToF Sensor & `craterslab` Integration**: 3-pass surface reconstruction with
  sphere-cap optical occlusion, Archimedes spherical-slice `(r, z)` subsurface fields,
  structured outputs (`data/outputs/{loose,compact}/<case>/`), and `.npz` export.

### ❌ Known Architectural Boundaries (Scheduled for Phases 2–4)
- **OpenMP & CUDA Backends**: `OpenMpSolver` and `CudaSolver` interfaces are wired into
  `backends/ContactSolvers.hpp` and validated by CLI/tests, awaiting kernel implementation
  in Phases 2 and 3.
- **Softened Contact Modulus (`E = 10⁸ Pa`)**: Standard DEM stiffness softening is used to
  keep the Rayleigh time step (`Δt = 2–5 μs`) tractable in serial runs; full `10⁶`-grain
  millimetric beds (`r = 0.4–0.6 mm`) require the Phase 3 CUDA backend.
