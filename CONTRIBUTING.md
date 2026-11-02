# Contributing to GranImpact

Thank you for contributing to **GranImpact**. This document summarizes the build workflow,
testing requirements, architectural invariants, and guidelines for adding new physics or
parallel backends.

---

## 1. Build & Verification Workflow

### Configure, Build & Run All 7 Test Suites

```bash
cmake -S . -B build/release -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DGRANIMPACT_BUILD_TESTS=ON \
      -DGRANIMPACT_WARNINGS_AS_ERRORS=ON
cmake --build build/release -j$(nproc)
ctest --test-dir build/release --output-on-failure
```

### Python Post-Processing & Unit Tests

```bash
pip install -e ./python
python3 -m unittest discover -s python/tests
```

### End-to-End Smoke Run & Visualization Generation

```bash
./scripts/run_quick.sh
python3 -m granimpact_post visualize-all \
    data/outputs/loose/default/default_crater_evolution.csv
```

---

## 2. Engineering Ground Rules

1. **No Invented Numbers**: Every metric in `README.md`, `docs/`, or code comments must be
   either a documented input parameter or a measurement from a reproducible simulation run.
2. **Strict Determinism**: All stochastic sampling (polydisperse grain radii and lattice
   jitter) must flow exclusively through the seeded `std::mt19937_64` in `SimConfig`.
3. **Physics Changes Require Unit Tests**: Any modification to `src/physics/`,
   `src/generation/`, or `src/analysis/` must include a corresponding unit test in `tests/`
   (e.g., `test_contact_model.cpp`, `test_crater_analysis.cpp`).
4. **Zero Compiler Warnings**: The codebase compiles cleanly with `-Wall -Wextra -Wpedantic
   -Wshadow -Wconversion -Wsign-conversion -Wdouble-promotion -Werror`.
5. **Respect Module Layering**:
   `core ← physics ← backends`, `core ← generation`, `core ← analysis ← io`, and `sim`
   orchestrates (`docs/architecture.md`). Upward includes across layers are not permitted.
6. **Preserve Output Contracts**:
   - Simulation outputs are organized under `data/outputs/loose/<case>/` and
     `data/outputs/compact/<case>/` via `SimConfig::bedRegimeFolder()`.
   - `<case>_crater_evolution.csv`, `<case>_depth_map.csv`, `<case>_subsurface_chains.csv`,
     and `<case>_radial_evolution.csv` are consumed by `python/granimpact_post` and
     `tests/test_simulation_smoke.cpp`. Changing column names requires updating both C++
     writers and Python parsers.

---

## 3. Adding Parallel Backends (Phase 2 OpenMP / Phase 3 CUDA)

Extension points are detailed in `docs/architecture.md` and `src/backends/README.md`:
1. Implement `ContactSolver` in `src/backends/OpenMpSolver.cpp` or `src/backends/CudaSolver.cu`.
2. Register the solver in `src/backends/Backend.cpp` behind `GRANIMPACT_ENABLE_OPENMP` or
   `GRANIMPACT_ENABLE_CUDA`.
3. Add a CTest equivalence check verifying that the new backend reproduces `--backend serial`
   within `< 10⁻⁶` relative tolerance on the CSV time series for a fixed seed.
