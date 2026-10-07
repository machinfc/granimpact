# GranImpact

DEM (Discrete Element Method) simulator for the experimental study of impact craters in
granular media. The project grows in **stages over a single skeleton**: each stage adds a
compute backend, not a rewrite.

- **Current status:** stage 3 finished (serial core + analysis + outputs). The OpenMP and
  CUDA backends exist as interfaces only: the binary says so explicitly and the tests
  verify it.
- **Numerical method:** Hertz–Mindlin contact (Hertz normal with Di Renzo damping +
  Mindlin tangential with a Coulomb cap), Velocity Verlet integration, cell-list
  neighbour search, bed generated as an expanded nested lattice and settled under gravity.
- **No invented numbers:** every value in this README and in `docs/` is either a
  measurement from a real run or an input parameter of the experimental case. When
  something has not been measured, it is said so.
- **Known, measured limitation:** the current bed preparation over-compacts (φ ≈ 0.82
  measured against the 0.524 target), so the quick case produces a crater smaller than
  the analytical reference. It is documented in `docs/physics.md` §9.1 and it is the
  first pending task in `docs/roadmap.md` — the skeleton is ready; that part is bed
  physics, not architecture.

## Quick start

```bash
./scripts/install.sh              # dependencies, build, tests and install
./scripts/run_quick.sh            # quick case: validates the full pipeline in minutes
```

If you downloaded the project as a ZIP and the scripts do not start, it is just the
executable bit: `chmod +x scripts/*` (or run them as `bash scripts/install.sh`).

`install.sh` is idempotent: it can be re-run after every code change.

## Usage

```bash
granimpact configs/default.json --out data/outputs          # quick case (minutes)
granimpact configs/loose_vertical.json --no-vtk             # ~53k grains
granimpact --backends                                       # backends in this build
granimpact --help
```

Per-case outputs, in `--out`:

| File | Content |
|---|---|
| `<case>_crater_evolution.csv` | time series: energy, contacts, crater morphometry |
| `<case>_<step>.vtu` + `<case>.pvd` | per-particle fields for ParaView |
| `<case>_summary.txt` | human-readable summary: parameters, result, densities and timings |

## Benchmarks

Measured on a 2-core container, g++ 14, Release, serial backend.

| Case | Grains | Steps | Wall time | Notes |
|---|---|---|---|---|
| `default` (with VTK) | 22 003 | 51 251 | **470.7 s – 577.7 s** | two runs, same 2-core container; `.vtu` every 2 000 steps |
| `default` (`--no-vtk`) | 22 003 | 51 251 | **461.0 s** | the reference number for this case |
| `loose_vertical` | 52 822 | — | hours (serial) | why stage 4 (OpenMP) exists |

Wall time moves with machine load; the physics result does not (seed 20240517 reproduces
the same crater bit for bit).

Determinism check: re-running `default` with seed 20240517 reproduces the same result
(D = 26.5165 mm, d_exc = 11.1357 mm, ε = 0.9354, Z/D = 0.4199). Test suite: `ctest`
7/7 plus 8 Python unit tests.

## Stages

| Stage | Particles | Backend | Status |
|---|---|---|---|
| 1–2 | — | core, physics, analysis | done |
| 3 | 21 904 (`default`) and 52 822 (`loose_vertical`) | serial | done (the ~53k one is hours in serial) |
| 4 | 540 225 | OpenMP | pending (interface ready) |
| 5 | 1 182 447 | CUDA | pending (interface ready) |

The grain size of each stage is not arbitrary: it is chosen so the bed fills the same
thesis sandbox (0.45 × 0.45 × 0.25 m) with that stage's particle count while keeping the
target density. See `docs/configuration.md` and `docs/roadmap.md`.

## Documentation

- `docs/physics.md` — physical model, bed generation, stop criteria, observables.
- `docs/architecture.md` — layers, interfaces and recipes for adding a stage.
- `docs/configuration.md` — full JSON schema, the six cases and the CLI.
- `docs/roadmap.md` — stage plan with what is missing in each one.

## Post-processing

```bash
python3 -m granimpact_post summary data/outputs/*.csv
python3 -m granimpact_post plot data/outputs/default_crater_evolution.csv -o crater.png
python3 -m granimpact_post sweep data/outputs/energy_sweep_*.csv -o D_vs_E.png
```

## Requirements

CMake ≥ 3.16, a C++17 compiler (GCC 9+ / Clang 10+ / MSVC 2019+), git and curl (used to
fetch nlohmann/json and GoogleTest with FetchContent on the first configure). Ninja and
Python 3.9+ are optional (Ninja speeds up the build; Python enables post-processing).

## License

MIT — see [LICENSE](LICENSE).

## Citation

If you use this code in academic work, see [CITATION.cff](CITATION.cff).

