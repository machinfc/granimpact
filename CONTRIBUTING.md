# Contributing

Thanks for your interest in GranImpact. This document describes how to build the
project, run the tests and propose changes.

## Build and test

```bash
./scripts/install.sh          # dependencies, build, ctest, install
./scripts/run_quick.sh        # end-to-end sanity check (minutes)
```

Manual equivalent:

```bash
cmake --preset release
cmake --build --preset release -j
ctest --preset release --output-on-failure
```

Python post-processing:

```bash
python3 -m pip install -e python/
python3 -m unittest discover -s python/tests
```

## Ground rules

1. **No invented numbers.** Any figure in the README, `docs/` or the code comments is
   either a measurement from a run you can reproduce or an input parameter. If it is not
   measured, write "pending".
2. **Determinism.** All randomness flows through the seeded `std::mt19937_64` in
   `SimConfig`. A change that breaks bitwise reproducibility for a fixed seed must be
   called out in the pull request.
3. **Physics changes come with a test.** Touching `physics/` or `generation/` requires a
   ctest case that pins the behaviour (see `tests/`).
4. **Warning-free.** The project builds with `-Wall -Wextra -Wpedantic -Wshadow
   -Wconversion -Wsign-conversion -Wdouble-promotion`; CI runs with
   `-DGRANIMPACT_WARNINGS_AS_ERRORS=ON`.
5. **Respect the layering.** `core ← physics ← backends`, `core ← generation`,
   `core ← analysis ← io`, `sim` orchestrates. An include that crosses the layers
   upwards is a design bug (see `docs/architecture.md`).
6. **Keep the CSV contract.** The header line and the final `# result:` block are read by
   `python/granimpact_post` and by `tests/test_simulation_smoke.cpp`; changing keys means
   changing both.

## Adding a new stage (OpenMP / CUDA)

The extension points are described in `docs/architecture.md` §2. In short: implement
`ContactSolver` in `src/backends/`, register it in `CMakeLists.txt` behind its option,
wire it into the backend factory, and add the equivalence test (same case, serial vs new
backend, relative difference < 1e-6 on the CSV series).

## Reporting bugs

Open an issue with: the config used (or the exact command line), the seed, the expected
behaviour and the actual output. If it is a numerical issue, attach the CSV or the
`# result:` line — it makes the report reproducible.
