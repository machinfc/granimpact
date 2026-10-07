# Backends

The contact solver is the only piece that changes between stages. Everything else
(configuration, bed generation, integration, analysis, outputs) is shared.

| File | Stage | Status |
|---|---|---|
| `SerialSolver.hpp` + `Backend.cpp` | 3 | implemented |
| `OpenMpSolver.cpp` | 4 | pending (skeleton ready) |
| `CudaSolver.cu` | 5 | pending (skeleton ready) |

## How to add the OpenMP stage

1. Create `src/backends/OpenMpSolver.cpp` with
   `std::unique_ptr<ContactSolver> makeOpenMpContactSolver(const core::SimConfig&)`.
   The brute-force loop of `physics::ContactModel` is the starting point: parallelize
   the pair loop and accumulate forces per thread (avoid races with a per-thread
   `std::vector<Vec3>`, not with `#pragma omp atomic`).
2. Build with `-DGRANIMPACT_ENABLE_OPENMP=ON`.
3. `ctest` already runs the same suites for `serial` and `openmp` (the cases are
   parameterized by backend): if the results do not agree within tolerance, the test
   fails.

## How to add the CUDA stage

Same, but `makeCudaContactSolver` in `src/backends/CudaSolver.cu`, with the memory
policy on the GPU and only body exchange over PCIe.
