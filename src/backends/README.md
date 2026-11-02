# Contact Solver Backends (`src/backends/`)

The `ContactSolver` interface (`include/granimpact/sim/ContactSolvers.hpp`) isolates the
pair-force evaluation kernel from the rest of the simulator. Configuration parsing,
bed/projectile generation, 6-DOF integration, Virtual ToF crater morphometry, Archimedes
subsurface fields, and CSV/VTK/`.npz` writers are shared across all backends.

| File | Roadmap Phase | Status | Description |
| :--- | :---: | :---: | :--- |
| `SerialSolver.hpp` + `Backend.cpp` | Phase 1 | ✅ Implemented | Single-thread 6-DOF CPU reference solver (`--backend serial`) |
| `OpenMpSolver.cpp` | Phase 2 | ⏳ Planned | Shared-memory multi-core CPU solver (`--backend openmp`) |
| `CudaSolver.cu` | Phase 3 | ⏳ Planned | NVIDIA CUDA GPU solver (`--backend cuda`) |

---

## 1. How `ContactModel` Is Structured for Lock-Free Parallelism

`physics::ContactModel` (`include/granimpact/physics/ContactModel.hpp`) is already designed
for direct OpenMP and CUDA parallelization:
- **6-DOF SoA Particle Arrays**: `ParticleSystem` stores `position_`, `velocity_`, `omega_`,
  `force_`, `torque_`, `normal_force_sum_`, and `coordination_` in contiguous vectors.
- **Flat 16-Slot Per-Particle Contact History (`kMaxContactsPerParticle = 16`)**: Tangential
  Mindlin displacements (`delta_t`) and breakable cohesive bond flags (`cohesive`) are stored
  in a contiguous array `history_slots_` of size `N × 16` indexed by the lower particle index
  `i < j`. No `std::unordered_map` or dynamic heap allocation occurs during `computeForces()`.

---

## 2. Adding Phase 2 (`OpenMpSolver.cpp`)

1. Implement `std::unique_ptr<ContactSolver> makeOpenMpContactSolver(const core::SimConfig&)`
   in `src/backends/OpenMpSolver.cpp`.
2. Parallelize the cell/pair traversal using thread-private force/torque reduction buffers
   (or 8-color spatial cell coloring) while updating `history_slots_[i * 16 + slot]`
   without locks (since each particle `i` owns its 16 slots for `i < j`).
3. Configure with `-DGRANIMPACT_ENABLE_OPENMP=ON` and verify bitwise/numerical equivalence
   against `--backend serial` (`relative difference < 10⁻⁶`) via `ctest`.

---

## 3. Adding Phase 3 (`CudaSolver.cu`)

1. Implement `std::unique_ptr<ContactSolver> makeCudaContactSolver(const core::SimConfig&)`
   in `src/backends/CudaSolver.cu`.
2. Mirror the SoA `ParticleSystem` arrays and flat `N × 16` `ContactSlot` table in GPU device
   memory, build cell lists on device via radix sort, and transfer particle states back to
   host only on `output_every` / `vtk_every` steps.
3. Configure with `-DGRANIMPACT_ENABLE_CUDA=ON` and run the backend equivalence suite.
