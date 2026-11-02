# Engineering Roadmap & Scaling Ladder

**GranImpact** is structured around four incremental HPC & physics milestones. Because the
core engine uses contiguous 6-DOF Structure-of-Arrays (`ParticleSystem`), zero-allocation
CSR cell lists (`CellList`), and a flat 16-slot per-particle contact table (`ContactModel`),
advancing from Phase 1 to Phases 2–4 only requires adding backend solvers in `src/backends/`.

---

## 1. Polydisperse Particle Resolution Ladder

All stages model the same `450 × 450 × 250 mm` experimental sandbox (`60%` initial fill height)
and `D_p = 70.7 mm` (`M_p = 243.46 g`, `σᵧ = 17.15 kPa`) projectile. Varying the polydisperse
grain radius range `[r_min, r_max]` (`r_max / r_min = 1.35`) scales the total particle count
`N ∝ 1 / <r³>` and the Rayleigh timestep `Δt ∝ r_min`:

| Roadmap Phase | Radii `[r_min, r_max]` | Bed Grains (Loose / Compact) | Projectile Grains | Rayleigh `Δt` | Target Compute Backend |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **Phase 1 (Fast Baseline)** | `4.10–5.50 mm` | `25,857` / `29,348` | `103` | `5.0 μs` | Serial CPU (`granimpact_serial`) |
| **Phase 1–2 (Thesis Replica)** | `3.10–4.20 mm` | `58,926` / `66,735` | `233` | `2.0 μs` | Serial CPU → OpenMP Multi-Core |
| **Phase 2 (Intermediate)** | `1.50–2.00 mm` | `≈ 520,000` | `≈ 2,100` | `1.0 μs` | OpenMP Shared-Memory (`8–32` cores) |
| **Phase 3 (Fine Grain)** | `1.00–1.35 mm` | `≈ 1,750,000` | `≈ 7,000` | `0.6 μs` | NVIDIA CUDA GPU (`sm_75+`) |
| **Phase 4 (Full Sand Scale)** | `0.40–0.60 mm` | `≈ 24,000,000` | `≈ 95,000` | `0.25 μs` | Multi-GPU / AWS Cloud HPC |

---

## 2. Status by Engineering Phase

### Phase 1 — Serial 6-DOF Physics & `craterslab` Validation Baseline (Completed ✅)

Delivered and verified across all `7` CTest suites:
- **6-DOF Rotational & Contact Mechanics**: Symplectic 6-DOF Velocity Verlet, Hertz–Tsuji
  normal viscoelasticity, Mindlin–Deresiewicz tangential history, Coulomb sliding, and
  Ai et al. (2011) rolling resistance (`μ_r = 0.15–0.18`).
- **Consolidated Projectile Comminution (`σᵧ = 17.15 kPa`)**: Breakable cement bridges
  capturing the transition from Regime 1 (`P₀ < σᵧ`, `Simple` crater) to Regime 2
  (`P₀ ≥ σᵧ`, `Complex (central peak)` porous piston).
- **Loose vs. Pre-Compacted Beds & Planetary Extensions**: Uniaxial top-lid press (`4.46 kPa`)
  during settling (`φ₀ = 0.527` loose vs. `0.640` compact), oblique ballistic lead targeting
  (`45°`, `75°`), lunar gravity (`gz = -1.62 m/s²`), capillary bed cohesion (`bed_cohesion_kPa`),
  and dipole–dipole magnetic sand forces (`magnetic_susceptibility`, `magnetic_field_T`).
- **Virtual ToF Camera, Archimedes Subsurface Fields & `craterslab` Export**: 3-pass
  optical sphere-cap surface reconstruction (`96 × 96`), analytical spherical-slice `(r, z)`
  packing/stress integration, structured outputs (`data/outputs/{loose,compact}/<case>/`),
  3 automated scientific figures, and direct `<case>_craterslab.npz` export.

### Phase 2 — OpenMP Shared-Memory Multi-Core Parallelism (Next ⏳)

Implementation steps:
1. Implement `makeOpenMpContactSolver` in `src/backends/OpenMpSolver.cpp` using the flat
   16-slot SoA `history_slots_` table (each particle `i` exclusively owns its `16` slots for
   `i < j`, eliminating hash-table locks) and thread-private force/torque reduction buffers
   or spatial cell coloring.
2. Enable with `-DGRANIMPACT_ENABLE_OPENMP=ON`.
3. Add a deterministic equivalence test in `tests/` requiring relative difference `< 10⁻⁶`
   between `--backend serial` and `--backend openmp` on identical seeds.
4. Benchmark strong scaling (`1, 2, 4, 8, 16` threads) on `loose_vertical.json` (`58,926`
   grains) and `compact_vertical.json` (`66,735` grains).

### Phase 3 — NVIDIA CUDA GPU Acceleration (Planned ⏳)

Implementation steps:
1. Implement `makeCudaContactSolver` in `src/backends/CudaSolver.cu` uploading the SoA
   `ParticleSystem` vectors and flat `N × 16` contact history array to device memory.
2. Build cell lists on GPU via radix sort / CSR cell-start arrays (`thrust` / `CUB`) and
   evaluate 6-DOF contacts in warp-coalesced device kernels.
3. Keep particle arrays resident on the GPU across timesteps, transferring only when
   `output_every` or `vtk_every` triggers host analysis/IO.
4. Validate numerical equivalence against serial/OpenMP and execute full `10⁶`-particle
   energy sweeps (`scripts/energy_sweep.py`).

### Phase 4 — Cloud HPC / AWS Batch & Full Millimetric Sand Scale (Planned ⏳)

Implementation steps:
1. Containerize the CUDA build (`Apptainer` / `Docker`) for AWS GPU instances (`G4dn` / `G5` / `P4d`).
2. Automate multi-parameter sweeps (`v₀`, `angle_deg`, `φ₀`, `bed_cohesion_kPa`, `B₀`) with
   automated `craterslab` `.npz` archiving and binary VTK (`.pvtu`) output.

---

## 3. Measured Serial Optimization Progression

Measured on Linux (`g++ 14 -O3 -DNDEBUG`, single core, `N ≈ 22k–26k` grains, `51,251` total
steps including `t_settle` + `t_end` and VTK output every `2,000` steps):

| Engine Revision | Contact History & Spatial Search | Particle Count (N) | Step Cost | Total Wall Time | Speedup |
| :--- | :--- | :---: | :---: | :---: | :---: |
| **Initial Serial Baseline** | `std::unordered_map` hash table + `O(N)` max-radius scan | `22,003` | `11.3 ms` | `577.7 s` | `1.00×` |
| **Optimized Serial SoA** | Flat 16-slot SoA `history_slots_` + `O(1)` cached `max_radius_` + precomputed damping | `25,960` | `8.9 ms` | `304.9 s` | **`1.89×`** |
