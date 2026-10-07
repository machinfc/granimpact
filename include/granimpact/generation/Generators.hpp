#pragma once
// Granular bed and projectile generation.
//
// Bed method: cubic lattice with alternate layers shifted by half a cell.
// Why this and not a "rain" (random deposition):
// * NO overlaps by construction -> no false elastic energies appear;
// * it preserves the full radius distribution [r_min, r_max];
// * it is O(N) and deterministic for a fixed seed, so cases are repeatable;
// * the density it leaves is a realistic loose bed (phi ~ 0.45-0.52), and the
// settling phase of Simulation compacts it to its physical value.
// The final density is therefore a MEASURED RESULT (reported in the log and in
// the CSV), not an input datum: `packing_fraction` in the config is the TARGET
// that sets the particle count and the seeded region.
#include <random>

#include "granimpact/core/ParticleSystem.hpp"
#include "granimpact/core/SimConfig.hpp"

namespace granimpact::generation {

struct BedInfo {
    std::size_t count{0};        // particles actually seeded
    core::Real region_height{0};  // height of the seeded region [m]
    core::Real target_phi{0};     // target from the config
    core::Real seeded_phi{0};     // measured after seeding (before settling)
    core::Real lattice_spacing{0}; // horizontal lattice spacing [m]
    core::Real layer_spacing{0};  // vertical spacing between layers [m]
};

// Seeds the bed in the region [box.min.z, box.min.z + bed_height_fraction * H].
BedInfo generateBed(core::ParticleSystem& particles, const core::SimConfig& config,
                    std::mt19937_64& rng);

// Granular projectile: a sphere of radius R filled with grains at initial velocity
// (speed, angle_deg). The same lattice is reused with spherical sampling.
std::size_t generateProjectile(core::ParticleSystem& particles, const core::SimConfig& config,
                               std::mt19937_64& rng);

// Measured density of a set of grains inside a given volume.
core::Real measuredPacking(const core::ParticleSystem& particles, core::Real volume);

}  // namespace granimpact::generation
