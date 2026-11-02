#pragma once
// Hertz-Mindlin-Deresiewicz 6-DOF contact model + Ai et al. (2011) rolling resistance
// + breakable cohesive bonds for weakly consolidated granular projectiles.
// References:
// - Hertz (1882), Mindlin & Deresiewicz (1953), Tsuji et al. (1992)
// - Di Renzo & Di Maio (2004), Luding (2008)
// - Ai, Chen, Rotter & Ooi, Powder Technology 206, 269-282 (2011)
// - Corrales-Machin et al., Phys. Rev. E 109, 054907 (2024)
#include <cstdint>
#include <vector>

#include "granimpact/core/ParticleSystem.hpp"
#include "granimpact/core/SimConfig.hpp"
#include "granimpact/physics/CellList.hpp"

namespace granimpact::physics {

using core::Index;
using core::Real;
using core::Vec3;

struct ContactStats {
    std::size_t candidates{0};
    std::size_t touching{0};
    std::size_t sliding{0};
    std::size_t active_bonds{0};
    std::size_t broken_bonds{0};
    double max_overlap{0.0};
    double max_normal_force{0.0};
    double dissipated_energy{0.0};
};

// Effective properties of a pair (i, j), exposed so they can be tested alone.
struct EffectiveProperties {
    Real e_star{0};   // effective Young's modulus
    Real g_star{0};   // effective shear modulus
    Real r_star{0};   // effective radius
    Real m_eff{0};    // effective mass
};

class ContactModel {
public:
    // Fixed number of contact history slots per particle (stored on the lower index i < j).
    // In 3D sphere packings the kissing number is 12 (~6 with i < j), so 16 slots per
    // particle provides generous headroom while mapping 1:1 to OpenMP and CUDA SoA arrays.
    static constexpr std::size_t kMaxContactsPerParticle = 16;

    explicit ContactModel(const core::SimConfig& config);

    // Computes and accumulates every contact force and torque on the particles.
    // `dt` is used to integrate the tangential displacement (Mindlin) and rolling resistance.
    ContactStats computeForces(core::ParticleSystem& particles, const CellList& cells, Real dt);

    // Registers breakable cohesive cement bonds between neighbouring Projectile grains.
    void registerCohesiveBonds(const core::ParticleSystem& particles, Real max_gap_fraction = 0.05);

    // Exposed for the unit tests.
    [[nodiscard]] static EffectiveProperties effectiveProperties(
        Real ri, Real rj, Real ei, Real ej, Real ni, Real nj, Real mi, Real mj);

    // Hertz normal force (undamped): F = (4/3) E* sqrt(R*) delta^(3/2)
    [[nodiscard]] static Real hertzNormalForce(const EffectiveProperties& props, Real overlap);

    // Normal damping coefficient from the restitution coefficient.
    [[nodiscard]] static Real normalDamping(Real restitution, Real m_eff, Real stiffness);

    [[nodiscard]] const ContactStats& stats() const { return stats_; }
    [[nodiscard]] double dissipatedEnergy() const { return dissipated_energy_; }
    [[nodiscard]] std::size_t activeBonds() const { return active_bonds_; }
    [[nodiscard]] std::size_t brokenBonds() const { return broken_bonds_; }
    void clearHistory();

private:
    struct ContactSlot {
        core::Index partner{-1};
        std::uint32_t last_step{0};
        core::Vec3 delta_t{};
        std::uint8_t cohesive{0};
    };

    void ensureCapacity(std::size_t n_particles);
    ContactSlot* findOrInsertSlot(std::size_t i, core::Index j, bool create_if_missing);

    Real restitution_{0};
    Real damping_factor_{0};              // precomputed -2 ln(e) / sqrt(pi^2 + ln^2(e))
    bool tangential_enabled_{true};
    bool rolling_friction_enabled_{false};
    Real rolling_mu_{0.15};
    Real tangential_ratio_{1.0};
    bool cohesive_projectile_{true};
    Real yield_stress_Pa_{0};             // projectile compressive/tensile yield stress [Pa]
    Real bed_cohesion_Pa_{0};             // capillary / moisture / van der Waals cohesion in the bed [Pa]
    Real magnetic_chi_{0};                // magnetic susceptibility chi_m
    Real magnetic_b_T_{0};                // ambient magnetic field B_0 [T]
    bool periodic_xy_{false};
    core::Vec3 domain_size_{};

    std::uint32_t step_stamp_{0};
    double dissipated_energy_{0.0};
    std::size_t active_bonds_{0};
    std::size_t broken_bonds_{0};
    bool bonds_initialized_{false};

    ContactStats stats_{};
    // Flat per-particle contact history table of size N * kMaxContactsPerParticle.
    std::vector<ContactSlot> history_slots_;
    std::vector<std::uint8_t> has_cohesive_;
};

}  // namespace granimpact::physics
