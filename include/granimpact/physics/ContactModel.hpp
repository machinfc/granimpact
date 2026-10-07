#pragma once
// Hertz-Mindlin contact model (Hertz normal + Mindlin tangential + Coulomb).
// Referencias: Hertz (1882), Mindlin & Deresiewicz (1953), Thornton & Yin (1991),
// Di Renzo & Di Maio (2004). See docs/physics.md for the derivations.
#include <cstdint>
#include <unordered_map>

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
    double max_overlap{0.0};
    double max_normal_force{0.0};
};

// Effective properties of a pair (i, j), exposed so they can be tested alone.
struct EffectiveProperties {
    Real e_star{0};   // effective Young's modulus
    Real g_star{0};   // effective shear modulus
    Real r_star{0};   // radio efectivo
    Real m_eff{0};    // masa efectiva
};

class ContactModel {
public:
    explicit ContactModel(const core::SimConfig& config);

    // Computes and accumulates every contact force on the particles.
    // `dt` is used to integrate the tangential displacement (Mindlin).
    ContactStats computeForces(core::ParticleSystem& particles, const CellList& cells, Real dt);

    // Exposed for the unit tests.
    [[nodiscard]] static EffectiveProperties effectiveProperties(
        Real ri, Real rj, Real ei, Real ej, Real ni, Real nj, Real mi, Real mj);

    // Hertz normal force (undamped): F = (4/3) E* sqrt(R*) delta^(3/2)
    [[nodiscard]] static Real hertzNormalForce(const EffectiveProperties& props, Real overlap);

    // Normal damping coefficient from the restitution coefficient.
    [[nodiscard]] static Real normalDamping(Real restitution, Real m_eff, Real stiffness);

    [[nodiscard]] const ContactStats& stats() const { return stats_; }
    void clearHistory() { tangential_history_.clear(); }

private:
    // Pair key i<j for the tangential history (which does have memory).
    [[nodiscard]] static std::uint64_t pairKey(core::Index i, core::Index j) {
        const std::uint64_t a = static_cast<std::uint64_t>(i < j ? i : j);
        const std::uint64_t b = static_cast<std::uint64_t>(i < j ? j : i);
        return (a << 32) | b;
    }

    Real restitution_{0};
    Real friction_{0};
    bool tangential_enabled_{true};
    Real tangential_ratio_{1.0};
    ContactStats stats_{};
    // Tangential displacement accumulated per contact: it is what gives Mindlin
    // its hysteresis. Cleared when contacts close or open (documented limitation).
    std::unordered_map<std::uint64_t, core::Vec3> tangential_history_;
};

}  // namespace granimpact::physics
