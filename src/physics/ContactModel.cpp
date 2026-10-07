#include "granimpact/physics/ContactModel.hpp"

#include <algorithm>
#include <cmath>

namespace granimpact::physics {

using core::Index;
using core::kPi;
using core::Real;
using core::Vec3;

ContactModel::ContactModel(const core::SimConfig& config)
    : restitution_(config.particles.restitution),
      tangential_enabled_(config.contact.tangential),
      tangential_ratio_(config.contact.tangential_stiffness_ratio) {}

EffectiveProperties ContactModel::effectiveProperties(
    Real ri, Real rj, Real ei, Real ej, Real ni, Real nj, Real mi, Real mj) {
    EffectiveProperties p;
    // 1/E* = (1-nu_i^2)/E_i + (1-nu_j^2)/E_j
    const Real inv_e = (1.0 - ni * ni) / ei + (1.0 - nj * nj) / ej;
    p.e_star = 1.0 / inv_e;
    // Shear moduli and 1/G* = (2-nu_i)/G_i + (2-nu_j)/G_j
    const Real gi = ei / (2.0 * (1.0 + ni));
    const Real gj = ej / (2.0 * (1.0 + nj));
    p.g_star = 1.0 / ((2.0 - ni) / gi + (2.0 - nj) / gj);
    p.r_star = (ri * rj) / (ri + rj);
    p.m_eff = (mi * mj) / (mi + mj);
    return p;
}

Real ContactModel::hertzNormalForce(const EffectiveProperties& props, Real overlap) {
    if (overlap <= 0) return 0;
    return (Real(4.0 / 3.0) * props.e_star) * std::sqrt(props.r_star)
           * std::pow(overlap, Real(1.5));
}

Real ContactModel::normalDamping(Real restitution, Real m_eff, Real stiffness) {
    // gamma = -2 ln(e) sqrt(m_eff * k) / sqrt(pi^2 + ln^2(e))   (Tsuji / Di Renzo)
    if (restitution >= 1.0) return 0.0;              // perfectly elastic collision
    const Real e = std::max(restitution, Real(1e-6));
    const Real ln_e = std::log(e);
    return -2.0 * ln_e * std::sqrt(m_eff * stiffness) / std::sqrt(kPi * kPi + ln_e * ln_e);
}

ContactStats ContactModel::computeForces(core::ParticleSystem& particles, const CellList& cells,
                                        Real dt) {
    stats_ = ContactStats{};

    auto& force = particles.force();
    auto& torque = particles.torque();
    const auto& pos = particles.position();
    const auto& vel = particles.velocity();
    const auto& rad = particles.radius();
    const auto& mass = particles.mass();
    const auto& young = particles.young();
    const auto& poisson = particles.poisson();
    const auto& friction = particles.friction();

    // Cell-list indices are `Index` (signed); the containers are indexed with
    // std::size_t. The conversion is done once per pair, explicitly, so that the
    // loop stays clean and free of implicit sign conversions.
    cells.forEachPair([&](Index i_signed, Index j_signed) {
        const auto i = static_cast<std::size_t>(i_signed);
        const auto j = static_cast<std::size_t>(j_signed);
        ++stats_.candidates;
        const Vec3 d = pos[j] - pos[i];
        const Real dist = d.norm();
        const Real r_sum = rad[i] + rad[j];
        const Real overlap = r_sum - dist;
        if (overlap <= 0 || dist <= 0) return;       // no contact

        ++stats_.touching;
        stats_.max_overlap = std::max(stats_.max_overlap, static_cast<double>(overlap));

        const Vec3 n = d / dist;                     // normal from i towards j
        const EffectiveProperties props = effectiveProperties(
            rad[i], rad[j], young[i], young[j], poisson[i], poisson[j], mass[i], mass[j]);

        // --- Normal: Hertz + damping that depends on the restitution coefficient
        const Real k_hertz = Real(4.0 / 3.0) * props.e_star * std::sqrt(props.r_star);
        const Real k_linear = 2.0 * props.e_star * std::sqrt(props.r_star * overlap);
        const Real gamma = normalDamping(restitution_, props.m_eff, k_linear);

        const Vec3 rel_vel = vel[j] - vel[i];
        const Real v_n = rel_vel.dot(n);
        Real f_n = k_hertz * std::pow(overlap, Real(1.5)) - gamma * v_n;
        f_n = std::max(f_n, Real(0));                // the contact never pulls, it only pushes
        stats_.max_normal_force = std::max(stats_.max_normal_force, static_cast<double>(f_n));

        // --- Tangential: Mindlin with history + Coulomb cap
        Vec3 f_t{};
        if (tangential_enabled_) {
            const Vec3 v_t = rel_vel - n * v_n;
            const std::uint64_t key = pairKey(i_signed, j_signed);
            Vec3& delta_t = tangential_history_[key];
            delta_t += v_t * dt;

            const Real k_s = 8.0 * props.g_star * std::sqrt(props.r_star * overlap)
                             * tangential_ratio_;
            Vec3 f_spring = delta_t * (-k_s);
            const Real mu = std::min(friction[i], friction[j]);
            const Real f_t_norm = f_spring.norm();
            const Real limit = mu * f_n;
            if (f_t_norm > limit && f_t_norm > 0) {
                // Sliding: the force is clipped and the history is corrected
                // so that it does not grow artificially.
                f_t = f_spring * (limit / f_t_norm);
                delta_t = delta_t * (limit / f_t_norm);
                ++stats_.sliding;
            } else {
                f_t = f_spring;
            }
        }

        // `n` goes from i to j, so the repulsion on i is -n f_n and on j is
        // +n f_n: the contact pushes, never pulls. This is the classic sign error.
        const Vec3 total = n * f_n + f_t;
        force[i] -= total;
        force[j] += total;

        // tau = r x F_t, with r from each grain's centre to the contact point.
        // Both torques carry the same sign: opposite forces applied at the
        // same point with opposite lever arms.
        const Vec3 lever_factor = n.cross(f_t);
        torque[i] -= lever_factor * rad[i];
        torque[j] -= lever_factor * rad[j];
    });

    return stats_;
}

}  // namespace granimpact::physics
