#include "granimpact/physics/ContactModel.hpp"

#include <algorithm>
#include <cmath>

namespace granimpact::physics {

using core::Index;
using core::kPi;
using core::ParticleKind;
using core::Real;
using core::Vec3;

namespace {

Real computeDampingFactor(Real restitution) {
    if (restitution >= 1.0) return 0.0;
    const Real e = std::max(restitution, Real(1e-6));
    const Real ln_e = std::log(e);
    return -2.0 * ln_e / std::sqrt(kPi * kPi + ln_e * ln_e);
}

}  // namespace

ContactModel::ContactModel(const core::SimConfig& config)
    : damping_factor_(computeDampingFactor(config.particles.restitution)),
      tangential_enabled_(config.contact.tangential),
      rolling_friction_enabled_(config.contact.rolling_friction),
      rolling_mu_(config.contact.rolling_friction_coeff),
      tangential_ratio_(config.contact.tangential_stiffness_ratio),
      cohesive_projectile_(config.contact.cohesive_projectile),
      yield_stress_Pa_(config.projectile.yield_kPa * Real(1000.0)),
      bed_cohesion_Pa_(config.contact.bed_cohesion_kPa * Real(1000.0)),
      magnetic_chi_(config.contact.magnetic_susceptibility),
      magnetic_b_T_(config.contact.magnetic_field_T),
      periodic_xy_(config.boundaries.periodic_xy),
      domain_size_{config.domain.x, config.domain.y, config.domain.z} {}

void ContactModel::clearHistory() {
    history_slots_.clear();
    has_cohesive_.clear();
    step_stamp_ = 0;
    dissipated_energy_ = 0.0;
    active_bonds_ = 0;
    broken_bonds_ = 0;
    bonds_initialized_ = false;
}

void ContactModel::ensureCapacity(std::size_t n_particles) {
    const std::size_t needed = n_particles * kMaxContactsPerParticle;
    if (history_slots_.size() < needed) {
        history_slots_.resize(needed);
    }
    if (has_cohesive_.size() < n_particles) {
        has_cohesive_.resize(n_particles, 0);
    }
}

ContactModel::ContactSlot* ContactModel::findOrInsertSlot(
    std::size_t i, Index j, bool create_if_missing) {
    ContactSlot* base = &history_slots_[i * kMaxContactsPerParticle];
    ContactSlot* free_slot = nullptr;

    for (std::size_t s = 0; s < kMaxContactsPerParticle; ++s) {
        ContactSlot& slot = base[s];
        if (slot.partner == j) {
            // If this contact separated in an earlier step and had no cohesive bond,
            // reset its accumulated tangential displacement upon re-contact.
            if (slot.cohesive == 0 && slot.last_step != step_stamp_
                && slot.last_step + 1 != step_stamp_) {
                slot.delta_t = Vec3{};
            }
            slot.last_step = step_stamp_;
            return &slot;
        }
        if (free_slot == nullptr && slot.cohesive == 0) {
            if (slot.partner < 0
                || (slot.last_step != step_stamp_ && slot.last_step + 1 != step_stamp_)) {
                free_slot = &slot;
            }
        }
    }

    if (!create_if_missing) return nullptr;
    if (free_slot == nullptr) {
        // Fallback: evict the last non-cohesive slot if all 16 slots were active.
        for (std::size_t s = kMaxContactsPerParticle; s > 0; --s) {
            if (base[s - 1].cohesive == 0) {
                free_slot = &base[s - 1];
                break;
            }
        }
        if (free_slot == nullptr) free_slot = &base[kMaxContactsPerParticle - 1];
    }

    free_slot->partner = j;
    free_slot->last_step = step_stamp_;
    free_slot->delta_t = Vec3{};
    free_slot->cohesive = 0;
    return free_slot;
}

void ContactModel::registerCohesiveBonds(
    const core::ParticleSystem& particles, Real max_gap_fraction) {
    if (!cohesive_projectile_ || yield_stress_Pa_ <= 0) return;

    const std::size_t n = particles.size();
    ensureCapacity(n);
    const auto& pos = particles.position();
    const auto& rad = particles.radius();
    const auto& kinds = particles.kind();

    // Collect indices of projectile particles to avoid O(N_bed^2) work.
    std::vector<std::size_t> proj_indices;
    for (std::size_t i = 0; i < n; ++i) {
        if (kinds[i] == ParticleKind::Projectile) {
            proj_indices.push_back(i);
        }
    }
    if (proj_indices.size() < 2) return;

    ++step_stamp_;
    for (std::size_t a = 0; a < proj_indices.size(); ++a) {
        const std::size_t i = proj_indices[a];
        for (std::size_t b = a + 1; b < proj_indices.size(); ++b) {
            const std::size_t j = proj_indices[b];
            const Vec3 d = pos[j] - pos[i];
            const Real r_sum = rad[i] + rad[j];
            const Real max_dist = r_sum * (Real(1.0) + max_gap_fraction);
            if (d.norm2() <= max_dist * max_dist) {
                ContactSlot* slot = findOrInsertSlot(i, static_cast<Index>(j), true);
                if (slot != nullptr && slot->cohesive == 0) {
                    slot->cohesive = 1;
                    has_cohesive_[i] = 1;
                    ++active_bonds_;
                }
            }
        }
    }
    bonds_initialized_ = true;
}

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
           * (overlap * std::sqrt(overlap));
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
    ensureCapacity(particles.size());
    if (!bonds_initialized_ && cohesive_projectile_ && yield_stress_Pa_ > 0) {
        // Auto-initialize cohesive bonds if projectile grains are already present.
        const auto& kinds = particles.kind();
        for (std::size_t k = 0; k < particles.size(); ++k) {
            if (kinds[k] == ParticleKind::Projectile) {
                registerCohesiveBonds(particles);
                break;
            }
        }
    }

    ++step_stamp_;
    stats_ = ContactStats{};

    auto& force = particles.force();
    auto& torque = particles.torque();
    auto& normal_force_sum = particles.normalForceSum();
    auto& coordination = particles.coordination();
    const auto& pos = particles.position();
    const auto& vel = particles.velocity();
    const auto& omega = particles.omega();
    const auto& rad = particles.radius();
    const auto& mass = particles.mass();
    const auto& young = particles.young();
    const auto& poisson = particles.poisson();
    const auto& friction = particles.friction();
    const auto& kinds = particles.kind();
    const bool has_magnetic = (magnetic_chi_ > 0 && magnetic_b_T_ > 0);
    // Precompute magnetic dipole prefactor: (3 * mu_0 / 4 pi) * (4 pi / mu_0 * chi/(chi+3) * B0)^2
    // = 12 * pi / mu_0 * (chi / (chi + 3) * B0)^2, with mu_0 = 4 pi * 1e-7 H/m -> 3e7 * (chi/(chi+3)*B0)^2.
    const Real chi_eff = has_magnetic ? (magnetic_chi_ / (magnetic_chi_ + Real(3.0))) : Real(0);
    const Real mag_prefactor = has_magnetic
                                   ? (Real(3.0e7) * (chi_eff * magnetic_b_T_) * (chi_eff * magnetic_b_T_))
                                   : Real(0);

    const Real half_lx = Real(0.5) * domain_size_.x;
    const Real half_ly = Real(0.5) * domain_size_.y;

    // Cell-list indices are `Index` (signed), with i_signed < j_signed guaranteed by CellList.
    cells.forEachPair([&](Index i_signed, Index j_signed) {
        const auto i = static_cast<std::size_t>(i_signed);
        const auto j = static_cast<std::size_t>(j_signed);
        ++stats_.candidates;

        Vec3 d = pos[j] - pos[i];
        if (periodic_xy_) {
            if (d.x > half_lx) d.x -= domain_size_.x;
            else if (d.x < -half_lx) d.x += domain_size_.x;
            if (d.y > half_ly) d.y -= domain_size_.y;
            else if (d.y < -half_ly) d.y += domain_size_.y;
        }

        const Real dist2 = d.norm2();
        const Real r_sum = rad[i] + rad[j];

        // Fast squared-distance rejection before std::sqrt for non-cohesive, non-magnetic pairs.
        ContactSlot* cohesive_slot = nullptr;
        if (dist2 >= r_sum * r_sum) {
            if (has_cohesive_[i] != 0) {
                cohesive_slot = findOrInsertSlot(i, j_signed, false);
            }
            if ((cohesive_slot == nullptr || cohesive_slot->cohesive == 0) && !has_magnetic) {
                return;
            }
        }
        if (dist2 <= 0) return;

        const Real dist = std::sqrt(dist2);
        const Real overlap = r_sum - dist;
        const Vec3 n = d * (Real(1.0) / dist);       // unit normal from i towards j

        // Magnetostatic dipole-dipole interaction (when enabled for magnetic sand / regolith)
        if (has_magnetic) {
            const Real eff_d = std::max(dist, Real(0.8) * r_sum);
            const Real r3_ij = (rad[i] * rad[i] * rad[i]) * (rad[j] * rad[j] * rad[j]);
            const Real d4 = (eff_d * eff_d) * (eff_d * eff_d);
            const Real f_mag_scale = mag_prefactor * r3_ij / d4;
            const Real nz = n.z;
            // F_mag on j = scale * ((1 - 5 nz^2) n + 2 nz z_hat)
            const Vec3 f_mag{f_mag_scale * (Real(1.0) - Real(5.0) * nz * nz) * n.x,
                             f_mag_scale * (Real(1.0) - Real(5.0) * nz * nz) * n.y,
                             f_mag_scale * ((Real(1.0) - Real(5.0) * nz * nz) * nz + Real(2.0) * nz)};
            force[i] -= f_mag;
            force[j] += f_mag;
            if (overlap <= 0 && (cohesive_slot == nullptr || cohesive_slot->cohesive == 0)) {
                return;
            }
        }

        const EffectiveProperties props = effectiveProperties(
            rad[i], rad[j], young[i], young[j], poisson[i], poisson[j], mass[i], mass[j]);

        const Real r_min_pair = std::min(rad[i], rad[j]);
        const Real bond_area = kPi * r_min_pair * r_min_pair;
        const Real max_bond_force = yield_stress_Pa_ * bond_area;

        // Check cohesive status if not already looked up.
        ContactSlot* slot = cohesive_slot;
        if (slot == nullptr && (tangential_enabled_ || has_cohesive_[i] != 0)) {
            slot = findOrInsertSlot(i, j_signed, tangential_enabled_);
        }
        const bool is_bonded = (slot != nullptr && slot->cohesive != 0);

        const Vec3 rel_vel_trans = vel[j] - vel[i];
        const Real v_n = rel_vel_trans.dot(n);

        Real f_n = 0;
        Real gamma_n = 0;
        if (overlap > 0) {
            ++stats_.touching;
            stats_.max_overlap = std::max(stats_.max_overlap, static_cast<double>(overlap));

            const Real sqrt_r_star = std::sqrt(props.r_star);
            const Real sqrt_overlap = std::sqrt(overlap);
            const Real k_hertz = Real(4.0 / 3.0) * props.e_star * sqrt_r_star;
            const Real k_linear = Real(2.0) * props.e_star * sqrt_r_star * sqrt_overlap;
            gamma_n = damping_factor_ * std::sqrt(props.m_eff * k_linear);

            const Real f_bed_coh = (bed_cohesion_Pa_ > 0 && kinds[i] == ParticleKind::Bed
                                    && kinds[j] == ParticleKind::Bed)
                                       ? (bed_cohesion_Pa_ * kPi * props.r_star * props.r_star)
                                       : Real(0);
            f_n = k_hertz * (overlap * sqrt_overlap) - gamma_n * v_n - f_bed_coh;
            if (!is_bonded && f_bed_coh <= 0) {
                f_n = std::max(f_n, Real(0));        // cohesionless contact never pulls
            } else if (!is_bonded && f_bed_coh > 0) {
                f_n = std::max(f_n, -f_bed_coh);     // capillary / van der Waals pull-off limit
            } else if (f_n < -max_bond_force) {
                // Tensile fracture of the cement bridge
                slot->cohesive = 0;
                if (active_bonds_ > 0) --active_bonds_;
                ++broken_bonds_;
                f_n = 0;
            }
        } else if (is_bonded) {
            // Cohesive bond under pure tension (overlap <= 0, gap = -overlap > 0)
            const Real k_bond = props.e_star * r_min_pair;
            const Real f_tensile = k_bond * overlap - damping_factor_ * std::sqrt(props.m_eff * k_bond) * v_n;
            if (-f_tensile > max_bond_force) {
                // Tensile rupture of the cement bond
                slot->cohesive = 0;
                if (active_bonds_ > 0) --active_bonds_;
                ++broken_bonds_;
                return;
            }
            f_n = f_tensile;
        } else {
            return;
        }

        const Real abs_fn = std::abs(f_n);
        stats_.max_normal_force = std::max(stats_.max_normal_force, static_cast<double>(abs_fn));
        normal_force_sum[i] += abs_fn;
        normal_force_sum[j] += abs_fn;
        ++coordination[i];
        ++coordination[j];

        // --- Tangential: 6-DOF Mindlin-Deresiewicz with rotational slip + Coulomb / cohesive cap
        Vec3 f_t{};
        if (tangential_enabled_ && slot != nullptr) {
            // Surface relative tangential velocity including grain rotation:
            // v_t = (v_j - v_i - v_n n) - (r_i omega_i + r_j omega_j) x n
            const Vec3 v_rot = (omega[i] * rad[i] + omega[j] * rad[j]).cross(n);
            const Vec3 v_t = (rel_vel_trans - n * v_n) - v_rot;

            // Project previous tangential displacement onto current tangent plane and advance
            Vec3& delta_t = slot->delta_t;
            delta_t -= n * delta_t.dot(n);
            delta_t += v_t * dt;

            const Real eff_overlap = std::max(overlap, Real(1e-6) * props.r_star);
            const Real k_s = Real(8.0) * props.g_star * std::sqrt(props.r_star * eff_overlap)
                             * tangential_ratio_;
            const Real gamma_t = damping_factor_ * std::sqrt(props.m_eff * k_s);
            const Vec3 f_spring = delta_t * (-k_s);
            const Vec3 f_trial = f_spring - v_t * gamma_t;

            const Real mu = std::min(friction[i], friction[j]);
            const Real f_t_norm = f_trial.norm();
            const Real f_bed_coh = (bed_cohesion_Pa_ > 0 && kinds[i] == ParticleKind::Bed
                                    && kinds[j] == ParticleKind::Bed)
                                       ? (bed_cohesion_Pa_ * kPi * props.r_star * props.r_star)
                                       : Real(0);
            Real limit = mu * std::max(f_n + f_bed_coh, Real(0)) + Real(0.5) * f_bed_coh;
            if (slot->cohesive != 0) {
                limit += max_bond_force;
            }

            if (f_t_norm > limit && f_t_norm > 0) {
                if (slot->cohesive != 0) {
                    // Shear comminution of the cement bridge (Regime 2: P_0 > sigma_y)
                    slot->cohesive = 0;
                    if (active_bonds_ > 0) --active_bonds_;
                    ++broken_bonds_;
                    limit = mu * std::max(f_n, Real(0));
                }
                f_t = f_trial * (limit / f_t_norm);
                if (k_s > 0) {
                    delta_t = f_t * (Real(-1.0) / k_s);
                }
                ++stats_.sliding;
                dissipated_energy_ += static_cast<double>(limit * v_t.norm() * dt);
            } else {
                f_t = f_trial;
                dissipated_energy_ += static_cast<double>(gamma_t * v_t.norm2() * dt);
            }
        }

        if (gamma_n > 0) {
            dissipated_energy_ += static_cast<double>(gamma_n * v_n * v_n * dt);
        }

        // `n` goes from i to j: repulsion on i is -n f_n and on j is +n f_n.
        const Vec3 total = n * f_n + f_t;
        force[i] -= total;
        force[j] += total;

        // Tangential torque: tau = r_c x F_t. Both torques carry the same sign.
        const Vec3 lever_factor = n.cross(f_t);
        torque[i] -= lever_factor * rad[i];
        torque[j] -= lever_factor * rad[j];

        // --- Rolling resistance torque: Ai et al. (2011) / Luding (2008)
        // Crucial for reproducing the experimental angle of repose (40-44 deg) and rim stability.
        if (rolling_friction_enabled_ && f_n > 0) {
            const Vec3 omega_rel = omega[i] - omega[j];
            const Vec3 omega_roll = omega_rel - n * omega_rel.dot(n);
            const Real omega_roll_norm = omega_roll.norm();
            if (omega_roll_norm > Real(1e-12)) {
                const Real tau_r_max = rolling_mu_ * props.r_star * f_n;
                const Real eff_overlap = std::max(overlap, Real(1e-6) * props.r_star);
                const Real k_r = Real(16.0) * props.g_star * std::sqrt(props.r_star * eff_overlap)
                                 * props.r_star * props.r_star;
                const Real tau_r_mag = std::min(tau_r_max, k_r * omega_roll_norm * dt);
                const Vec3 tau_r = omega_roll * (tau_r_mag / omega_roll_norm);
                torque[i] -= tau_r;
                torque[j] += tau_r;
                dissipated_energy_ += static_cast<double>(tau_r_mag * omega_roll_norm * dt);
            }
        }
    });

    stats_.active_bonds = active_bonds_;
    stats_.broken_bonds = broken_bonds_;
    stats_.dissipated_energy = dissipated_energy_;
    return stats_;
}

}  // namespace granimpact::physics
