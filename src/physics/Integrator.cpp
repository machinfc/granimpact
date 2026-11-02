#include "granimpact/physics/Integrator.hpp"

#include <algorithm>
#include <cmath>

namespace granimpact::physics {

using core::Index;
using core::Real;
using core::Vec3;

VelocityVerlet::VelocityVerlet(const core::ParticleSystem& particles, Real dt, Real gz)
    : dt_(dt), gz_(gz) {
    (void)particles;
}

void VelocityVerlet::firstHalfKick(core::ParticleSystem& particles) const {
    auto& vel = particles.velocity();
    auto& omega = particles.omega();
    const auto& force = particles.force();
    const auto& torque = particles.torque();
    const auto& inv_m = particles.invMass();
    const auto& inv_i = particles.invInertia();
    const Real half = Real(0.5) * dt_;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        vel[i] += force[i] * (inv_m[i] * half);
        omega[i] += torque[i] * (inv_i[i] * half);
    }
}

void VelocityVerlet::drift(core::ParticleSystem& particles) const {
    auto& pos = particles.position();
    const auto& vel = particles.velocity();
    for (std::size_t i = 0; i < particles.size(); ++i) {
        pos[i] += vel[i] * dt_;
    }
}

void VelocityVerlet::secondHalfKick(core::ParticleSystem& particles) const {
    firstHalfKick(particles);   // symmetric: v += a dt/2, omega += alpha dt/2 with new forces/torques
}

void VelocityVerlet::freeStep(core::ParticleSystem& particles) const {
    auto& pos = particles.position();
    auto& vel = particles.velocity();
    const auto& inv_m = particles.invMass();
    // For a constant force this sequence is EXACT (not approximate): that is the
    // reason for using Verlet instead of Euler. test_integrator.cpp verifies it
    // against the analytical solution.
    for (std::size_t i = 0; i < particles.size(); ++i) {
        const Real a = gz_;                     // gravity only: a = g (mass cancels)
        pos[i].z += vel[i].z * dt_ + Real(0.5) * a * dt_ * dt_;
        vel[i].z += a * dt_;
        (void)inv_m;
    }
}

void applyBoundaries(core::ParticleSystem& particles, const core::Box& box,
                     const core::SimConfig& config) {
    auto& pos = particles.position();
    auto& vel = particles.velocity();
    auto& omega = particles.omega();
    auto& kinds = particles.kind();
    const auto& rad = particles.radius();
    const Real e_wall = config.boundaries.wall_restitution;
    const Real e_side = e_wall * (Real(1.0) - config.boundaries.wall_damping);
    const Real mu_wall = config.boundaries.wall_friction;
    const Real tang_decay = std::clamp(Real(1.0) - mu_wall * Real(0.15), Real(0.0), Real(1.0));
    const Real rot_decay = std::clamp(Real(1.0) - mu_wall * Real(0.10), Real(0.0), Real(1.0));

    for (std::size_t i = 0; i < particles.size(); ++i) {
        const Real r = rad[i];
        // Rigid floor at z = box.min.z
        if (config.boundaries.fixed_bottom && pos[i].z - r < box.min.z) {
            pos[i].z = box.min.z + r;
            if (vel[i].z < 0) {
                vel[i].z = -e_wall * vel[i].z;
                vel[i].x *= tang_decay;
                vel[i].y *= tang_decay;
                omega[i] *= rot_decay;
            }
        }
        // Side walls: periodic or rigid with optional acoustic damping and Coulomb friction
        for (int axis = 0; axis < 2; ++axis) {
            Real& p = (axis == 0) ? pos[i].x : pos[i].y;
            Real& v = (axis == 0) ? vel[i].x : vel[i].y;
            Real& vt1 = (axis == 0) ? vel[i].y : vel[i].x;
            const Real lo = (axis == 0) ? box.min.x : box.min.y;
            const Real hi = (axis == 0) ? box.max.x : box.max.y;
            const Real length = hi - lo;

            if (config.boundaries.periodic_xy) {
                if (p < lo) p += length;
                else if (p >= hi) p -= length;
            } else {
                if (p - r < lo) {
                    p = lo + r;
                    if (v < 0) {
                        v = -e_side * v;
                        vt1 *= tang_decay;
                        vel[i].z *= tang_decay;
                        omega[i] *= rot_decay;
                    }
                } else if (p + r > hi) {
                    p = hi - r;
                    if (v > 0) {
                        v = -e_side * v;
                        vt1 *= tang_decay;
                        vel[i].z *= tang_decay;
                        omega[i] *= rot_decay;
                    }
                }
            }
        }
        // No ceiling: ejecta may escape the domain (correct physics for the
        // experiment). It is flagged as Ejecta for the later analysis.
        if (pos[i].z > box.max.z && kinds[i] == core::ParticleKind::Bed) {
            kinds[i] = core::ParticleKind::Ejecta;
        }
    }
}

}  // namespace granimpact::physics
