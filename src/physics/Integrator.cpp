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
    const auto& acc = particles.force();
    const auto& inv_m = particles.invMass();
    const Real half = Real(0.5) * dt_;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        vel[i] += acc[i] * (inv_m[i] * half);
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
    firstHalfKick(particles);   // symmetric: v += a dt/2 with the new forces
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
    auto& kinds = particles.kind();
    const auto& rad = particles.radius();
    const Real e_wall = config.boundaries.wall_restitution;
    const Real mu_wall = config.particles.friction;

    for (std::size_t i = 0; i < particles.size(); ++i) {
        const Real r = rad[i];
        // Rigid floor at z=0
        if (config.boundaries.fixed_bottom && pos[i].z - r < box.min.z) {
            pos[i].z = box.min.z + r;
            if (vel[i].z < 0) {
                vel[i].z = -e_wall * vel[i].z;
                vel[i].x *= (1.0 - mu_wall * 0.1);   // tangential loss to friction
                vel[i].y *= (1.0 - mu_wall * 0.1);
            }
        }
        // Side walls: periodic or rigid
        for (int axis = 0; axis < 2; ++axis) {
            Real& p = (axis == 0) ? pos[i].x : pos[i].y;
            Real& v = (axis == 0) ? vel[i].x : vel[i].y;
            const Real lo = (axis == 0) ? box.min.x : box.min.y;
            const Real hi = (axis == 0) ? box.max.x : box.max.y;
            const Real length = hi - lo;

            if (config.boundaries.periodic_xy) {
                if (p < lo) p += length;
                else if (p >= hi) p -= length;
            } else {
                if (p - r < lo) { p = lo + r; if (v < 0) v = -e_wall * v; }
                else if (p + r > hi) { p = hi - r; if (v > 0) v = -e_wall * v; }
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
