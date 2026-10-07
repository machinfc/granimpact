#pragma once
// Velocity Verlet integrator in separable steps (kick-drift-kick).
//
// Why in two halves: between the drift and the second half the forces must be
// recomputed. That is exactly where the OpenMP and CUDA kernels get inserted
// later without touching the integrator or the physics.
#include "granimpact/core/ParticleSystem.hpp"
#include "granimpact/core/SimConfig.hpp"

namespace granimpact::physics {

using core::Real;

class VelocityVerlet {
public:
    VelocityVerlet(const core::ParticleSystem& particles, core::Real dt, core::Real gz);

    // v += a * dt/2 (first half of the kick)
    void firstHalfKick(core::ParticleSystem& particles) const;

    // x += v * dt     (drift)
    void drift(core::ParticleSystem& particles) const;

    // v += a * dt/2 (second half of the kick, with the forces already recomputed)
    void secondHalfKick(core::ParticleSystem& particles) const;

    [[nodiscard]] core::Real dt() const { return dt_; }

    // Full step without contacts (free fall): useful for tests and for settling.
    void freeStep(core::ParticleSystem& particles) const;

private:
    core::Real dt_{0};
    core::Real gz_{0};
};

// Boundary conditions: rigid walls with restitution and, optionally,
// periodicity in x/y (which avoids edge effects in small boxes).
void applyBoundaries(core::ParticleSystem& particles, const core::Box& box,
                     const core::SimConfig& config);

}  // namespace granimpact::physics
