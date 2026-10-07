// Verifies the integrator against analytical solutions. If this passes, the errors
// that show up in a simulation come from the physics or the collisions, not from
// the time-stepping scheme.
#include <gtest/gtest.h>

#include <cmath>
#include <random>

#include "granimpact/core/ParticleSystem.hpp"
#include "granimpact/physics/Integrator.hpp"

using namespace granimpact;

namespace {

core::ParticleSystem twoParticles() {
    core::ParticleSystem ps;
    ps.add({0, 0, 0.5}, 1e-3, 2650.0, core::ParticleKind::Bed);
    ps.add({0, 0, 0.4}, 1e-3, 2650.0, core::ParticleKind::Bed);
    ps.setMaterial(70e9, 0.3, 0.4, 0.5);
    return ps;
}

}  // namespace

TEST(Integrator, FreeFallMatchesAnalyticSolution) {
    auto ps = twoParticles();
    const core::Real dt = 1e-5;
    const core::Real gz = -9.81;
    physics::VelocityVerlet integrator(ps, dt, gz);

    const std::size_t steps = 500;
    for (std::size_t s = 0; s < steps; ++s) integrator.freeStep(ps);

    const core::Real t = dt * static_cast<core::Real>(steps);
    // z(t) = z0 - 0.5 g t^2 (starting from rest)
    const core::Real expected = 0.5 - 0.5 * 9.81 * t * t;
    EXPECT_NEAR(ps.position()[0].z, expected, 1e-12);

    // v(t) = -g t
    EXPECT_NEAR(ps.velocity()[0].z, -9.81 * t, 1e-12);
}

TEST(Integrator, VerletConservesEnergyInConstantField) {
    auto ps = twoParticles();
    const core::Real dt = 1e-5;
    physics::VelocityVerlet integrator(ps, dt, -9.81);

    const core::Real e0 = ps.kineticEnergy() + ps.gravitationalPotentialEnergy(-9.81);
    for (std::size_t s = 0; s < 2000; ++s) integrator.freeStep(ps);
    const core::Real e1 = ps.kineticEnergy() + ps.gravitationalPotentialEnergy(-9.81);

    // Under a constant force the scheme is exact: total energy must not drift.
    EXPECT_NEAR(e1, e0, std::abs(e0) * 1e-12 + 1e-20);
}

TEST(Integrator, KickDriftKickPhasesMatchFreeStep) {
    auto a = twoParticles();
    auto b = twoParticles();
    const core::Real dt = 1e-5;
    physics::VelocityVerlet integrator(a, dt, -9.81);

    // The kick-drift-kick path consumes the ACCUMULATED FORCES, so gravity has to be
    // applied to it the same way Simulation does at every step.
    for (std::size_t s = 0; s < 100; ++s) {
        b.clearForces();
        b.applyGravity(-9.81);
        integrator.firstHalfKick(b);
        integrator.drift(b);
        b.clearForces();
        b.applyGravity(-9.81);
        integrator.secondHalfKick(b);
    }
    for (std::size_t s = 0; s < 100; ++s) integrator.freeStep(a);

    EXPECT_NEAR(a.position()[0].z, b.position()[0].z, 1e-14);
    EXPECT_NEAR(a.velocity()[0].z, b.velocity()[0].z, 1e-14);
}

TEST(Integrator, PeriodicBoundariesWrapParticles) {
    auto ps = twoParticles();
    core::SimConfig cfg;
    cfg.domain = {0.1, 0.1, 0.1};
    cfg.boundaries.periodic_xy = true;
    cfg.boundaries.fixed_bottom = false;

    ps.position()[0] = {0.11, -0.01, 0.05};   // out to the right and below
    physics::applyBoundaries(ps, cfg.domain.box(), cfg);

    EXPECT_GE(ps.position()[0].x, 0.0);
    EXPECT_LT(ps.position()[0].x, 0.1);
    EXPECT_GE(ps.position()[0].y, 0.0);
    EXPECT_LT(ps.position()[0].y, 0.1);
}

TEST(Integrator, RigidWallsBounceWithRestitution) {
    auto ps = twoParticles();
    core::SimConfig cfg;
    cfg.domain = {0.1, 0.1, 0.1};
    cfg.boundaries.periodic_xy = false;
    cfg.boundaries.wall_restitution = 0.5;

    ps.position()[0] = {0.05, 0.05, -0.01};
    ps.velocity()[0] = {0, 0, -1.0};
    physics::applyBoundaries(ps, cfg.domain.box(), cfg);

    EXPECT_NEAR(ps.velocity()[0].z, 0.5, 1e-12);   // e = 0.5
    EXPECT_GE(ps.position()[0].z, 0.0);
}
