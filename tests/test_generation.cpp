// Generation fixes the initial conditions. What is required here is what
// holds up all the physics afterwards: no overlaps, inside the box, with the
// requested radius distribution and reproducible with the same seed.
#include <gtest/gtest.h>

#include <cmath>
#include <random>

#include "granimpact/core/SimConfig.hpp"
#include "granimpact/generation/Generators.hpp"

using namespace granimpact;

namespace {

core::SimConfig testConfig() {
    core::SimConfig cfg;
    cfg.domain = {0.06, 0.06, 0.05};
    cfg.particles.r_min = 1.0e-3;
    cfg.particles.r_max = 1.5e-3;
    cfg.particles.packing_fraction = 0.524;
    cfg.projectile.radius = 8.0e-3;
    cfg.projectile.start_z = 0.04;
    cfg.simulation.seed = 20240517;
    return cfg;
}

}  // namespace

TEST(Generators, BedIsSeededInsideRegionWithoutOverlaps) {
    auto cfg = testConfig();
    core::ParticleSystem ps;
    std::mt19937_64 rng(cfg.simulation.seed);

    const auto info = generation::generateBed(ps, cfg, rng);
    const core::Box box = cfg.domain.box();

    ASSERT_GT(info.count, 1000u);
    EXPECT_EQ(ps.size(), info.count);

    for (std::size_t i = 0; i < ps.size(); ++i) {
        const auto& p = ps.position()[i];
        const core::Real r = ps.radius()[i];
        EXPECT_GE(ps.radius()[i], cfg.particles.r_min - 1e-12);
        EXPECT_LE(ps.radius()[i], cfg.particles.r_max + 1e-12);
        EXPECT_GE(p.x - r, box.min.x - 1e-9);
        EXPECT_LE(p.x + r, box.max.x + 1e-9);
        EXPECT_GE(p.y - r, box.min.y - 1e-9);
        EXPECT_LE(p.y + r, box.max.y + 1e-9);
        EXPECT_GE(p.z - r, box.min.z - 1e-9);
        EXPECT_LE(p.z + r, info.region_height + 1e-9);
    }
}

TEST(Generators, AnyTwoGrainsAreSeparated) {
    auto cfg = testConfig();
    core::ParticleSystem ps;
    std::mt19937_64 rng(4242);
    generation::generateBed(ps, cfg, rng);

    std::size_t overlaps = 0;
    for (std::size_t i = 0; i < ps.size(); ++i) {
        for (std::size_t j = i + 1; j < ps.size(); ++j) {
            const core::Real dist = (ps.position()[i] - ps.position()[j]).norm();
            if (dist < ps.radius()[i] + ps.radius()[j] - 1e-12) ++overlaps;
        }
    }
    EXPECT_EQ(overlaps, 0u);
}

TEST(Generators, NarrowRadiusSpreadReachesTargetDensity) {
    // Generator contract: the nested lattice is expanded to land EXACTLY on the
    // target density, and that is only possible if the dense lattice is denser than
    // the target. With a narrow radius range (1.10:1 here) it is, and the seeded
    // density must match the config target.
    auto cfg = testConfig();
    cfg.particles.r_min = 2.0e-3;
    cfg.particles.r_max = 2.2e-3;
    core::ParticleSystem ps;
    std::mt19937_64 rng(7);
    const auto info = generation::generateBed(ps, cfg, rng);

    EXPECT_NEAR(info.seeded_phi, cfg.particles.packing_fraction, 0.02);
    // And the lattice is still a lattice: the horizontal spacing covers the maximum
    // diameter and the vertical one is smaller (the layers are nested, not stacked).
    EXPECT_GT(info.lattice_spacing, 2 * cfg.particles.r_max);
    EXPECT_LT(info.layer_spacing, info.lattice_spacing);
    EXPECT_GT(info.count, 100u);
}

TEST(Generators, WideRadiusSpreadDoesNotCompressTheLattice) {
    // With a wide range (1.5:1, as in the smoke test) the dense lattice falls below
    // the target. The design decision is NOT to compress - compressing would create
    // overlaps - and to report the density obtained: the bed ends up looser and
    // settling is what takes it to its physical density.
    auto cfg = testConfig();
    core::ParticleSystem ps;
    std::mt19937_64 rng(7);
    const auto info = generation::generateBed(ps, cfg, rng);

    EXPECT_LT(info.seeded_phi, cfg.particles.packing_fraction);
    EXPECT_GT(info.lattice_spacing, 2 * cfg.particles.r_max);
    // Even so it has to be a recognisable bed, not a handful of grains.
    EXPECT_GT(info.seeded_phi, 0.25);
    EXPECT_GT(info.count, 100u);
}

TEST(Generators, BedParticlesCarryConfigMaterial) {
    auto cfg = testConfig();
    cfg.particles.young = 1.0e7;
    cfg.particles.friction = 0.37;
    core::ParticleSystem ps;
    std::mt19937_64 rng(11);
    generation::generateBed(ps, cfg, rng);

    ASSERT_GT(ps.size(), 0u);
    EXPECT_NEAR(ps.young()[0], 1.0e7, 1e-3);
    EXPECT_NEAR(ps.friction()[0], 0.37, 1e-12);
    EXPECT_EQ(ps.kind()[0], core::ParticleKind::Bed);
    // Mass comes from density and radius, not from a fixed value.
    const core::Real r = ps.radius()[0];
    const core::Real expected = cfg.particles.rho_grain * 4.0 / 3.0 * core::kPi * r * r * r;
    EXPECT_NEAR(ps.mass()[0], expected, expected * 1e-12);
}

TEST(Generators, ProjectileIsAGrainSphereWithRequestedSpeed) {
    auto cfg = testConfig();
    cfg.projectile.speed = 2.0;
    cfg.projectile.angle_deg = 90.0;

    core::ParticleSystem ps;
    std::mt19937_64 rng(3);
    const auto count = generation::generateProjectile(ps, cfg, rng);

    ASSERT_GT(count, 20u);
    EXPECT_EQ(ps.size(), count);

    const core::Vec3 center{cfg.domain.x / 2, cfg.domain.y / 2, cfg.projectile.start_z};
    core::Real max_dist = 0;
    for (std::size_t i = 0; i < ps.size(); ++i) {
        const core::Real dist = (ps.position()[i] - center).norm();
        // All grains fall inside the projectile sphere...
        EXPECT_LE(dist + ps.radius()[i], cfg.projectile.radius + 1e-12);
        max_dist = std::max(max_dist, dist);
        EXPECT_EQ(ps.kind()[i], core::ParticleKind::Projectile);
        // ...and move with the requested velocity (90 degrees = downwards).
        EXPECT_NEAR(ps.velocity()[i].z, -2.0, 1e-12);
        EXPECT_NEAR(ps.velocity()[i].x, 0.0, 1e-12);
        EXPECT_NEAR(ps.velocity()[i].y, 0.0, 1e-12);
    }
    // The sphere is filled, not just hollow.
    EXPECT_GT(max_dist, 0.5 * cfg.projectile.radius);
}

TEST(Generators, ObliqueAngleSplitsSpeedInXY) {
    auto cfg = testConfig();
    cfg.projectile.speed = 2.0;
    cfg.projectile.angle_deg = 45.0;

    core::ParticleSystem ps;
    std::mt19937_64 rng(5);
    generation::generateProjectile(ps, cfg, rng);

    ASSERT_GT(ps.size(), 0u);
    EXPECT_NEAR(ps.velocity()[0].x, 2.0 * std::cos(core::kPi / 4), 1e-12);
    EXPECT_NEAR(ps.velocity()[0].z, -2.0 * std::sin(core::kPi / 4), 1e-12);
}

TEST(Generators, TwoGenerationsWithSameSeedAreIdentical) {
    auto cfg = testConfig();
    core::ParticleSystem a, b;
    std::mt19937_64 ra(cfg.simulation.seed), rb(cfg.simulation.seed);
    generation::generateBed(a, cfg, ra);
    generation::generateBed(b, cfg, rb);
    ASSERT_EQ(a.size(), b.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        EXPECT_NEAR(a.position()[i].x, b.position()[i].x, 1e-15);
        EXPECT_NEAR(a.radius()[i], b.radius()[i], 1e-15);
    }
}

TEST(Generators, MeasuredPackingFractionUsesSolidVolume) {
    core::ParticleSystem ps;
    ps.add({0, 0, 0}, 1e-3, 2650.0, core::ParticleKind::Bed);
    ps.add({0, 0, 0.01}, 1e-3, 2650.0, core::ParticleKind::Bed);
    const core::Real expected = 2 * (4.0 / 3.0) * core::kPi * 1e-9;
    EXPECT_NEAR(ps.solidsVolume(), expected, expected * 1e-12);
    // Volume such that phi = 0.5 -> solid / volume = 0.5.
    EXPECT_NEAR(generation::measuredPacking(ps, expected / 0.5), 0.5, 1e-12);
}
