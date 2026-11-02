// The contact tests are the project's safety net: if Hertz, Mindlin or
// Coulomb break, everything computed afterwards (and the thesis figures)
// stops being valid.
#include <gtest/gtest.h>

#include <cmath>
#include <random>

#include "granimpact/core/ParticleSystem.hpp"
#include "granimpact/core/SimConfig.hpp"
#include "granimpact/physics/CellList.hpp"
#include "granimpact/physics/ContactModel.hpp"
#include "granimpact/physics/Integrator.hpp"

using namespace granimpact;
using physics::ContactModel;
using physics::EffectiveProperties;

namespace {

// Two identical grains whose centres are 2r - overlap apart: the overlap is `overlap`.
core::ParticleSystem overlappingPair(core::Real overlap, core::Real v_rel = 0) {
    core::ParticleSystem ps;
    const core::Real r = 1e-3;
    ps.add({0, 0, r}, r, 2650.0, core::ParticleKind::Bed);                    // below
    ps.add({0, 0, r + 2 * r - overlap}, r, 2650.0, core::ParticleKind::Bed); // above
    ps.setMaterial(70e9, 0.3, 0.4, 0.5);
    ps.velocity()[1] = {0, 0, v_rel};
    return ps;
}

}  // namespace

TEST(ContactModel, EffectivePropertiesEqualForIdenticalSpheres) {
    const auto p = ContactModel::effectiveProperties(1e-3, 1e-3, 70e9, 70e9, 0.3, 0.3,
                                                     2.0e-5, 2.0e-5);
    // E* = E / (2 (1 - nu^2))
    EXPECT_NEAR(p.e_star, 70e9 / (2.0 * (1.0 - 0.09)), 1e6);
    // R* = r/2
    EXPECT_NEAR(p.r_star, 5e-4, 1e-18);
    // m* = m/2
    EXPECT_NEAR(p.m_eff, 1.0e-5, 1e-18);
    // For identical materials E* = E/(2(1-nu^2)) = 0.53 E and G* ~= 0.22 E*, which
    // is why Mindlin's tangential stiffness is lower than the normal one.
    EXPECT_LT(p.g_star, p.e_star);
    EXPECT_GT(p.g_star, 0.15 * p.e_star);
    EXPECT_LT(p.g_star, 0.30 * p.e_star);
}

TEST(ContactModel, NormalForceFollowsHertzLaw) {
    const auto p = ContactModel::effectiveProperties(1e-3, 1e-3, 70e9, 70e9, 0.3, 0.3, 1e-5, 1e-5);
    const core::Real delta = 1e-7;
    const core::Real f = ContactModel::hertzNormalForce(p, delta);
    const core::Real expected = (4.0 / 3.0) * p.e_star * std::sqrt(p.r_star) * std::pow(delta, 1.5);
    EXPECT_NEAR(f, expected, expected * 1e-12);
    // F(4 delta) = 8 F(delta): that nonlinearity is precisely the point of the Hertz model.
    EXPECT_NEAR(ContactModel::hertzNormalForce(p, 4 * delta), 8 * f, 8 * f * 1e-9);
    // No penetration, no force.
    EXPECT_EQ(ContactModel::hertzNormalForce(p, 0.0), 0.0);
}

TEST(ContactModel, DampingDecreasesWithRestitutionCoefficient) {
    const auto p = ContactModel::effectiveProperties(1e-3, 1e-3, 70e9, 70e9, 0.3, 0.3, 1e-5, 1e-5);
    const core::Real k = 1e5;
    const core::Real g_elastic = ContactModel::normalDamping(1.0, p.m_eff, k);
    const core::Real g_mid = ContactModel::normalDamping(0.5, p.m_eff, k);
    const core::Real g_inelastico = ContactModel::normalDamping(0.1, p.m_eff, k);

    EXPECT_EQ(g_elastic, 0.0);              // e = 1 -> no dissipation
    EXPECT_GT(g_mid, 0.0);
    EXPECT_GT(g_inelastico, g_mid);        // lower restitution -> more dissipation
}

TEST(ContactModel, TwoOverlappingParticlesRepelWithOppositeForces) {
    auto ps = overlappingPair(1e-7);
    core::SimConfig cfg;
    cfg.particles.r_min = 1e-3;
    cfg.particles.r_max = 1e-3;
    cfg.particles.young = 70e9;
    cfg.particles.poisson = 0.3;
    cfg.particles.friction = 0.4;
    cfg.particles.restitution = 0.5;

    physics::ContactModel model(cfg);
    physics::CellList cells;
    cells.build(ps, cfg.domain.box(), physics::CellList::suggestCellSize(1e-3));
    const auto stats = model.computeForces(ps, cells, 1e-7);

    ASSERT_EQ(stats.touching, 1u);
    const core::Vec3 f0 = ps.force()[0];
    const core::Vec3 f1 = ps.force()[1];

    EXPECT_LT(f0.z, 0.0);          // the lower particle pushed towards -z
    EXPECT_GT(f1.z, 0.0);          // the upper one towards +z
    EXPECT_NEAR(f0.x + f1.x, 0.0, 1e-18);
    EXPECT_NEAR(f0.y + f1.y, 0.0, 1e-18);
    EXPECT_NEAR(f0.z + f1.z, 0.0, 1e-18);   // Newton's third law

    // The force must be close to pure Hertz, because the relative velocities are zero.
    const auto props = ContactModel::effectiveProperties(1e-3, 1e-3, 70e9, 70e9, 0.3, 0.3,
                                                        ps.mass()[0], ps.mass()[1]);
    EXPECT_NEAR(f1.z, ContactModel::hertzNormalForce(props, 1e-7),
                ContactModel::hertzNormalForce(props, 1e-7) * 1e-9);
}

TEST(ContactModel, CollisionDissipatesEnergyAndRespectsRestitution) {
    // Head-on collision of two identical particles with e = 0.8: the relative velocity
    // must land around e * v0 (the damping model is approximate, hence the
    // tolerance).
    const core::Real r = 5e-4;
    const core::Real v0 = 0.5;
    core::SimConfig cfg;
    cfg.gravity.gz = 0.0;
    cfg.particles.r_min = r;
    cfg.particles.r_max = r;
    cfg.particles.young = 70e9;
    cfg.particles.poisson = 0.3;
    cfg.particles.friction = 0.0;          // no tangential: normal only
    cfg.particles.restitution = 0.8;
    cfg.contact.tangential = false;
    cfg.boundaries.fixed_bottom = false;
    cfg.boundaries.periodic_xy = false;
    cfg.domain = {0.01, 0.01, 0.01};

    core::ParticleSystem ps;
    ps.add({0.005, 0.005, 0.005 - r + 1e-9}, r, 2650.0, core::ParticleKind::Bed);
    ps.add({0.005, 0.005, 0.005 + r - 1e-9}, r, 2650.0, core::ParticleKind::Bed);
    ps.setMaterial(70e9, 0.3, 0.0, 0.8);
    ps.velocity()[0] = {0, 0, v0};          // they approach each other
    ps.velocity()[1] = {0, 0, -v0};

    physics::ContactModel model(cfg);
    physics::VelocityVerlet integrator(ps, 1e-8, 0.0);

    const core::Real e0 = ps.kineticEnergy();
    // The contact lasts ~1e-7 s (tens of steps); the rest of the time the two
    // particles do not interact. 400 steps are plenty to measure the bounce.
    for (int s = 0; s < 400; ++s) {
        integrator.firstHalfKick(ps);
        integrator.drift(ps);
        ps.clearForces();
        physics::CellList cells;
        cells.build(ps, cfg.domain.box(), physics::CellList::suggestCellSize(r));
        model.computeForces(ps, cells, 1e-8);
        integrator.secondHalfKick(ps);
    }

    // v_rel = v_below - v_above: they were approaching (v_rel > 0) and must
    // separate now (v_rel < 0), with |v_rel| ~ e * 2v0.
    const core::Real v_rel = ps.velocity()[0].z - ps.velocity()[1].z;
    EXPECT_LT(v_rel, 0.0);
    EXPECT_NEAR(std::abs(v_rel), 0.8 * 2 * v0, 0.10 * 2 * v0);
    EXPECT_LT(ps.kineticEnergy(), e0);                       // and it dissipates energy
}

TEST(ContactModel, TangentialForceIsLimitedByCoulomb) {
    // With friction 0 there can be no tangential force no matter how large the
    // tangential slip is: the limit mu*Fn cancels it.
    auto ps = overlappingPair(1e-7);
    core::SimConfig cfg;
    cfg.particles.r_min = 1e-3;
    cfg.particles.r_max = 1e-3;
    cfg.particles.young = 70e9;
    cfg.particles.poisson = 0.3;
    cfg.particles.friction = 0.0;
    cfg.particles.restitution = 0.5;
    cfg.contact.tangential = true;

    // Friction is a property of each particle (the model uses the ParticleSystem
    // materials, not the config ones), so it has to be set to 0 here.
    ps.setMaterial(70e9, 0.3, 0.0, 0.5);
    ps.velocity()[0] = {1e-3, 0, 0};
    ps.velocity()[1] = {-1e-3, 0, 0};

    physics::ContactModel model(cfg);
    physics::CellList cells;
    cells.build(ps, cfg.domain.box(), physics::CellList::suggestCellSize(1e-3));
    model.computeForces(ps, cells, 1e-6);

    EXPECT_NEAR(ps.force()[0].x, 0.0, 1e-18);
    EXPECT_NEAR(ps.force()[1].x, 0.0, 1e-18);
}

TEST(ContactModel, CellListFindsAllContacts) {
    // With random particles, the number of detected contacts must match a
    // brute-force search exactly.
    std::mt19937_64 rng(20240517);
    std::uniform_real_distribution<core::Real> uni(0.0, 1.0);
    const core::Real r = 2e-3;

    core::ParticleSystem ps;
    core::SimConfig cfg;
    cfg.domain = {0.05, 0.05, 0.05};
    for (int i = 0; i < 400; ++i) {
        ps.add({0.001 + 0.048 * uni(rng), 0.001 + 0.048 * uni(rng), 0.001 + 0.048 * uni(rng)},
               r, 2650.0, core::ParticleKind::Bed);
    }
    ps.setMaterial(70e9, 0.3, 0.4, 0.5);

    physics::ContactModel model(cfg);
    physics::CellList cells;
    cells.build(ps, cfg.domain.box(), physics::CellList::suggestCellSize(r));
    const auto stats = model.computeForces(ps, cells, 1e-8);

    std::size_t brute = 0;
    for (std::size_t i = 0; i < ps.size(); ++i) {
        for (std::size_t j = i + 1; j < ps.size(); ++j) {
            if ((ps.position()[i] - ps.position()[j]).norm() < 2 * r) ++brute;
        }
    }
    EXPECT_EQ(stats.touching, brute);
}

TEST(ContactModel, RollingResistanceProducesOpposingTorque) {
    auto ps = overlappingPair(1e-6);
    core::SimConfig cfg;
    cfg.particles.r_min = 1e-3;
    cfg.particles.r_max = 1e-3;
    cfg.particles.young = 70e9;
    cfg.particles.poisson = 0.3;
    cfg.particles.friction = 0.4;
    cfg.particles.restitution = 0.5;
    cfg.contact.tangential = false;
    cfg.contact.rolling_friction = true;
    cfg.contact.rolling_friction_coeff = 0.15;

    // Spin particle 0 in +y so rolling velocity omega_i - omega_j != 0
    ps.omega()[0] = {0.0, 10.0, 0.0};
    ps.omega()[1] = {0.0, 0.0, 0.0};

    physics::ContactModel model(cfg);
    physics::CellList cells;
    cells.build(ps, cfg.domain.box(), physics::CellList::suggestCellSize(1e-3));
    model.computeForces(ps, cells, 1e-7);

    // Rolling resistance torque on particle 0 must oppose its +y rotation, and be equal/opposite on 1
    EXPECT_LT(ps.torque()[0].y, 0.0);
    EXPECT_GT(ps.torque()[1].y, 0.0);
    EXPECT_NEAR(ps.torque()[0].y + ps.torque()[1].y, 0.0, 1e-18);
}

TEST(ContactModel, CohesiveProjectileBondsResistTensionAndBreakAtYield) {
    core::SimConfig cfg;
    cfg.domain = {0.02, 0.02, 0.02};
    cfg.particles.r_min = 1e-3;
    cfg.particles.r_max = 1e-3;
    cfg.particles.young = 1e7;
    cfg.particles.poisson = 0.25;
    cfg.particles.friction = 0.5;
    cfg.particles.restitution = 0.3;
    cfg.projectile.yield_kPa = 17.15;
    cfg.contact.cohesive_projectile = true;

    core::ParticleSystem ps;
    const core::Real r = 1e-3;
    // Two projectile grains initially touching at distance 2r
    ps.add({0.01, 0.01, 0.01}, r, 1500.0, core::ParticleKind::Projectile);
    ps.add({0.01, 0.01, 0.01 + 2.0 * r}, r, 1500.0, core::ParticleKind::Projectile);
    ps.setMaterial(1e7, 0.25, 0.5, 0.3);

    physics::ContactModel model(cfg);
    model.registerCohesiveBonds(ps);
    EXPECT_EQ(model.activeBonds(), 1u);

    // Small separation within tensile strength -> attractive force pulls them back together
    ps.position()[1].z = 0.01 + 2.0 * r + 1e-7;
    ps.clearForces();
    physics::CellList cells;
    cells.build(ps, cfg.domain.box(), physics::CellList::suggestCellSize(r) * 1.2);
    model.computeForces(ps, cells, 1e-6);
    EXPECT_GT(ps.force()[0].z, 0.0);   // pulled upwards toward particle 1
    EXPECT_LT(ps.force()[1].z, 0.0);   // pulled downwards toward particle 0
    EXPECT_EQ(model.activeBonds(), 1u);

    // Large separation beyond fracture threshold -> bond breaks
    ps.position()[1].z = 0.01 + 2.0 * r + 1e-4;
    ps.clearForces();
    cells.build(ps, cfg.domain.box(), physics::CellList::suggestCellSize(r) * 1.2);
    model.computeForces(ps, cells, 1e-6);
    EXPECT_EQ(model.activeBonds(), 0u);
    EXPECT_EQ(model.brokenBonds(), 1u);
}

TEST(ContactModel, PeriodicBoundaryWrapsContactsAcrossBoxEdges) {
    core::SimConfig cfg;
    cfg.domain = {0.02, 0.02, 0.02};
    cfg.boundaries.periodic_xy = true;
    const core::Real r = 1e-3;
    cfg.particles.r_min = r;
    cfg.particles.r_max = r;

    core::ParticleSystem ps;
    // One particle near x = 0.0005, one near x = 0.0195 -> wrapped distance is 0.001 < 2r = 0.002
    ps.add({0.0005, 0.01, 0.01}, r, 2650.0, core::ParticleKind::Bed);
    ps.add({0.0195, 0.01, 0.01}, r, 2650.0, core::ParticleKind::Bed);
    ps.setMaterial(1e7, 0.25, 0.5, 0.3);

    physics::ContactModel model(cfg);
    physics::CellList cells;
    cells.build(ps, cfg.domain.box(), physics::CellList::suggestCellSize(r), /*periodic_xy=*/true);
    const auto stats = model.computeForces(ps, cells, 1e-6);

    EXPECT_EQ(stats.touching, 1u);
    EXPECT_GT(ps.force()[0].x, 0.0);   // pushed into +x (away from x=0 edge)
    EXPECT_LT(ps.force()[1].x, 0.0);   // pushed into -x (away from x=Lx edge)
}

TEST(ContactModel, BedCohesionProducesCapillaryAttractionAndShearStrength) {
    auto ps = overlappingPair(1e-9);
    core::SimConfig cfg;
    cfg.particles.r_min = 1e-3;
    cfg.particles.r_max = 1e-3;
    cfg.particles.young = 1e7;
    cfg.particles.poisson = 0.25;
    cfg.particles.friction = 0.4;
    cfg.particles.restitution = 0.5;
    cfg.contact.bed_cohesion_kPa = 5.0;   // 5 kPa capillary/moisture cohesion
    ps.setMaterial(1e7, 0.25, 0.4, 0.5);

    physics::ContactModel model(cfg);
    physics::CellList cells;
    cells.build(ps, cfg.domain.box(), physics::CellList::suggestCellSize(1e-3));
    model.computeForces(ps, cells, 1e-6);

    // At tiny overlap (1e-9 m), capillary cohesion dominates Hertz repulsion and pulls grains together
    EXPECT_GT(ps.force()[0].z, 0.0);
    EXPECT_LT(ps.force()[1].z, 0.0);
}

TEST(ContactModel, MagneticDipolesAttractVerticallyAndRepelHorizontally) {
    core::SimConfig cfg;
    cfg.domain = {0.02, 0.02, 0.02};
    cfg.particles.r_min = 1e-3;
    cfg.particles.r_max = 1e-3;
    cfg.contact.magnetic_susceptibility = 1.5;
    cfg.contact.magnetic_field_T = 0.05;

    // 1) Vertically aligned pair (head-to-tail along B_z) -> attractive
    core::ParticleSystem ps_vert;
    ps_vert.add({0.01, 0.01, 0.008}, 1e-3, 2650.0, core::ParticleKind::Bed);
    ps_vert.add({0.01, 0.01, 0.0101}, 1e-3, 2650.0, core::ParticleKind::Bed);
    ps_vert.setMaterial(1e7, 0.25, 0.4, 0.5);

    physics::ContactModel model(cfg);
    physics::CellList cells;
    cells.build(ps_vert, cfg.domain.box(), physics::CellList::suggestCellSize(1e-3) * 1.2);
    model.computeForces(ps_vert, cells, 1e-6);
    EXPECT_GT(ps_vert.force()[0].z, 0.0);   // lower grain pulled up toward upper grain
    EXPECT_LT(ps_vert.force()[1].z, 0.0);   // upper grain pulled down toward lower grain

    // 2) Horizontally aligned pair (side-by-side perpendicular to B_z) -> repulsive
    core::ParticleSystem ps_horiz;
    ps_horiz.add({0.008, 0.01, 0.01}, 1e-3, 2650.0, core::ParticleKind::Bed);
    ps_horiz.add({0.0101, 0.01, 0.01}, 1e-3, 2650.0, core::ParticleKind::Bed);
    ps_horiz.setMaterial(1e7, 0.25, 0.4, 0.5);
    cells.build(ps_horiz, cfg.domain.box(), physics::CellList::suggestCellSize(1e-3) * 1.2);
    model.computeForces(ps_horiz, cells, 1e-6);
    EXPECT_LT(ps_horiz.force()[0].x, 0.0);  // left grain pushed -x
    EXPECT_GT(ps_horiz.force()[1].x, 0.0);  // right grain pushed +x
}

