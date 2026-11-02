// The analyzer is tested with synthetic craters of known geometry: it is the
// only way to know that d_exc, D and epsilon measure what they claim to measure.
#include <gtest/gtest.h>

#include <cmath>

#include "granimpact/analysis/CraterAnalyzer.hpp"
#include "granimpact/core/SimConfig.hpp"

using namespace granimpact;

namespace {

core::SimConfig syntheticConfig() {
    core::SimConfig cfg;
    cfg.domain = {0.2, 0.2, 0.1};
    cfg.particles.r_min = 1.5e-3;
    cfg.particles.r_max = 1.5e-3;
    cfg.particles.packing_fraction = 0.5;
    // The analyzer estimates drift on the outer ring (beyond 4 R_proj): the test
    // projectile is made small so that ring exists in a 20 cm box and the synthetic
    // crater stays inside the perturbed zone.
    cfg.projectile.radius = 0.005;
    return cfg;
}

// Flat bed on a grid; if `depth` > 0 an elliptical depression is excavated.
core::ParticleSystem bedWithDepression(const core::SimConfig& cfg, core::Real z_ref,
                                       core::Real depth, core::Real rx, core::Real ry) {
    core::ParticleSystem ps;
    const core::Box box = cfg.domain.box();
    const int nx = 80, ny = 80;
    for (int iy = 0; iy < ny; ++iy) {
        for (int ix = 0; ix < nx; ++ix) {
            const core::Real x = box.min.x + (ix + 0.5) * box.size().x / nx;
            const core::Real y = box.min.y + (iy + 0.5) * box.size().y / ny;
            core::Real z = z_ref;
            if (depth > 0) {
                const core::Real dx = (x - box.size().x / 2) / rx;
                const core::Real dy = (y - box.size().y / 2) / ry;
                const core::Real s = dx * dx + dy * dy;
                if (s < 1.0) z -= depth * std::sqrt(1.0 - s);   // hemispherical cap
            }
            ps.add({x, y, z}, cfg.particles.r_min, 2650.0, core::ParticleKind::Bed);
        }
    }
    ps.setMaterial(1.0e7, 0.3, 0.4, 0.5);
    return ps;
}

}  // namespace

TEST(CraterModels, DiameterFollowsTheLogModelFromTheThesis) {
    // E = 1 kJ -> D = 0.089 * 3 = 0.267 m
    EXPECT_NEAR(analysis::models::diameterLogModel(1000.0), 0.089 * 3.0, 1e-12);
    EXPECT_GT(analysis::models::diameterLogModel(5000.0),
              analysis::models::diameterLogModel(1000.0));
    // Below 1 J the logarithmic law is not calibrated: it returns 0 instead
    // of a meaningless value (the thesis fitted it with energies of tens of J).
    EXPECT_EQ(analysis::models::diameterLogModel(0.5), 0.0);
}

TEST(CraterModels, UeharaDiameterScalesWithEnergyToTheOneFourth) {
    const core::Real d1 = analysis::models::diameterUeharaModel(10.0, 1500.0, 2650.0, 0.07);
    const core::Real d16 = analysis::models::diameterUeharaModel(160.0, 1500.0, 2650.0, 0.07);
    EXPECT_GT(d1, 0.0);
    EXPECT_NEAR(d16 / d1, 2.0, 1e-12);   // 16x energy -> 2x diameter
}

TEST(CraterModels, DepthFollowsHeckel) {
    EXPECT_NEAR(analysis::models::depthHeckelModel(0.05, 0.5, 0.5), 0.0, 1e-15);
    EXPECT_GT(analysis::models::depthHeckelModel(0.05, 0.524, 0.60), 0.0);
    EXPECT_NEAR(analysis::models::depthHeckelModel(0.05, 0.30, 0.60), 0.025, 1e-15);
}

TEST(CraterAnalyzer, RecoversGeometryOfAKnownDepression) {
    const auto cfg = syntheticConfig();
    const core::Real z_ref = 0.05;
    const core::Real depth = 0.02;
    const core::Real rx = 0.015;

    const auto flat = bedWithDepression(cfg, z_ref, 0.0, rx, rx);
    auto deformed = bedWithDepression(cfg, z_ref, depth, rx, rx);

    analysis::CraterAnalyzer analyzer(cfg, 128);
    analyzer.recordReferenceSurface(flat);
    const auto obs = analyzer.analyze(deformed, 5.0);

    EXPECT_NEAR(analyzer.referenceHeight(), z_ref, 1e-9);
    EXPECT_GT(obs.crater_cells, 50u);
    // The excavated depth reproduces the cap's (within the grid margin).
    EXPECT_NEAR(obs.d_exc, depth, 0.10 * depth);
    // The depression diameter is several times the nominal cap radius.
    EXPECT_GT(obs.D, rx);
    EXPECT_LT(obs.D, 4.0 * rx);
    EXPECT_LT(obs.epsilon, 0.2);                 // circular -> epsilon ~ 0
    EXPECT_GT(obs.v_in, 0.0);
    EXPECT_NEAR(obs.aspect_zd, obs.d_exc / obs.D, 1e-12);
    EXPECT_NEAR(obs.impact_energy, 5.0, 1e-12);  // the energy passed in, as is
    EXPECT_GT(obs.bed_height, 0.0);
    EXPECT_GT(obs.phi_bed, 0.0);
}

TEST(CraterAnalyzer, DetectsAnElongatedDepression) {
    const auto cfg = syntheticConfig();
    auto flat = bedWithDepression(cfg, 0.05, 0.0, 0.025, 0.025);
    auto elongated = bedWithDepression(cfg, 0.05, 0.01, 0.025, 0.010);

    analysis::CraterAnalyzer analyzer(cfg, 128);
    analyzer.recordReferenceSurface(flat);
    const auto obs = analyzer.analyze(elongated, 2.0);

    EXPECT_GT(obs.crater_cells, 20u);
    EXPECT_GT(obs.epsilon, 0.5);                 // clearly elliptical
    EXPECT_GT(obs.D, obs.d_minor);
}

TEST(CraterAnalyzer, AFlatSurfaceProducesNoCrater) {
    const auto cfg = syntheticConfig();
    auto flat = bedWithDepression(cfg, 0.05, 0.0, 0.02, 0.02);

    analysis::CraterAnalyzer analyzer(cfg, 64);
    analyzer.recordReferenceSurface(flat);
    const auto obs = analyzer.analyze(flat, 5.0);

    EXPECT_EQ(obs.crater_cells, 0u);
    EXPECT_EQ(obs.d_exc, 0.0);
    EXPECT_EQ(obs.classification, "Sand Mound / no crater");
    EXPECT_GT(obs.bed_height, 0.0);              // the bed is reported all the same
}

TEST(CraterAnalyzer, MarginFiltersGrainScaleNoise) {
    const auto cfg = syntheticConfig();
    auto flat = bedWithDepression(cfg, 0.05, 0.0, 0.02, 0.02);
    // Subsidence by half a particle: it must not be read as a crater.
    auto tiny = bedWithDepression(cfg, 0.05, 0.4 * cfg.particles.r_min, 0.02, 0.02);

    analysis::CraterAnalyzer analyzer(cfg, 64);
    analyzer.recordReferenceSurface(flat);
    const auto obs = analyzer.analyze(tiny, 1.0);
    EXPECT_EQ(obs.crater_cells, 0u);
}

TEST(CraterAnalyzer, ClassifiesByAspectRatioZ) {
    const auto cfg = syntheticConfig();
    // Deep, narrow crater -> high Z/D -> 'Simple (deep)'.
    auto flat = bedWithDepression(cfg, 0.08, 0.0, 0.008, 0.008);
    auto deep = bedWithDepression(cfg, 0.08, 0.03, 0.008, 0.008);

    analysis::CraterAnalyzer analyzer(cfg, 128);
    analyzer.recordReferenceSurface(flat);
    const auto obs = analyzer.analyze(deep, 20.0);

    ASSERT_GT(obs.crater_cells, 0u);
    const char* expected = obs.aspect_zd > 0.20 ? "Simple (deep)"
                           : (obs.aspect_zd > 0.05 ? "Simple" : "Sand Mound");
    EXPECT_EQ(obs.classification, expected);
}

TEST(CraterModels, UniaxialHeckelModelMatchesPhysRevEPaper) {
    // At P = 0, phi(0) == phi_0 and depth == 0
    EXPECT_NEAR(analysis::models::packingHeckelModel(0.524, 0.0), 0.524, 1e-12);
    EXPECT_NEAR(analysis::models::depthHeckelUniaxialModel(0.11, 0.524, 0.0), 0.0, 1e-12);

    // Loose bed (phi_L = 0.524) vs Compacted bed (phi_C = 0.574) under P_0 = 50 kPa (50000 Pa):
    // Loose bed must yield a deeper crater than compacted bed.
    const core::Real z_loose = analysis::models::depthHeckelUniaxialModel(0.11, 0.524, 50000.0);
    const core::Real z_comp  = analysis::models::depthHeckelUniaxialModel(0.11, 0.574, 50000.0);
    EXPECT_GT(z_loose, z_comp);
    EXPECT_GT(z_comp, 0.0);

    // High pressure asymptotically approaches the jamming limit: Z_max = H0 * (1 - phi_0 / phi_m)
    const core::Real z_sat = analysis::models::depthHeckelUniaxialModel(0.11, 0.524, 1e8);
    EXPECT_NEAR(z_sat, 0.11 * (1.0 - 0.524 / 0.64), 1e-6);
}

TEST(CraterModels, RimHeightFollowsDoubleLogLawFromThesis) {
    // Zero energy -> 0
    EXPECT_EQ(analysis::models::rimHeightLogLogModel(0.0), 0.0);
    // Above E0, h_rim grows double-logarithmically and saturates softly
    const core::Real h10 = analysis::models::rimHeightLogLogModel(10.0);
    const core::Real h50 = analysis::models::rimHeightLogLogModel(50.0);
    EXPECT_GT(h10, 0.0);
    EXPECT_GT(h50, h10);
}

TEST(CraterAnalyzer, ExtractsRadialProfileAndSubsurfaceField) {
    const auto cfg = syntheticConfig();
    const auto flat = bedWithDepression(cfg, 0.05, 0.0, 0.02, 0.02);
    const auto deformed = bedWithDepression(cfg, 0.05, 0.015, 0.02, 0.02);

    analysis::CraterAnalyzer analyzer(cfg, 64);
    analyzer.recordReferenceSurface(flat);

    const auto profile = analyzer.radialProfile(deformed, 16);
    ASSERT_EQ(profile.size(), 16u);
    // Center of crater (bin 0) has negative z_bed_mm relative to z_ref
    EXPECT_LT(profile[0].z_bed_mm, -5.0);

    const auto sub = analyzer.subsurfaceField(deformed, 8, 6);
    ASSERT_EQ(sub.size(), 48u);
    bool has_particles = false;
    for (const auto& cell : sub) {
        if (cell.phi_local > 0.0) has_particles = true;
    }
    EXPECT_TRUE(has_particles);
}

TEST(CraterAnalyzer, DetectsCentralPeakFromProjectileRemnants) {
    const auto cfg = syntheticConfig();
    auto ps = bedWithDepression(cfg, 0.05, 0.015, 0.025, 0.025);
    // Add a cluster of Projectile remnants at the center forming a porous piston / central peak
    for (int ix = -2; ix <= 2; ++ix) {
        for (int iy = -2; iy <= 2; ++iy) {
            ps.add({0.10 + ix * 0.002, 0.10 + iy * 0.002, 0.045},
                   cfg.particles.r_min, 1500.0, core::ParticleKind::Projectile);
        }
    }

    analysis::CraterAnalyzer analyzer(cfg, 64);
    const auto flat = bedWithDepression(cfg, 0.05, 0.0, 0.025, 0.025);
    analyzer.recordReferenceSurface(flat);
    const auto obs = analyzer.analyze(ps, 10.0);

    EXPECT_GT(obs.h_cp, 0.005);
    EXPECT_GT(obs.v_cp, 0.0);
}

