// End-to-end smoke test: if this passes, the full pipeline (bed ->
// settling -> projectile -> impact -> analysis -> files) really works.
// The case is small and with reduced stiffness on purpose: the thesis cases
// (thousands of particles, real sand) take hours and are not tests.
#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>

#include "granimpact/core/SimConfig.hpp"
#include "granimpact/sim/Simulation.hpp"

using namespace granimpact;

namespace {

core::SimConfig smokeConfig(const std::filesystem::path& out_dir) {
    core::SimConfig cfg;
    cfg.simulation.name = "smoke";
    cfg.simulation.dt = 2.0e-6;
    cfg.simulation.t_settle = 0.015;
    // t_end must cover the projectile's free fall down to the bed (about 3 mm of
    // gap) plus the impact dynamics; hence 0.01 s with v = 5 m/s.
    cfg.simulation.t_end = 0.010;
    cfg.simulation.output_every = 50;
    cfg.simulation.seed = 12345;

    cfg.domain = {0.04, 0.04, 0.03};
    cfg.particles.r_min = 0.8e-3;
    cfg.particles.r_max = 1.2e-3;
    cfg.particles.packing_fraction = 0.524;
    cfg.particles.young = 1.0e7;          // rigidez reducida -> dt razonable
    cfg.projectile.radius = 6.0e-3;
    cfg.projectile.start_z = 0.024;
    cfg.projectile.speed = 5.0;
    cfg.projectile.young = 1.0e7;

    cfg.output.dir = out_dir.string();
    cfg.output.csv = true;
    cfg.output.vtk = false;
    return cfg;
}

}  // namespace

TEST(Simulation, ArrancaSimulaYProduceResultados) {
    const auto out = std::filesystem::temp_directory_path() / "granimpact_smoke";
    std::filesystem::remove_all(out);

    auto cfg = smokeConfig(out);
    cfg.validate();

    sim::Simulation simulation(cfg);
    std::ostringstream log;
    const auto result = simulation.run(log);

    const std::string text = log.str();
    EXPECT_NE(text.find("Bed"), std::string::npos);
    EXPECT_NE(text.find("Settling"), std::string::npos);
    EXPECT_NE(text.find("Projectile"), std::string::npos);
    EXPECT_NE(text.find("Impact"), std::string::npos);

    EXPECT_GT(result.bed_particles, 500u);
    EXPECT_GT(result.projectile_particles, 10u);
    EXPECT_EQ(simulation.particles().size(), result.bed_particles + result.projectile_particles);
    // Settling may be cut earlier by energy convergence, so the total time lies
    // between t_end (if it converged immediately) and t_settle + t_end.
    EXPECT_GT(result.simulated_time, cfg.simulation.t_end);
    EXPECT_LE(result.simulated_time, cfg.simulation.t_settle + cfg.simulation.t_end + 1e-9);
    EXPECT_GT(result.settle_steps, 0u);

    // Measured impact energy: 1/2 m v^2 with m > 0 and v = 2 m/s.
    EXPECT_GT(result.impact_energy, 0.0);

    // Densities: `seeded_phi` measures the density INSIDE the seeded lattice and
    // `settled_phi` that of the settled bed over the occupied volume. They are two
    // different volumes (see docs/physics.md §5), so here it is only required that
    // both lie within the physical range of a granular bed.
    EXPECT_GT(result.seeded_phi, 0.30);
    EXPECT_LT(result.seeded_phi, 0.70);
    EXPECT_GT(result.settled_phi, 0.30);
    EXPECT_LT(result.settled_phi, 0.70);
    // What must hold: settling compacts the bed vertically, that is,
    // the final height ends up below the seeded region.
    EXPECT_LT(result.crater.bed_height, cfg.particles.bed_height_fraction * cfg.domain.z);

    // Crater geometry: finite and not absurd.
    EXPECT_GE(result.crater.crater_cells, 3u) << "the projectile must excavate something";
    EXPECT_GT(result.crater.d_exc, 0.0);
    EXPECT_GT(result.crater.D, 0.0);
    EXPECT_LT(result.crater.d_exc, cfg.domain.z);
    EXPECT_GE(result.crater.epsilon, 0.0);
    EXPECT_LT(result.crater.epsilon, 1.0);
    EXPECT_GT(result.crater.bed_height, 0.0);

    // CSV: header comment + header + rows.
    ASSERT_TRUE(std::filesystem::exists(result.csv_path));
    std::ifstream csv(result.csv_path);
    std::string first, header;
    std::getline(csv, first);
    std::getline(csv, header);
    EXPECT_EQ(first.front(), '#');
    EXPECT_NE(header.find("d_exc"), std::string::npos);
    std::size_t rows = 0;
    std::string line;
    std::string last;
    while (std::getline(csv, line)) { ++rows; last = line; }
    EXPECT_GT(rows, 3u);
    // The CSV ends with the measured-results block (real impact velocity,
    // settling steps, densities): the file explains itself.
    EXPECT_EQ(last.rfind("# result:", 0), 0u);
    EXPECT_NE(last.find("impact_speed_m_s="), std::string::npos);
    EXPECT_NE(last.find("settled_phi="), std::string::npos);

    std::filesystem::remove_all(out);
}

TEST(Simulation, ElCriterioDeEstabilidadDeRayleighAvisaDeUnDtGrande) {
    const auto out = std::filesystem::temp_directory_path() / "granimpact_rayleigh";
    auto cfg = smokeConfig(out);
    sim::Simulation simulation(cfg);

    const auto [dt_max, ok] = simulation.rayleighStabilityCheck();
    EXPECT_GT(dt_max, 0.0);
    EXPECT_TRUE(ok);

    cfg.simulation.dt = 1.0e-2;      // absurd on purpose
    sim::Simulation bad(cfg);
    const auto [dt_max_bad, ok_bad] = bad.rayleighStabilityCheck();
    EXPECT_GT(dt_max_bad, 0.0);
    EXPECT_FALSE(ok_bad);
}

TEST(Simulation, EscribeVtkYLaColeccionParaView) {
    const auto out = std::filesystem::temp_directory_path() / "granimpact_vtk";
    std::filesystem::remove_all(out);
    auto cfg = smokeConfig(out);
    cfg.output.vtk = true;
    cfg.output.vtk_every = 2000;   // just a few frames: writing 20,000 files adds nothing

    sim::Simulation simulation(cfg);
    std::ostringstream log;
    const auto result = simulation.run(log);

    EXPECT_TRUE(std::filesystem::exists(result.vtk_collection));
    std::ifstream pvd(result.vtk_collection);
    const std::string content((std::istreambuf_iterator<char>(pvd)),
                              std::istreambuf_iterator<char>());
    EXPECT_NE(content.find("<VTKFile"), std::string::npos);
    EXPECT_NE(content.find(".vtu"), std::string::npos);

    std::size_t vtu = 0;
    for (const auto& entry : std::filesystem::directory_iterator(out)) {
        if (entry.path().extension() == ".vtu") ++vtu;
    }
    EXPECT_GT(vtu, 1u);

    std::filesystem::remove_all(out);
}

TEST(Simulation, ElBackendAceleradoNoSeFingeSiNoEstaCompilado) {
    const auto out = std::filesystem::temp_directory_path() / "granimpact_backend";
    auto cfg = smokeConfig(out);
    const auto build_openmp = [&]() { sim::Simulation s(cfg, "openmp"); };
    const auto build_unknown = [&]() { sim::Simulation s(cfg, "quantum"); };
#ifndef GRANIMPACT_ENABLE_OPENMP
    EXPECT_THROW(build_openmp(), std::invalid_argument);
#else
    EXPECT_NO_THROW(build_openmp());
#endif
    EXPECT_THROW(build_unknown(), std::invalid_argument);
}

TEST(SimConfig, LosValoresPorDefectoDelConfigSonLosDeLaTesis) {
    const core::SimConfig cfg;
    EXPECT_NEAR(cfg.domain.x, 0.45, 1e-12);
    EXPECT_NEAR(cfg.domain.z, 0.25, 1e-12);
    EXPECT_NEAR(cfg.particles.r_min, 4e-4, 1e-15);
    EXPECT_NEAR(cfg.particles.r_max, 6e-4, 1e-15);
    EXPECT_NEAR(cfg.particles.rho_grain, 2650.0, 1e-9);
    EXPECT_NEAR(cfg.particles.packing_fraction, 0.524, 1e-12);
    EXPECT_NEAR(cfg.particles.bed_height_fraction, 0.60, 1e-12);
    EXPECT_NEAR(cfg.projectile.radius, 0.0354, 1e-12);
    EXPECT_NEAR(cfg.projectile.density, 1310.0, 1e-9);
    EXPECT_NEAR(cfg.projectile.yield_kPa, 17.15, 1e-9);
    EXPECT_NEAR(cfg.projectile.speed, 2.0, 1e-12);
    EXPECT_NEAR(cfg.gravity.gz, -9.81, 1e-12);
    EXPECT_NEAR(cfg.settle.gravity_boost, 1.0, 1e-12);
    EXPECT_EQ(cfg.simulation.seed, 20240517u);
}

TEST(SimConfig, ElJsonParcialConservaLosValoresPorDefecto) {
    const std::string text = R"({
        "simulation": { "name": "parcial", "dt": 1e-6 },
        "particles": { "packing_fraction": 0.574, "r_mean_placeholder": 0 },
        "projectile": { "angle_deg": 45.0 },
        "settle": { "gravity_boost": 2.5 }
    })";
    const auto cfg = core::SimConfig::fromString(text);

    EXPECT_EQ(cfg.simulation.name, "parcial");
    EXPECT_NEAR(cfg.simulation.dt, 1e-6, 1e-18);
    EXPECT_NEAR(cfg.particles.packing_fraction, 0.574, 1e-12);
    EXPECT_NEAR(cfg.projectile.angle_deg, 45.0, 1e-12);
    EXPECT_NEAR(cfg.settle.gravity_boost, 2.5, 1e-12);
    // Not specified -> thesis default value
    EXPECT_NEAR(cfg.domain.x, 0.45, 1e-12);
    EXPECT_NEAR(cfg.projectile.speed, 2.0, 1e-12);
    EXPECT_NEAR(cfg.gravity.gz, -9.81, 1e-12);
    EXPECT_NEAR(cfg.particles.r_mean(), 5.0e-4, 1e-15);
}

TEST(SimConfig, LaValidacionRechazaConfiguracionesImposibles) {
    core::SimConfig cfg;
    cfg.particles.packing_fraction = 0.9;         // more than the maximum packing
    EXPECT_THROW(cfg.validate(), std::invalid_argument);

    core::SimConfig ok;
    EXPECT_NO_THROW(ok.validate());

    ok.projectile.angle_deg = 0.0;
    EXPECT_THROW(ok.validate(), std::invalid_argument);

    ok = core::SimConfig{};
    ok.particles.r_max = 1e-5;                    // r_max < r_min
    EXPECT_THROW(ok.validate(), std::invalid_argument);

    ok = core::SimConfig{};
    ok.particles.bed_height_fraction = 0.99;      // leaves the box with no flight space
    EXPECT_THROW(ok.validate(), std::invalid_argument);

    ok = core::SimConfig{};
    ok.settle.gravity_boost = 0.0;
    EXPECT_THROW(ok.validate(), std::invalid_argument);
}

TEST(SimConfig, ElJsonSePuedeIdaYVuelta) {
    core::SimConfig cfg;
    cfg.simulation.name = "roundtrip";
    cfg.particles.friction = 0.37;
    const auto again = core::SimConfig::fromString(cfg.toJson());

    EXPECT_EQ(again.simulation.name, "roundtrip");
    EXPECT_NEAR(again.particles.friction, 0.37, 1e-15);
    EXPECT_NEAR(again.projectile.speed, cfg.projectile.speed, 1e-15);
    EXPECT_NEAR(again.particles.packing_fraction, cfg.particles.packing_fraction, 1e-15);
    EXPECT_EQ(again.boundaries.periodic_xy, cfg.boundaries.periodic_xy);
    EXPECT_NEAR(again.settle.gravity_boost, cfg.settle.gravity_boost, 1e-15);
}
