#pragma once
 // Orchestrator: bed -> settling -> impact -> analysis -> outputs.
 // It is the only place where the temporal flow lives; everything else is loose
 // and separately testable components.
#include <filesystem>
#include <iosfwd>
#include <memory>
#include <string>

#include "granimpact/analysis/CraterAnalyzer.hpp"
#include "granimpact/sim/ContactSolvers.hpp"
#include "granimpact/core/SimConfig.hpp"
#include "granimpact/io/Writers.hpp"
#include "granimpact/physics/CellList.hpp"

namespace granimpact::sim {

struct SimulationResult {
    analysis::CraterObservables crater;
    io::StepRecord last_step;
    core::Real simulated_time{0};
    std::size_t bed_particles{0};
    std::size_t projectile_particles{0};
    core::Real impact_energy{0};       // projectile kinetic energy at release [J]
    core::Real impact_speed{0};        // |v| of the centre of mass at first contact [m/s]
    std::uint64_t settle_steps{0};     // settling steps actually used
    bool settle_converged{false};      // true if the energy tolerance was reached
    core::Real seeded_phi{0};          // density when the bed is seeded
    core::Real settled_phi{0};         // density after settling
    std::filesystem::path csv_path;
    std::filesystem::path summary_path;
    std::filesystem::path vtk_collection;
};

class Simulation {
public:
    Simulation(core::SimConfig config, std::string backend = "serial");

    SimulationResult run(std::ostream& log);

    // One full step: kick-drift-kick + boundaries + contacts.
    // `gravity_scale` allows settling under increased gravity (compaction).
    io::StepRecord step(io::StepRecord record, bool measure_crater, core::Real gravity_scale = 1.0);

    [[nodiscard]] core::ParticleSystem& particles() { return particles_; }
    [[nodiscard]] const core::SimConfig& config() const { return config_; }

    // Stability criterion: dt << Hertz Rayleigh time. Returns
    // the recommended maximum dt and whether the config respects it.
    [[nodiscard]] std::pair<core::Real, bool> rayleighStabilityCheck() const;

private:
    core::SimConfig config_;
    core::ParticleSystem particles_;
    std::unique_ptr<backends::ContactSolver> solver_;   // declared first: init order
    physics::CellList cells_;
    analysis::CraterAnalyzer analyzer_;
    std::unique_ptr<io::CsvWriter> csv_;
    std::unique_ptr<io::VtkWriter> vtk_;
    core::Real impact_energy_{0};
    core::Real impact_speed_{0};
    bool first_contact_seen_{false};
    core::Real z_surface_{0};
    std::uint64_t settle_steps_used_{0};
    bool settle_converged_{false};
};

}  // namespace granimpact::sim
