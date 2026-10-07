#pragma once
// Outputs: time-series CSV (analysis) and VTK .vtu + .pvd (visualization).
// VTK XML ASCII on purpose: it is readable by eye, diffable and does not force
// compiling the VTK library.
#include <filesystem>
#include <fstream>
#include <string>

#include "granimpact/analysis/CraterAnalyzer.hpp"
#include "granimpact/core/ParticleSystem.hpp"
#include "granimpact/core/SimConfig.hpp"

namespace granimpact::io {

struct StepRecord {
    core::Real time{0};
    std::size_t particles{0};
    std::size_t contacts{0};
    core::Real kinetic_energy{0};
    core::Real potential_energy{0};
    core::Real max_speed{0};
    bool has_crater{false};
    analysis::CraterObservables crater{};
};

// Plain-text case summary: what was run, with which parameters, what was measured
// and how long it took. It is the file cited in the thesis next to the CSV: it
// survives deletion of the .vtu files and needs no tools to be read.
struct CaseSummary {
    std::string backend;
    std::size_t bed_particles{0};
    std::size_t projectile_particles{0};
    std::uint64_t settle_steps{0};
    bool settle_converged{false};
    core::Real impact_speed{0};
    core::Real impact_energy{0};
    core::Real seeded_phi{0};
    core::Real settled_phi{0};
    core::Real simulated_time{0};
    core::Real wall_seconds{0};
};

class CsvWriter {
public:
    CsvWriter(const std::filesystem::path& path, const core::SimConfig& config);
    void writeHeader(const core::SimConfig& config);
    void writeRow(const StepRecord& record);
    // Final block: closes the CSV with the measured run parameters (impact speed,
    // settling steps, measured density). That way the CSV explains itself.
    void writeFooter(const core::SimConfig& config, const analysis::CraterObservables& crater,
                     const CaseSummary& summary);
    [[nodiscard]] const std::filesystem::path& path() const { return path_; }

private:
    std::filesystem::path path_;
    std::ofstream out_;
};

class VtkWriter {
public:
    VtkWriter(const std::filesystem::path& directory, std::string case_name);
    // Writes data/<name>_<step>.vtu and returns the path of the file written.
    std::filesystem::path write(const core::ParticleSystem& particles, std::uint64_t step,
                                core::Real time);
    // .pvd collection: open THIS file in ParaView to see the full animation.
    void writeCollection(const core::SimConfig& config) const;
    [[nodiscard]] std::size_t written() const { return written_; }

private:
    std::filesystem::path directory_;
    std::string case_name_;
    std::vector<std::pair<std::uint64_t, core::Real>> steps_;
    std::size_t written_{0};
};

void writeCaseSummary(const std::filesystem::path& path, const core::SimConfig& config,
                      const analysis::CraterObservables& crater, const CaseSummary& summary);

}  // namespace granimpact::io
