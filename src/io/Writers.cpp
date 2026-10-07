#include "granimpact/io/Writers.hpp"

#include <iomanip>
#include <sstream>

namespace granimpact::io {

using core::Real;

CsvWriter::CsvWriter(const std::filesystem::path& path, const core::SimConfig& config)
    : path_(path) {
    std::filesystem::create_directories(path.parent_path());
    out_.open(path_);
    if (!out_) throw std::runtime_error("Cannot write the CSV file: " + path_.string());
    writeHeader(config);
}

void CsvWriter::writeHeader(const core::SimConfig& config) {
    out_ << "# GranImpact " << config.simulation.name
         << " | dt=" << config.simulation.dt
         << " t_settle=" << config.simulation.t_settle
         << " t_end=" << config.simulation.t_end
         << " phi=" << config.particles.packing_fraction
         << " v=" << config.projectile.speed
         << " angle=" << config.projectile.angle_deg
         << " gz=" << config.gravity.gz << "\n";
    out_ << "t,particles,contacts,kinetic_energy,potential_energy,max_speed,"
            "d_exc,D,epsilon,h_rim,aspect_zd,v_in,ejecta_mass,classification\n";
}

void CsvWriter::writeRow(const StepRecord& r) {
    out_ << std::setprecision(10) << r.time << ',' << r.particles << ',' << r.contacts << ','
         << r.kinetic_energy << ',' << r.potential_energy << ',' << r.max_speed << ',';
    if (r.has_crater) {
        out_ << r.crater.d_exc << ',' << r.crater.D << ',' << r.crater.epsilon << ','
             << r.crater.h_rim << ',' << r.crater.aspect_zd << ',' << r.crater.v_in << ','
             << r.crater.ejecta_mass << ',' << r.crater.classification;
    } else {
        out_ << ",,,,,,,";
    }
    out_ << '\n';
    // Deliberate flush: a multi-hour run can be interrupted (Ctrl-C, OOM, cluster
    // queue) and the last rows have to be on disk, not in the buffer.
    out_.flush();
}

VtkWriter::VtkWriter(const std::filesystem::path& directory, std::string case_name)
    : directory_(directory), case_name_(std::move(case_name)) {
    std::filesystem::create_directories(directory_);
}

std::filesystem::path VtkWriter::write(const core::ParticleSystem& particles, std::uint64_t step,
                                       Real time) {
    std::ostringstream name;
    name << case_name_ << '_' << std::setw(6) << std::setfill('0') << step << ".vtu";
    const auto path = directory_ / name.str();

    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write the VTK file: " + path.string());

    const std::size_t n = particles.size();
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <Piece NumberOfPoints=\"" << n << "\" NumberOfCells=\"" << n << "\">\n";
    out << "      <FieldData>\n";
    out << "        <DataArray type=\"Float64\" Name=\"time\" NumberOfTuples=\"1\" format=\"ascii\">"
        << time << "</DataArray>\n";
    out << "      </FieldData>\n";
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = particles.position()[i];
        out << "          " << p.x << ' ' << p.y << ' ' << p.z << '\n';
    }
    out << "        </DataArray>\n      </Points>\n";
    out << "      <PointData>\n";
    out << "        <DataArray type=\"Float64\" Name=\"radius\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < n; ++i) out << "          " << particles.radius()[i] << '\n';
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Float64\" Name=\"speed\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < n; ++i) out << "          " << particles.velocity()[i].norm() << '\n';
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"type\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < n; ++i) {
        out << "          " << static_cast<int>(particles.kind()[i]) << '\n';
    }
    out << "        </DataArray>\n      </PointData>\n";
    out << "      <Cells>\n";
    out << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < n; ++i) out << "          " << i << '\n';
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < n; ++i) out << "          " << (i + 1) << '\n';
    out << "        </DataArray>\n";
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    // 1 = VTK_VERTEX: ParaView draws it with Glyph -> Sphere
    for (std::size_t i = 0; i < n; ++i) out << "          1\n";
    out << "        </DataArray>\n      </Cells>\n";
    out << "    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";

    steps_.emplace_back(step, time);
    ++written_;
    return path;
}

void VtkWriter::writeCollection(const core::SimConfig& config) const {
    const auto path = directory_ / (case_name_ + ".pvd");
    std::ofstream out(path);
    if (!out) return;
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <Collection>\n";
    for (const auto& [step, time] : steps_) {
        std::ostringstream name;
        name << case_name_ << '_' << std::setw(6) << std::setfill('0') << step << ".vtu";
        out << "    <DataSet timestep=\"" << time << "\" group=\"\" part=\"0\" file=\"" << name.str()
            << "\"/>\n";
    }
    out << "  </Collection>\n</VTKFile>\n";
    (void)config;
}

void CsvWriter::writeFooter(const core::SimConfig& config,
                            const analysis::CraterObservables& crater,
                            const CaseSummary& summary) {
    if (!out_) return;
    out_ << "# result:"
         << " case=" << config.simulation.name
         << " launch_speed_m_s=" << config.projectile.speed
         << " angle_deg=" << config.projectile.angle_deg
         << " gz=" << config.gravity.gz
         << " target_phi=" << config.particles.packing_fraction
         << " impact_speed_m_s=" << summary.impact_speed
         << " launch_energy_J=" << summary.impact_energy
         << " settle_steps=" << summary.settle_steps
         << " settle_converged=" << (summary.settle_converged ? "1" : "0")
         << " seeded_phi=" << summary.seeded_phi
         << " settled_phi=" << summary.settled_phi
         << " bed_phi=" << crater.phi_bed
         << " bed_shift_m=" << crater.bed_shift
         << " D_model_uehara_m=" << crater.d_uehara_model
         << " D_model_log_m=" << crater.d_log_model
         << " d_model_heckel_m=" << crater.d_heckel_model
         << " wall_time_s=" << summary.wall_seconds
         << " backend=" << summary.backend << "\n";
    out_.flush();
}

void writeCaseSummary(const std::filesystem::path& path, const core::SimConfig& config,
                      const analysis::CraterObservables& crater, const CaseSummary& summary) {
    std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path);
    if (!out) return;
    out << "GranImpact - case summary\n";
    out << "=========================\n\n";
    out << "Case              : " << config.simulation.name << '\n';
    out << "Backend           : " << summary.backend << '\n';
    out << "Seed              : " << config.simulation.seed << '\n';
    out << "dt                : " << config.simulation.dt << " s\n";
    out << "Simulated time    : " << summary.simulated_time << " s\n";
    out << "Wall time         : " << summary.wall_seconds << " s\n\n";

    out << "Bed\n";
    out << "  particles       : " << summary.bed_particles << '\n';
    out << "  target phi      : " << config.particles.packing_fraction << '\n';
    out << "  seeded phi      : " << summary.seeded_phi << "  (expanded nested lattice)\n";
    out << "  settled phi     : " << summary.settled_phi << "  (measured, occupied volume)\n";
    out << "  settling steps  : " << summary.settle_steps
        << (summary.settle_converged ? "  (converged)" : "  (NOT converged: raise t_settle)") << '\n';
    out << "  final height    : " << crater.bed_height * 1000.0 << " mm\n\n";

    out << "Projectile\n";
    out << "  grains          : " << summary.projectile_particles << '\n';
    out << "  radius / density: " << config.projectile.radius << " m / "
        << config.projectile.density << " kg/m3\n";
    out << "  launch speed    : " << config.projectile.speed << " m/s at "
        << config.projectile.angle_deg << " degrees\n";
    out << "  impact speed    : " << summary.impact_speed << " m/s (center of mass, first contact)\n";
    out << "  launch energy   : " << summary.impact_energy << " J\n\n";

    out << "Crater morphometry\n";
    out << "  D (diameter)    : " << crater.D * 1000.0 << " mm\n";
    out << "  d_exc (max excavation): " << crater.d_exc * 1000.0 << " mm\n";
    out << "  h_rim           : " << crater.h_rim * 1000.0 << " mm\n";
    out << "  epsilon (ellipse): " << crater.epsilon << '\n';
    out << "  Z/D             : " << crater.aspect_zd << '\n';
    out << "  V_in            : " << crater.v_in * 1e6 << " cm3\n";
    out << "  class           : " << crater.classification << '\n';
    out << "  cells           : " << crater.crater_cells << '\n';
    out << "  bed drift       : " << crater.bed_shift * 1000.0 << " mm\n\n";

    out << "Ejecta\n";
    out << "  ejected mass    : " << crater.ejecta_mass * 1000.0 << " g\n";
    out << "  mean speed      : " << crater.ejecta_speed_mean << " m/s\n\n";

    out << "Analytical models (reference, not fits)\n";
    out << "  D ~ E^(1/4) (Uehara) : " << crater.d_uehara_model * 1000.0 << " mm\n";
    out << "  D ~ log(E)  (thesis) : " << crater.d_log_model * 1000.0 << " mm\n";
    out << "  d = Heckel           : " << crater.d_heckel_model * 1000.0 << " mm\n";
    out << "\nNote: values marked 'measured' come from this run. The models are\n"
           "analytical reference predictions, not experimental data.\n";
}

}  // namespace granimpact::io
