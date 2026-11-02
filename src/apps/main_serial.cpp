// Reference CLI. All the physics lives in granimpact_core; this file only
// parses arguments and decides what to show. It is the pattern the openmp and
// cuda executables will follow (they share this same core).
// Unbuffered logging: in multi-hour runs the progress has to be visible live (and
// in cluster queues, in the output file, without waiting for the run to finish).
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

#include "granimpact/sim/ContactSolvers.hpp"
#include "granimpact/core/SimConfig.hpp"
#include "granimpact/sim/Simulation.hpp"

namespace {

void printUsage() {
    std::cout <<
        "GranImpact · DEM granular impact simulator\n"
        "\n"
        "Usage:\n"
        "  granimpact_serial <config.json> [options]\n"
        "\n"
        "Options:\n"
        "  --backend <name>         serial | openmp | cuda  (default: serial)\n"
        "  --out <dir>              output folder (overrides the config)\n"
        "  --name <name>            case name (file prefix)\n"
        "  --dt <seconds>           integration step\n"
        "  --t-end <seconds>        impact duration\n"
        "  --t-settle <s>           settling duration\n"
        "  --seed <int>             random generator seed\n"
        "  -v, --speed <m/s>        projectile speed (energy sweeps)\n"
        "  -a, --angle <deg>        impact angle in (0, 90] (90 = vertical)\n"
        "  -p, --phi <0..0.64>      target bed packing fraction\n"
        "  -g, --gravity <m/s2>     vertical gravity gz (e.g. -9.81 Earth, -1.62 Moon)\n"
        "  --box-x <m>              container width X [m]\n"
        "  --box-y <m>              container length Y [m]\n"
        "  --box-z <m>              container height Z [m]\n"
        "  --bed-height-frac <0..1> fraction of Z where the bed is seeded\n"
        "  --r-min <m>              minimum grain radius [m] (controls particle count N)\n"
        "  --r-max <m>              maximum grain radius [m] (controls particle count N)\n"
        "  --proj-radius <m>        projectile sphere radius R_p [m]\n"
        "  --proj-phi <0..0.64>     projectile internal packing fraction phi_p\n"
        "  --proj-density <kg/m3>   projectile bulk density rho_p [kg/m3]\n"
        "  --proj-yield <kPa>       projectile cohesive yield strength sigma_y [kPa]\n"
        "  --compaction-kpa <kPa>   uniaxial lid compaction pressure during settling [kPa]\n"
        "  --gravity-boost <factor> settling gravity multiplier\n"
        "  --bed-cohesion <kPa>     capillary / moisture bed cohesion [kPa]\n"
        "  --mag-chi <val>          magnetic susceptibility chi_m for magnetic sand\n"
        "  --mag-b <Tesla>          ambient vertical magnetic field B_0 [Tesla]\n"
        "  --no-vtk                 do not write VTK files (much faster)\n"
        "  --backends               list available backends and exit\n"
        "  -h, --help               this help\n"
        "\n"
        "Example:\n"
        "  granimpact_serial configs/default.json --out data/outputs\n";
}

int listBackends() {
    std::cout << "Contact backends\n";
    std::cout << "---------------------\n";
    for (const auto& status : granimpact::backends::backendStatus()) {
        std::cout << " " << status.name << "\t[" << (status.available ? "ready" : "not compiled")
                  << "]\t" << status.description << "\n";
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    using namespace granimpact;

    std::string config_path;
    std::string backend = "serial";
    std::string out_dir;
    std::string case_name;
    bool has_dt = false, has_t_end = false, has_settle = false, has_seed = false, has_speed = false;
    bool has_angle = false, has_phi = false, has_gravity = false;
    bool has_box_x = false, has_box_y = false, has_box_z = false, has_bed_hfrac = false;
    bool has_r_min = false, has_r_max = false;
    bool has_proj_rad = false, has_proj_phi = false, has_proj_rho = false, has_proj_yield = false;
    bool has_comp_kpa = false, has_g_boost = false;
    bool has_bed_coh = false, has_mag_chi = false, has_mag_b = false;
    double dt = 0, t_end = 0, t_settle = 0, speed = 0, angle_deg = 0, phi = 0, gz = 0;
    double box_x = 0, box_y = 0, box_z = 0, bed_hfrac = 0, r_min = 0, r_max = 0;
    double proj_rad = 0, proj_phi = 0, proj_rho = 0, proj_yield = 0;
    double comp_kpa = 0, g_boost = 0, bed_coh = 0, mag_chi = 0, mag_b = 0;
    std::uint64_t seed = 0;
    bool no_vtk = false;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        const auto next = [&](const char* what) -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "Missing the value of " << what << "\n";
                std::exit(2);
            }
            return argv[++i];
        };
        // Numeric values are checked here, during parsing: a typo must give a
        // clean message, not an uncaught std::stod exception.
        const auto number = [&](const char* what) -> double {
            const std::string value = next(what);
            try {
                return std::stod(value);
            } catch (const std::exception&) {
                // handled below, out of the catch, so std::exit is not "returning"
            }
            std::cerr << "Invalid number for " << what << ": '" << value << "'\n";
            std::exit(2);
        };
        const auto integer = [&](const char* what) -> std::uint64_t {
            const std::string value = next(what);
            try {
                return std::stoull(value);
            } catch (const std::exception&) {
                // handled below
            }
            std::cerr << "Invalid integer for " << what << ": '" << value << "'\n";
            std::exit(2);
        };
        if (arg == "-h" || arg == "--help") { printUsage(); return 0; }
        if (arg == "--backends") return listBackends();
        else if (arg == "-c" || arg == "--config") config_path = next("--config");
        else if (arg == "--backend") backend = next("--backend");
        else if (arg == "-o" || arg == "--out" || arg == "--output") out_dir = next("--out");
        else if (arg == "-n" || arg == "--name") case_name = next("--name");
        else if (arg == "--dt") { dt = number("--dt"); has_dt = true; }
        else if (arg == "--t-end") { t_end = number("--t-end"); has_t_end = true; }
        else if (arg == "--t-settle") { t_settle = number("--t-settle"); has_settle = true; }
        else if (arg == "--seed") { seed = integer("--seed"); has_seed = true; }
        else if (arg == "-v" || arg == "--speed" || arg == "--velocity") { speed = number("--speed"); has_speed = true; }
        else if (arg == "-a" || arg == "--angle") { angle_deg = number("--angle"); has_angle = true; }
        else if (arg == "-p" || arg == "--phi") { phi = number("--phi"); has_phi = true; }
        else if (arg == "-g" || arg == "--gravity") { gz = number("--gravity"); has_gravity = true; }
        else if (arg == "--box-x") { box_x = number("--box-x"); has_box_x = true; }
        else if (arg == "--box-y") { box_y = number("--box-y"); has_box_y = true; }
        else if (arg == "--box-z") { box_z = number("--box-z"); has_box_z = true; }
        else if (arg == "--bed-height-frac") { bed_hfrac = number("--bed-height-frac"); has_bed_hfrac = true; }
        else if (arg == "--r-min") { r_min = number("--r-min"); has_r_min = true; }
        else if (arg == "--r-max") { r_max = number("--r-max"); has_r_max = true; }
        else if (arg == "--proj-radius") { proj_rad = number("--proj-radius"); has_proj_rad = true; }
        else if (arg == "--proj-phi") { proj_phi = number("--proj-phi"); has_proj_phi = true; }
        else if (arg == "--proj-density") { proj_rho = number("--proj-density"); has_proj_rho = true; }
        else if (arg == "--proj-yield") { proj_yield = number("--proj-yield"); has_proj_yield = true; }
        else if (arg == "--compaction-kpa") { comp_kpa = number("--compaction-kpa"); has_comp_kpa = true; }
        else if (arg == "--gravity-boost") { g_boost = number("--gravity-boost"); has_g_boost = true; }
        else if (arg == "--bed-cohesion") { bed_coh = number("--bed-cohesion"); has_bed_coh = true; }
        else if (arg == "--mag-chi") { mag_chi = number("--mag-chi"); has_mag_chi = true; }
        else if (arg == "--mag-b") { mag_b = number("--mag-b"); has_mag_b = true; }
        else if (arg == "--no-vtk") no_vtk = true;
        else if (!arg.empty() && arg[0] == '-') {
            std::cerr << "Unknown option: " << arg << "\n";
            printUsage();
            return 2;
        } else if (config_path.empty()) {
            config_path = arg;
        }
    }

    // Unbuffered stdout: in multi-hour runs the progress must be visible live
    // (cluster queues, output files) without waiting for the run to finish.
    std::cout << std::unitbuf;

    if (config_path.empty()) {
        std::cerr << "Missing the configuration file.\n\n";
        printUsage();
        return 2;
    }

    try {
        auto config = core::SimConfig::fromFile(config_path);
        if (has_dt) config.simulation.dt = dt;
        if (has_t_end) config.simulation.t_end = t_end;
        if (has_settle) config.simulation.t_settle = t_settle;
        if (has_seed) config.simulation.seed = seed;
        if (has_speed) config.projectile.speed = speed;
        if (has_angle) config.projectile.angle_deg = angle_deg;
        if (has_phi) config.particles.packing_fraction = phi;
        if (has_gravity) config.gravity.gz = gz;
        if (has_box_x) config.domain.x = box_x;
        if (has_box_y) config.domain.y = box_y;
        if (has_box_z) config.domain.z = box_z;
        if (has_bed_hfrac) config.particles.bed_height_fraction = bed_hfrac;
        if (has_r_min) config.particles.r_min = r_min;
        if (has_r_max) config.particles.r_max = r_max;
        if (has_proj_rad) config.projectile.radius = proj_rad;
        if (has_proj_phi) config.projectile.phi = proj_phi;
        if (has_proj_rho) config.projectile.density = proj_rho;
        if (has_proj_yield) config.projectile.yield_kPa = proj_yield;
        if (has_comp_kpa) config.settle.compaction_pressure_kPa = comp_kpa;
        if (has_g_boost) config.settle.gravity_boost = g_boost;
        if (has_bed_coh) config.contact.bed_cohesion_kPa = bed_coh;
        if (has_mag_chi) config.contact.magnetic_susceptibility = mag_chi;
        if (has_mag_b) config.contact.magnetic_field_T = mag_b;
        if (!out_dir.empty()) config.output.dir = out_dir;
        if (!case_name.empty()) config.simulation.name = case_name;
        if (no_vtk) config.output.vtk = false;
        config.validate();

        sim::Simulation simulation(std::move(config), backend);
        const auto result = simulation.run(std::cout);

        std::cout << "\nFiles:\n";
        if (!result.csv_path.empty()) std::cout << "  CSV: " << result.csv_path << "\n";
        if (!result.summary_path.empty()) std::cout << "  Summary: " << result.summary_path << "\n";
        if (!result.vtk_collection.empty())
            std::cout << "  ParaView: paraview " << result.vtk_collection << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
