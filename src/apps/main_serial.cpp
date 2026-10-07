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
        "  --backend <name>     serial | openmp | cuda  (default: serial)\n"
        "  --out <dir>          output folder (overrides the config)\n"
        "  --name <name>        case name (file prefix)\n"
        "  --dt <seconds>       integration step\n"
        "  --t-end <seconds>    impact duration\n"
        "  --t-settle <s>       settling duration\n"
        "  --seed <int>         random generator seed\n"
        "  --speed <m/s>        projectile speed (energy sweeps)\n"
        "  --no-vtk             do not write VTK files (much faster)\n"
        "  --backends           list available backends and exit\n"
        "  -h, --help           this help\n"
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
    double dt = 0, t_end = 0, t_settle = 0, speed = 0;
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
        else if (arg == "--backend") backend = next("--backend");
        else if (arg == "--out") out_dir = next("--out");
        else if (arg == "--name") case_name = next("--name");
        else if (arg == "--dt") { dt = number("--dt"); has_dt = true; }
        else if (arg == "--t-end") { t_end = number("--t-end"); has_t_end = true; }
        else if (arg == "--t-settle") { t_settle = number("--t-settle"); has_settle = true; }
        else if (arg == "--seed") { seed = integer("--seed"); has_seed = true; }
        else if (arg == "--speed") { speed = number("--speed"); has_speed = true; }
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
