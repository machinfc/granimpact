#pragma once
// Simulation parameters. Maps 1:1 to the JSON files in configs/.
// The defaults reproduce the experimental setup of the thesis
// (sand d<=1 mm, loose bed phi=0.524, 7.07 cm / 242.5 g projectile).
#include <string>
#include <filesystem>
#include <optional>

#include "granimpact/core/Types.hpp"

namespace granimpact::core {

struct DomainConfig {
    Real x{0.45}, y{0.45}, z{0.25};   // experimental box 45x45x25 cm
    [[nodiscard]] Box box() const { return {{0, 0, 0}, {x, y, z}}; }
};

struct ParticleConfig {
    Real r_min{0.0004}, r_max{0.0006};   // 0.4-0.6 mm radius (d <= 1 mm)
    Real rho_grain{2650.0};              // cuarzo 2.65 g/cm3
    Real packing_fraction{0.524};        // loose-bed TARGET (the thesis measures
                                         // 0.524; the value actually reached is measured)
    Real bed_height_fraction{0.60};      // fraction of the box height where
                                         // the bed is seeded (the rest is flight space)
    Real young{70.0e9};                  // grain Young's modulus [Pa]
    Real poisson{0.25};
    Real restitution{0.55};              // grain-grain restitution coefficient
    Real friction{0.40};                 // Coulomb friction

    // Nominal mean radius: used to scale the grid resolution and the crater
    // detection margins.
    [[nodiscard]] Real r_mean() const { return 0.5 * (r_min + r_max); }
};

struct ProjectileConfig {
    Real radius{0.0354};        // 3.54 cm -> D = 7.07 cm
    Real density{1310.0};       // 1.31 g/cm3
    Real phi{0.50};             // internal packing fraction
    Real young{5.0e6};          // effective modulus of the granular projectile [Pa]
    Real poisson{0.20};
    Real yield_kPa{17.15};      // sigma_y measured in the thesis
    Real speed{2.0};            // |v| en m/s
    Real angle_deg{90.0};       // 90 = vertical
    Real start_z{0.20};         // projectile centre height
};

struct ContactConfig {
    bool tangential{true};       // Mindlin + Coulomb
    bool rolling_friction{false}; // reserved for the extended stage
    Real tangential_stiffness_ratio{1.0};
};

struct GravityConfig {
    Real gz{-9.81};             // -1.62 for the lunar analogue
};

struct BoundaryConfig {
    bool periodic_xy{false};
    Real wall_restitution{0.30};
    bool fixed_bottom{true};
};

struct SettleConfig {
    // Compaction during settling: temporary gravity multiplier.
    // 1.0 = normal settling; >1 compacts the bed (reproduces the compacted
    // bed of the thesis without touching the impact physics).
    Real gravity_boost{1.0};
    // Kinetic energy per particle below which the bed counts as settled and
    // settling is cut (J per particle).
    Real ke_tolerance{1.0e-8};
    // Maximum number of settling steps (0 = use t_settle/dt).
    std::uint64_t max_steps{0};
    // Minimum steps before allowing the convergence cut. Without this the cut
    // fires on the first step (the seeded lattice is at rest) and the bed
    // collapses during the impact phase, falsifying the whole analysis.
    std::uint64_t min_steps{200};
    // Artificial damping during settling [1/s]: v *= exp(-lambda dt).
    // Used so the seeded bed relaxes and jams quickly instead of bouncing for
    // tenths of a second. It does NOT affect the impact phase.
    Real velocity_damping{0.0};
};

struct SimulationConfig {
    Real dt{1.0e-6};
    Real t_settle{0.02};
    Real t_end{0.05};
    std::uint64_t output_every{500};
    std::uint64_t max_particles{2'000'000};
    std::uint64_t seed{20240517};
    std::string name{"default"};
};

struct OutputConfig {
    std::string dir{"data/outputs"};
    bool vtk{true};
    bool csv{true};
    std::uint64_t vtk_every{2000};
};

struct SimConfig {
    DomainConfig domain;
    ParticleConfig particles;
    ProjectileConfig projectile;
    ContactConfig contact;
    SettleConfig settle;
    GravityConfig gravity;
    BoundaryConfig boundaries;
    SimulationConfig simulation;
    OutputConfig output;

    // Load from JSON. Missing keys keep the default value, so a minimal
    // config ({"simulation":{"name":"test"}}) is valid.
    static SimConfig fromFile(const std::filesystem::path& path);
    static SimConfig fromString(const std::string& json);
    [[nodiscard]] std::string toJson() const;

    // Validation. Throws std::invalid_argument with an actionable message.
    void validate() const;
};

}  // namespace granimpact::core
