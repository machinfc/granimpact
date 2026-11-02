#include "granimpact/core/SimConfig.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

#include <nlohmann/json.hpp>

namespace granimpact::core {
namespace {

using json = nlohmann::json;

 // Reads a nested key without failing if it is absent: configs may be partial.
template <typename T>
void get(const json& node, const char* key, T& out) {
    if (const auto it = node.find(key); it != node.end() && !it->is_null()) {
        out = it->get<T>();
    }
}

const json* section(const json& root, const char* name) {
    const auto it = root.find(name);
    return (it != root.end() && it->is_object()) ? &(*it) : nullptr;
}

}  // namespace

SimConfig SimConfig::fromString(const std::string& text) {
    json root;
    try {
        root = json::parse(text);
    } catch (const json::parse_error& e) {
        throw std::invalid_argument(std::string("invalid JSON: ") + e.what());
    }

    SimConfig cfg;
    if (const auto* s = section(root, "domain")) {
        get(*s, "x", cfg.domain.x); get(*s, "y", cfg.domain.y); get(*s, "z", cfg.domain.z);
    }
    if (const auto* s = section(root, "particles")) {
        get(*s, "r_min", cfg.particles.r_min); get(*s, "r_max", cfg.particles.r_max);
        get(*s, "rho_grain", cfg.particles.rho_grain);
        get(*s, "packing_fraction", cfg.particles.packing_fraction);
        get(*s, "young", cfg.particles.young); get(*s, "poisson", cfg.particles.poisson);
        get(*s, "restitution", cfg.particles.restitution); get(*s, "friction", cfg.particles.friction);
        get(*s, "bed_height_fraction", cfg.particles.bed_height_fraction);
    }
    if (const auto* s = section(root, "projectile")) {
        get(*s, "radius", cfg.projectile.radius); get(*s, "density", cfg.projectile.density);
        get(*s, "phi", cfg.projectile.phi); get(*s, "young", cfg.projectile.young);
        get(*s, "poisson", cfg.projectile.poisson); get(*s, "yield_kPa", cfg.projectile.yield_kPa);
        get(*s, "speed", cfg.projectile.speed); get(*s, "angle_deg", cfg.projectile.angle_deg);
        get(*s, "start_z", cfg.projectile.start_z);
    }
    if (const auto* s = section(root, "contact")) {
        get(*s, "tangential", cfg.contact.tangential);
        get(*s, "rolling_friction", cfg.contact.rolling_friction);
        get(*s, "rolling_friction_coeff", cfg.contact.rolling_friction_coeff);
        get(*s, "tangential_stiffness_ratio", cfg.contact.tangential_stiffness_ratio);
        get(*s, "cohesive_projectile", cfg.contact.cohesive_projectile);
        get(*s, "bed_cohesion_kPa", cfg.contact.bed_cohesion_kPa);
        get(*s, "magnetic_susceptibility", cfg.contact.magnetic_susceptibility);
        get(*s, "magnetic_field_T", cfg.contact.magnetic_field_T);
    }
    if (const auto* s = section(root, "settle")) {
        get(*s, "gravity_boost", cfg.settle.gravity_boost);
        get(*s, "compaction_pressure_kPa", cfg.settle.compaction_pressure_kPa);
        get(*s, "ke_tolerance", cfg.settle.ke_tolerance);
        get(*s, "max_steps", cfg.settle.max_steps);
        get(*s, "min_steps", cfg.settle.min_steps);
        get(*s, "velocity_damping", cfg.settle.velocity_damping);
    }
    if (const auto* s = section(root, "gravity")) { get(*s, "gz", cfg.gravity.gz); }
    if (const auto* s = section(root, "boundaries")) {
        get(*s, "periodic_xy", cfg.boundaries.periodic_xy);
        get(*s, "wall_restitution", cfg.boundaries.wall_restitution);
        get(*s, "wall_friction", cfg.boundaries.wall_friction);
        get(*s, "wall_damping", cfg.boundaries.wall_damping);
        get(*s, "fixed_bottom", cfg.boundaries.fixed_bottom);
    }
    if (const auto* s = section(root, "simulation")) {
        get(*s, "dt", cfg.simulation.dt); get(*s, "t_settle", cfg.simulation.t_settle);
        get(*s, "t_end", cfg.simulation.t_end); get(*s, "output_every", cfg.simulation.output_every);
        get(*s, "max_particles", cfg.simulation.max_particles); get(*s, "seed", cfg.simulation.seed);
        get(*s, "name", cfg.simulation.name);
    }
    if (const auto* s = section(root, "output")) {
        get(*s, "dir", cfg.output.dir); get(*s, "vtk", cfg.output.vtk);
        get(*s, "csv", cfg.output.csv); get(*s, "vtk_every", cfg.output.vtk_every);
    }

    cfg.validate();
    return cfg;
}

SimConfig SimConfig::fromFile(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::invalid_argument("Cannot open the config file: " + path.string());
    }
    std::stringstream buffer;
    buffer << in.rdbuf();
    return fromString(buffer.str());
}

void SimConfig::validate() const {
    const auto fail = [](const std::string& msg) { throw std::invalid_argument("config: " + msg); };

    if (domain.x <= 0 || domain.y <= 0 || domain.z <= 0) fail("domain must be positive");
    if (particles.r_min <= 0 || particles.r_max < particles.r_min) fail("r_min/r_max are inconsistent");
    if (particles.packing_fraction <= 0.05 || particles.packing_fraction > 0.64)
        fail("packing_fraction out of physical range (0.05 - 0.64)");
    if (particles.restitution < 0 || particles.restitution > 1) fail("restitution out of [0,1]");
    if (particles.friction < 0) fail("friction must be >= 0");
    if (particles.bed_height_fraction <= 0.05 || particles.bed_height_fraction > 0.95)
        fail("bed_height_fraction out of (0.05, 0.95]");
    if (contact.rolling_friction_coeff < 0) fail("contact.rolling_friction_coeff must be >= 0");
    if (contact.bed_cohesion_kPa < 0) fail("contact.bed_cohesion_kPa must be >= 0");
    if (contact.magnetic_susceptibility < 0) fail("contact.magnetic_susceptibility must be >= 0");
    if (contact.magnetic_field_T < 0) fail("contact.magnetic_field_T must be >= 0");
    if (settle.gravity_boost <= 0) fail("settle.gravity_boost must be > 0");
    if (settle.compaction_pressure_kPa < 0) fail("settle.compaction_pressure_kPa cannot be negative");
    if (settle.ke_tolerance < 0) fail("settle.ke_tolerance cannot be negative");
    if (settle.velocity_damping < 0) fail("settle.velocity_damping cannot be negative");
    if (boundaries.wall_restitution < 0 || boundaries.wall_restitution > 1)
        fail("boundaries.wall_restitution out of [0,1]");
    if (boundaries.wall_friction < 0) fail("boundaries.wall_friction must be >= 0");
    if (boundaries.wall_damping < 0 || boundaries.wall_damping > 1)
        fail("boundaries.wall_damping out of [0,1]");
    if (projectile.radius <= 0) fail("projectile.radius must be positive");
    if (projectile.phi <= 0.05 || projectile.phi > 0.64)
        fail("projectile.phi out of physical range (0.05 - 0.64)");
    if (projectile.speed < 0) fail("projectile.speed must be >= 0");
    if (projectile.angle_deg <= 0 || projectile.angle_deg > 90) fail("angle_deg must be in (0, 90]");
    if (simulation.dt <= 0) fail("dt must be positive");
    if (simulation.t_end <= 0) fail("t_end must be positive");
    // Note: the Rayleigh stability criterion is not validated here because it depends
    // on each particle's material; Simulation::rayleighStabilityCheck() does that.
}

std::string SimConfig::toJson() const {
    json root;
    root["domain"] = {{"x", domain.x}, {"y", domain.y}, {"z", domain.z}};
    root["particles"] = {{"r_min", particles.r_min}, {"r_max", particles.r_max},
                         {"rho_grain", particles.rho_grain},
                         {"packing_fraction", particles.packing_fraction},
                         {"young", particles.young}, {"poisson", particles.poisson},
                         {"restitution", particles.restitution}, {"friction", particles.friction},
                         {"bed_height_fraction", particles.bed_height_fraction}};
    root["projectile"] = {{"radius", projectile.radius}, {"density", projectile.density},
                          {"phi", projectile.phi}, {"young", projectile.young},
                          {"poisson", projectile.poisson}, {"yield_kPa", projectile.yield_kPa},
                          {"speed", projectile.speed}, {"angle_deg", projectile.angle_deg},
                          {"start_z", projectile.start_z}};
    root["contact"] = {{"tangential", contact.tangential},
                       {"rolling_friction", contact.rolling_friction},
                       {"rolling_friction_coeff", contact.rolling_friction_coeff},
                       {"tangential_stiffness_ratio", contact.tangential_stiffness_ratio},
                       {"cohesive_projectile", contact.cohesive_projectile},
                       {"bed_cohesion_kPa", contact.bed_cohesion_kPa},
                       {"magnetic_susceptibility", contact.magnetic_susceptibility},
                       {"magnetic_field_T", contact.magnetic_field_T}};
    root["settle"] = {{"gravity_boost", settle.gravity_boost},
                      {"compaction_pressure_kPa", settle.compaction_pressure_kPa},
                      {"ke_tolerance", settle.ke_tolerance},
                      {"max_steps", settle.max_steps},
                      {"min_steps", settle.min_steps},
                      {"velocity_damping", settle.velocity_damping}};
    root["gravity"] = {{"gz", gravity.gz}};
    root["boundaries"] = {{"periodic_xy", boundaries.periodic_xy},
                          {"wall_restitution", boundaries.wall_restitution},
                          {"wall_friction", boundaries.wall_friction},
                          {"wall_damping", boundaries.wall_damping},
                          {"fixed_bottom", boundaries.fixed_bottom}};
    root["simulation"] = {{"dt", simulation.dt}, {"t_settle", simulation.t_settle},
                          {"t_end", simulation.t_end}, {"output_every", simulation.output_every},
                          {"max_particles", simulation.max_particles},
                          {"seed", simulation.seed}, {"name", simulation.name}};
    root["output"] = {{"dir", output.dir}, {"vtk", output.vtk}, {"csv", output.csv},
                      {"vtk_every", output.vtk_every}};
    return root.dump(2);
}

}  // namespace granimpact::core
