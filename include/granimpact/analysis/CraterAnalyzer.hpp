#pragma once
// Crater morphometric observables + predictions from the thesis models.
//
// How it measures: it stores a surface map of the ALREADY SETTLED bed (reference) and
// compares it with the current map. The per-cell difference is the deformation, which is
// what defines crater, rim and ejecta. Only the bed is used (not the projectile)
// so as not to mask the excavation where the projectile comes to rest.
//
// Documented limitation: this is a grid + moments estimator (cheap and robust).
// Least-squares ellipse fitting and the fine radial profile arrive with the
// fine-grained analysis layer (see docs/architecture.md).
#include <string>
#include <vector>

#include "granimpact/core/ParticleSystem.hpp"
#include "granimpact/core/SimConfig.hpp"

namespace granimpact::analysis {

struct CraterObservables {
    // Geometry
    core::Real d_exc{0};        // profundidad excavada [m]
    core::Real d_max{0};        // profundidad aparente = h_rim + d_exc [m]
    core::Real D{0};            // major diameter [m]
    core::Real d_minor{0};      // minor diameter [m]
    core::Real epsilon{0};      // excentricidad (0 = circular)
    core::Real theta_rad{0};    // major-axis orientation
    core::Real h_rim{0};        // rim height above the reference surface
    core::Real aspect_zd{0};    // d_exc / D (the thesis Z/D)
    core::Real v_in{0};         // cavity volume (cone) [m3]

    // Bed and energy
    core::Real z_reference{0};       // mean surface before impact [m]
    core::Real bed_height{0};        // height of the settled bed [m]
    core::Real phi_bed{0};           // MEASURED packing fraction [-]
    core::Real impact_energy{0};     // projectile kinetic energy at release [J]
    core::Real ejecta_mass{0};       // mass flagged as ejecta [kg]
    core::Real ejecta_speed_mean{0}; // mean ejecta speed [m/s]

    // Thesis model predictions (comparison, not results)
    core::Real d_log_model{0};      // Ec. 3.10: D ~ 0.089 log10(E)
    core::Real d_uehara_model{0};   // Uehara (2003): D ~ (rho_b/rho_g)^(1/4) D_b^(3/4) E^(1/4)
    core::Real d_heckel_model{0};   // Heckel: d = h (1 - phi0/phi_P)

    core::Real bed_shift{0};         // global bed subsidence (median of delta)
    std::string classification{"unknown"};   // Sand Mound / Simple / Complex
    std::size_t crater_cells{0};             // cells that form the cavity
};

class CraterAnalyzer {
public:
    // `margin_factor` multiplies r_max: a cell only counts as excavated if it lies
    // more than that margin below the reference (avoids reading grid noise).
    explicit CraterAnalyzer(const core::SimConfig& config, int grid = 96,
                            core::Real margin_factor = 1.0);

    // Reference: surface of the settled bed and its mean level.
    void recordReferenceSurface(const core::ParticleSystem& particles);

    // Observables of the current state. `impact_energy` is the real impact energy
    // (measured when the projectile is released, not estimated).
    [[nodiscard]] CraterObservables analyze(const core::ParticleSystem& particles,
                                            core::Real impact_energy) const;

    // Surface map (H x W) from z maxima. bed_only=true ignores the projectile.
    [[nodiscard]] std::vector<core::Real> surfaceMap(const core::ParticleSystem& particles,
                                                     bool bed_only = true) const;

    [[nodiscard]] core::Real referenceHeight() const { return z_reference_; }
    [[nodiscard]] bool hasReference() const { return has_reference_; }
    [[nodiscard]] int gridSize() const { return grid_; }

private:
    const core::SimConfig& config_;
    int grid_{96};
    core::Real margin_factor_{1.0};
    core::Real z_reference_{0};
    bool has_reference_{false};
    std::vector<core::Real> reference_map_;   // NaN = cell with no data
};

// Thesis models, exposed so they can be tested without running a simulation.
namespace models {
// Eq. 3.10 of the thesis: D_L ~ 0.089 * log10(E) [m]
core::Real diameterLogModel(core::Real energy_joules);
// Uehara et al. (2003): D ~ 1.84 (rho_b/rho_g)^(1/4) D_b^(3/4) E^(1/4)
core::Real diameterUeharaModel(core::Real energy_joules, core::Real bed_density,
                               core::Real grain_density, core::Real projectile_diameter);
// Heckel (1961) applied to depth: d = h (1 - phi0 / phi_P)
core::Real depthHeckelModel(core::Real h, core::Real phi0, core::Real phi_pressure);
}  // namespace models

}  // namespace granimpact::analysis
