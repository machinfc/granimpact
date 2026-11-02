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
    core::Real d_exc{0};        // excavated depth [m]
    core::Real d_max{0};        // apparent depth = h_rim + d_exc [m]
    core::Real D{0};            // major diameter [m]
    core::Real d_minor{0};      // minor diameter [m]
    core::Real epsilon{0};      // eccentricity (0 = circular)
    core::Real theta_rad{0};    // major-axis orientation
    core::Real h_rim{0};        // rim height above the reference surface [m]
    core::Real rim_slope{0};    // mean inner crater wall slope tan(theta) [-]
    core::Real aspect_zd{0};    // d_exc / D (the thesis Z/D)
    core::Real v_in{0};         // cavity volume (paraboloid/cone) [m3]
    core::Real v_in_num{0};     // numerical cavity volume integral [m3]
    core::Real v_ex{0};         // raised rim volume integral above reference surface [m3]
    core::Real v_exc{0};        // excess volume V_in_num - V_ex (uniaxial compression volume) [m3]
    core::Real h_cp{0};         // central peak height above cavity floor (porous piston) [m]
    core::Real v_cp{0};         // central peak volume [m3]

    // Bed, subsurface saturation, and energy
    core::Real z_reference{0};            // mean surface before impact [m]
    core::Real bed_height{0};             // height of the settled bed [m]
    core::Real phi_bed{0};                // MEASURED bulk bed packing fraction [-]
    core::Real phi_subsurface{0};         // MEASURED packing fraction in compressed column below crater [-]
    core::Real impact_energy{0};          // projectile kinetic energy at release [J]
    core::Real impact_pressure_kPa{0};    // dynamic stagnation pressure 0.5 * rho_p * v^2 [kPa]
    core::Real rotational_energy{0};      // total rotational kinetic energy [J]
    core::Real dissipated_energy{0};      // cumulative dissipated frictional/viscoelastic work [J]
    core::Real momentum_transfer_vh{0};   // horizontal-to-vertical momentum transfer ratio kappa_vh [-]
    core::Real mean_virial_stress_kPa{0}; // mean Love-Weber compressive stress beneath crater [kPa]
    core::Real max_virial_stress_kPa{0};  // peak Love-Weber compressive stress beneath crater [kPa]
    std::size_t active_bonds{0};          // intact cohesive bonds in the projectile
    std::size_t broken_bonds{0};          // fractured cohesive bonds in the projectile (comminution)
    core::Real ejecta_mass{0};            // mass flagged as ejecta [kg]
    core::Real ejecta_speed_mean{0};      // mean ejecta speed [m/s]

    // Thesis model predictions (comparison, not results)
    core::Real d_log_model{0};      // Ec. 3.10: D ~ 0.089 log10(E)
    core::Real d_uehara_model{0};   // Uehara (2003): D ~ (rho_b/(rho_g mu^2))^(1/4) D_b^(3/4) H^(1/4)
    core::Real d_heckel_model{0};   // Heckel (PRE 2024): d = h_col (1 - phi0/phi(P))
    core::Real h_rim_model{0};      // Double-logarithmic rim saturation model [m]

    core::Real bed_shift{0};         // global bed subsidence (median of delta)
    std::string classification{"unknown"};   // Sand Mound / Simple / Complex
    std::size_t crater_cells{0};             // cells that form the cavity
};

struct RadialProfileBin {
    core::Real r_m{0};          // radial distance from impact center [m]
    core::Real z_bed_mm{0};     // bed-only elevation relative to z_ref [mm]
    core::Real z_total_mm{0};   // total elevation (including projectile/central peak) [mm]
};

struct SubsurfaceCell {
    core::Real r_m{0};          // radial bin center [m]
    core::Real z_m{0};          // vertical bin center [m]
    core::Real phi_local{0};    // local solid packing fraction [-]
    core::Real stress_kPa{0};   // Love-Weber virial compressive stress [kPa]
    core::Real coordination{0}; // mean coordination number Z [-]
    core::Real speed_m_s{0};    // mean grain speed [m/s]
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

    // Drift-corrected elevation map delta_z(x,y) = z_curr - z_ref in [mm]
    // (+ = rim/central peak elevation, - = crater depth), directly compatible with craterslab.
    [[nodiscard]] std::vector<core::Real> elevationMapMm(const core::ParticleSystem& particles,
                                                         bool include_projectile = true) const;

    // Azimuthally averaged radial elevation profile z(r) [mm] relative to z_reference.
    [[nodiscard]] std::vector<RadialProfileBin> radialProfile(const core::ParticleSystem& particles,
                                                              int n_bins = 48) const;

    // 2D cylindrical (r, z) subsurface field of packing fraction phi(r,z) and force-chain
    // virial stress sigma(r,z) [kPa] to visualize subsurface saturation and the porous piston.
    [[nodiscard]] std::vector<SubsurfaceCell> subsurfaceField(const core::ParticleSystem& particles,
                                                              int nr = 28, int nz = 28) const;

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
// Uehara et al. (2003): D = 0.92 (rho_b / (rho_g mu^2))^(1/4) D_b^(3/4) H^(1/4)
// where H = E / (m_p g) is the effective drop height [m].
core::Real diameterUeharaModel(core::Real energy_joules, core::Real bed_density,
                               core::Real grain_density, core::Real projectile_diameter,
                               core::Real friction_mu = 0.40, core::Real gravity_g = 9.81);
// Heckel (1961) applied to depth: d = h (1 - phi0 / phi_P)
core::Real depthHeckelModel(core::Real h, core::Real phi0, core::Real phi_pressure);
// Uniaxial Heckel compression law (Corrales-Machin et al., PRE 2024):
// phi(P) = phi_max - (phi_max - phi0) * exp(-K * max(0, P - P0))
core::Real packingHeckelModel(core::Real phi0, core::Real pressure_Pa,
                              core::Real threshold_Pa = 4460.0, core::Real k_heckel = 2.0e-5,
                              core::Real phi_max = 0.64);
core::Real depthHeckelUniaxialModel(core::Real h_column, core::Real phi0, core::Real pressure_Pa,
                                    core::Real threshold_Pa = 4460.0, core::Real k_heckel = 2.0e-5,
                                    core::Real phi_max = 0.64);
// Double-logarithmic rim height saturation law (Thesis Ch. 4):
// h_rim(E) = h_scale * ln(1 + ln(1 + E / E0))
core::Real rimHeightLogLogModel(core::Real energy_joules, core::Real h_scale = 0.0038,
                                core::Real e0_joules = 0.20);
}  // namespace models

}  // namespace granimpact::analysis
