#include "granimpact/analysis/CraterAnalyzer.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

namespace granimpact::analysis {

using core::kPi;
using core::ParticleKind;
using core::Real;
using core::Vec3;

namespace models {

Real diameterLogModel(Real energy_joules) {
    if (energy_joules <= 1.0) return 0;   // the log law was calibrated for E > 1 J
    return 0.089 * std::log10(energy_joules);
}

Real diameterUeharaModel(Real energy_joules, Real bed_density, Real grain_density,
                         Real projectile_diameter, Real friction_mu, Real gravity_g) {
    if (energy_joules <= 0 || bed_density <= 0 || grain_density <= 0
        || projectile_diameter <= 0 || friction_mu <= 0 || gravity_g <= 0) {
        return 0;
    }
    // Effective drop height H = E / (m_p * g), with m_p = (pi/6) * rho_b * D_b^3.
    // Preserves exact D ~ E^(1/4) scaling while restoring physical dimensions [m].
    const Real m_proj = (kPi / Real(6.0)) * bed_density
                        * projectile_diameter * projectile_diameter * projectile_diameter;
    const Real h_eff = energy_joules / (m_proj * gravity_g);
    const Real ratio = bed_density / (grain_density * friction_mu * friction_mu);
    return Real(0.92) * std::pow(ratio, Real(0.25))
           * std::pow(projectile_diameter, Real(0.75))
           * std::pow(h_eff, Real(0.25));
}

Real depthHeckelModel(Real h, Real phi0, Real phi_pressure) {
    if (phi_pressure <= 0) return 0;
    const Real ratio = std::clamp(phi0 / phi_pressure, Real(0), Real(1));
    return h * (1.0 - ratio);
}

Real packingHeckelModel(Real phi0, Real pressure_Pa, Real threshold_Pa,
                        Real k_heckel, Real phi_max) {
    if (phi0 >= phi_max) return phi0;
    const Real eff_p = std::max(Real(0), pressure_Pa - threshold_Pa);
    return phi_max - (phi_max - phi0) * std::exp(-k_heckel * eff_p);
}

Real depthHeckelUniaxialModel(Real h_column, Real phi0, Real pressure_Pa,
                              Real threshold_Pa, Real k_heckel, Real phi_max) {
    const Real phi_p = packingHeckelModel(phi0, pressure_Pa, threshold_Pa, k_heckel, phi_max);
    return depthHeckelModel(h_column, phi0, phi_p);
}

Real rimHeightLogLogModel(Real energy_joules, Real h_scale, Real e0_joules) {
    if (energy_joules <= 0 || e0_joules <= 0) return 0;
    return h_scale * std::log(Real(1.0) + std::log(Real(1.0) + energy_joules / e0_joules));
}

}  // namespace models

namespace {

// Spatial smoothing of a surface map. Averaging over a neighbourhood of a few
// grains removes discrete grain-scale steps and recovers the continuous topography.
std::vector<core::Real> smoothMap(const std::vector<core::Real>& map, int grid, int radius) {
    std::vector<core::Real> out(map.size(), std::numeric_limits<core::Real>::quiet_NaN());
    for (int cy = 0; cy < grid; ++cy) {
        for (int cx = 0; cx < grid; ++cx) {
            core::Real sum = 0;
            int n = 0;
            for (int dy = -radius; dy <= radius; ++dy) {
                for (int dx = -radius; dx <= radius; ++dx) {
                    const int x = cx + dx, y = cy + dy;
                    if (x < 0 || y < 0 || x >= grid || y >= grid) continue;
                    const core::Real z = map[static_cast<std::size_t>(y)
                                            * static_cast<std::size_t>(grid)
                                            + static_cast<std::size_t>(x)];
                    if (std::isnan(z)) continue;
                    sum += z;
                    ++n;
                }
            }
            if (n > 0) {
                out[static_cast<std::size_t>(cy) * static_cast<std::size_t>(grid)
                    + static_cast<std::size_t>(cx)] = sum / n;
            }
        }
    }
    return out;
}

}  // namespace

CraterAnalyzer::CraterAnalyzer(const core::SimConfig& config, int grid, Real margin_factor)
    : config_(config), grid_(grid), margin_factor_(margin_factor) {}

std::vector<Real> CraterAnalyzer::surfaceMap(const core::ParticleSystem& particles,
                                             bool bed_only) const {
    const core::Box box = config_.domain.box();
    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    std::vector<Real> map(static_cast<std::size_t>(grid_) * static_cast<std::size_t>(grid_), nan);

    const Real cell_x = box.size().x / static_cast<Real>(grid_);
    const Real cell_y = box.size().y / static_cast<Real>(grid_);
    const auto& pos = particles.position();
    const auto& rad = particles.radius();
    const auto& kinds = particles.kind();

    // Pass 1: exact center cell projection.
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (bed_only && kinds[i] != ParticleKind::Bed) continue;
        if (!bed_only && kinds[i] == ParticleKind::Ejecta) continue;
        if (has_reference_ && pos[i].z > z_reference_ + Real(0.8) * config_.projectile.radius) {
            continue;   // ignore airborne grains high above the bed
        }
        const Vec3 p = pos[i];
        const auto cx = static_cast<int>((p.x - box.min.x) / box.size().x * grid_);
        const auto cy = static_cast<int>((p.y - box.min.y) / box.size().y * grid_);
        if (cx < 0 || cy < 0 || cx >= grid_ || cy >= grid_) continue;
        Real& cell = map[static_cast<std::size_t>(cy) * static_cast<std::size_t>(grid_)
                         + static_cast<std::size_t>(cx)];
        cell = std::isnan(cell) ? p.z : std::max(cell, p.z);
    }

    // Pass 2: physical sphere footprint occlusion (r_reach = 1.15 * r_i) for empty cells.
    std::vector<Real> env(map.size(), nan);
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (bed_only && kinds[i] != ParticleKind::Bed) continue;
        if (!bed_only && kinds[i] == ParticleKind::Ejecta) continue;
        if (has_reference_ && pos[i].z > z_reference_ + Real(0.8) * config_.projectile.radius) {
            continue;
        }
        const Vec3 p = pos[i];
        const Real r_reach = Real(1.15) * rad[i];
        if (r_reach < Real(0.45) * std::min(cell_x, cell_y)) continue;

        const int cx0 = static_cast<int>((p.x - box.min.x) / cell_x);
        const int cy0 = static_cast<int>((p.y - box.min.y) / cell_y);
        const int rx_cells = std::max(1, static_cast<int>(std::ceil(r_reach / cell_x)));
        const int ry_cells = std::max(1, static_cast<int>(std::ceil(r_reach / cell_y)));

        for (int dy = -ry_cells; dy <= ry_cells; ++dy) {
            const int cy = cy0 + dy;
            if (cy < 0 || cy >= grid_) continue;
            const Real yc = box.min.y + (static_cast<Real>(cy) + Real(0.5)) * cell_y;
            for (int dx = -rx_cells; dx <= rx_cells; ++dx) {
                const int cx = cx0 + dx;
                if (cx < 0 || cx >= grid_) continue;
                const Real xc = box.min.x + (static_cast<Real>(cx) + Real(0.5)) * cell_x;
                const Real d2 = (xc - p.x) * (xc - p.x) + (yc - p.y) * (yc - p.y);
                if (d2 > r_reach * r_reach) continue;

                Real& ecell = env[static_cast<std::size_t>(cy) * static_cast<std::size_t>(grid_)
                                  + static_cast<std::size_t>(cx)];
                ecell = std::isnan(ecell) ? p.z : std::max(ecell, p.z);
            }
        }
    }
    for (std::size_t c = 0; c < map.size(); ++c) {
        if (std::isnan(map[c])) {
            map[c] = env[c];
        }
    }

    // Pass 3: 3x3 pinhole pore filter — closes isolated inter-grain pores where a vertical
    // ray fell between touching top-layer spheres into deep bed layers, without bridging craters.
    std::vector<Real> filtered = map;
    const Real pore_drop = Real(2.0) * config_.particles.r_max;
    for (int cy = 0; cy < grid_; ++cy) {
        for (int cx = 0; cx < grid_; ++cx) {
            const std::size_t c = static_cast<std::size_t>(cy) * static_cast<std::size_t>(grid_)
                                  + static_cast<std::size_t>(cx);
            Real nb_sum = 0;
            int nb_count = 0;
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dx = -1; dx <= 1; ++dx) {
                    if (dx == 0 && dy == 0) continue;
                    const int nx = cx + dx;
                    const int ny = cy + dy;
                    if (nx < 0 || ny < 0 || nx >= grid_ || ny >= grid_) continue;
                    const Real z_nb = map[static_cast<std::size_t>(ny) * static_cast<std::size_t>(grid_)
                                          + static_cast<std::size_t>(nx)];
                    if (!std::isnan(z_nb)) {
                        nb_sum += z_nb;
                        ++nb_count;
                    }
                }
            }
            if (nb_count >= 5) {
                const Real nb_mean = nb_sum / static_cast<Real>(nb_count);
                if (std::isnan(map[c]) || (nb_mean - map[c] > pore_drop)) {
                    filtered[c] = nb_mean;
                }
            }
        }
    }
    return filtered;
}

void CraterAnalyzer::recordReferenceSurface(const core::ParticleSystem& particles) {
    has_reference_ = false;
    const int smooth_radius = std::max(1, static_cast<int>(std::lround(
        Real(2.0) * config_.particles.r_max / (config_.domain.x / grid_))));
    reference_map_ = smoothMap(surfaceMap(particles, /*bed_only=*/true), grid_, smooth_radius);

    Real sum = 0;
    std::size_t n = 0;
    for (const Real z : reference_map_) {
        if (std::isnan(z)) continue;
        sum += z;
        ++n;
    }
    z_reference_ = n > 0 ? sum / static_cast<Real>(n) : 0;
    has_reference_ = n > 0;
}

std::vector<Real> CraterAnalyzer::elevationMapMm(const core::ParticleSystem& particles,
                                                 bool include_projectile) const {
    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    std::vector<Real> elev(reference_map_.size(), nan);
    if (!has_reference_) return elev;

    const core::Box box = config_.domain.box();
    const Real cell_x = box.size().x / grid_;
    const Real cell_y = box.size().y / grid_;
    const int smooth_radius = std::max(
        1, static_cast<int>(std::lround(Real(2.0) * config_.particles.r_max
                                        / std::min(cell_x, cell_y))));
    const auto current = smoothMap(
        surfaceMap(particles, /*bed_only=*/!include_projectile), grid_, smooth_radius);

    const Vec3 impact{box.min.x + Real(0.5) * box.size().x,
                      box.min.y + Real(0.5) * box.size().y, 0};
    const Real exclude = Real(4.0) * config_.projectile.radius;
    const bool large_box = std::min(box.size().x, box.size().y) >= Real(10.0) * config_.projectile.radius;
    const Real wall_pad = large_box ? Real(2.5) * config_.particles.r_max : Real(0);
    const auto insideUsableBox = [&](int cx, int cy) {
        if (wall_pad <= 0) return true;
        const Real x = box.min.x + (cx + Real(0.5)) * cell_x;
        const Real y = box.min.y + (cy + Real(0.5)) * cell_y;
        return (x >= box.min.x + wall_pad && x <= box.max.x - wall_pad &&
                y >= box.min.y + wall_pad && y <= box.max.y - wall_pad);
    };

    double outer_sum = 0;
    std::size_t outer_n = 0;
    for (int cy = 0; cy < grid_; ++cy) {
        for (int cx = 0; cx < grid_; ++cx) {
            if (!insideUsableBox(cx, cy)) continue;
            const std::size_t c = static_cast<std::size_t>(cy) * static_cast<std::size_t>(grid_)
                                  + static_cast<std::size_t>(cx);
            if (std::isnan(reference_map_[c]) || std::isnan(current[c])) continue;
            const Real x = box.min.x + (cx + Real(0.5)) * cell_x;
            const Real y = box.min.y + (cy + Real(0.5)) * cell_y;
            if (std::hypot(x - impact.x, y - impact.y) < exclude) continue;
            outer_sum += static_cast<double>(reference_map_[c] - current[c]);
            ++outer_n;
        }
    }
    const Real drift = outer_n >= 10 ? static_cast<Real>(outer_sum / static_cast<double>(outer_n))
                                     : Real(0);

    for (int cy = 0; cy < grid_; ++cy) {
        for (int cx = 0; cx < grid_; ++cx) {
            const std::size_t c = static_cast<std::size_t>(cy) * static_cast<std::size_t>(grid_)
                                  + static_cast<std::size_t>(cx);
            if (!insideUsableBox(cx, cy)) {
                elev[c] = Real(0);
                continue;
            }
            if (std::isnan(reference_map_[c]) || std::isnan(current[c])) continue;
            // Elevation relative to drift-corrected reference surface in [mm]:
            // negative inside cavity, positive on rim and central peak.
            elev[c] = ((current[c] - reference_map_[c]) + drift) * Real(1000.0);
        }
    }
    return elev;
}

std::vector<RadialProfileBin> CraterAnalyzer::radialProfile(const core::ParticleSystem& particles,
                                                            int n_bins) const {
    std::vector<RadialProfileBin> bins(static_cast<std::size_t>( std::max(1, n_bins) ));
    if (!has_reference_ || n_bins <= 0) return bins;

    const auto elev_bed = elevationMapMm(particles, /*include_projectile=*/false);
    const auto elev_tot = elevationMapMm(particles, /*include_projectile=*/true);

    const core::Box box = config_.domain.box();
    const Real cell_x = box.size().x / grid_;
    const Real cell_y = box.size().y / grid_;
    const Real cx_center = box.min.x + Real(0.5) * box.size().x;
    const Real cy_center = box.min.y + Real(0.5) * box.size().y;
    const Real r_max_domain = Real(0.5) * std::min(box.size().x, box.size().y);
    const Real dr = r_max_domain / static_cast<Real>(n_bins);

    std::vector<double> sum_bed(static_cast<std::size_t>(n_bins), 0.0);
    std::vector<double> sum_tot(static_cast<std::size_t>(n_bins), 0.0);
    std::vector<std::size_t> count_bed(static_cast<std::size_t>(n_bins), 0);
    std::vector<std::size_t> count_tot(static_cast<std::size_t>(n_bins), 0);

    for (int cy = 0; cy < grid_; ++cy) {
        for (int cx = 0; cx < grid_; ++cx) {
            const std::size_t c = static_cast<std::size_t>(cy) * static_cast<std::size_t>(grid_)
                                  + static_cast<std::size_t>(cx);
            const Real x = box.min.x + (cx + Real(0.5)) * cell_x;
            const Real y = box.min.y + (cy + Real(0.5)) * cell_y;
            const Real r = std::hypot(x - cx_center, y - cy_center);
            const int b = static_cast<int>(r / dr);
            if (b < 0 || b >= n_bins) continue;
            const auto ub = static_cast<std::size_t>(b);
            if (!std::isnan(elev_bed[c])) {
                sum_bed[ub] += static_cast<double>(elev_bed[c]);
                ++count_bed[ub];
            }
            if (!std::isnan(elev_tot[c])) {
                sum_tot[ub] += static_cast<double>(elev_tot[c]);
                ++count_tot[ub];
            }
        }
    }

    for (int b = 0; b < n_bins; ++b) {
        const auto ub = static_cast<std::size_t>(b);
        bins[ub].r_m = (static_cast<Real>(b) + Real(0.5)) * dr;
        bins[ub].z_bed_mm = count_bed[ub] > 0
                                ? static_cast<Real>(sum_bed[ub] / static_cast<double>(count_bed[ub]))
                                : Real(0);
        bins[ub].z_total_mm = count_tot[ub] > 0
                                  ? static_cast<Real>(sum_tot[ub] / static_cast<double>(count_tot[ub]))
                                  : bins[ub].z_bed_mm;
    }
    return bins;
}

std::vector<SubsurfaceCell> CraterAnalyzer::subsurfaceField(const core::ParticleSystem& particles,
                                                            int nr, int nz) const {
    nr = std::max(1, nr);
    nz = std::max(1, nz);
    std::vector<SubsurfaceCell> grid(static_cast<std::size_t>(nr) * static_cast<std::size_t>(nz));

    const core::Box box = config_.domain.box();
    const Real cx_center = box.min.x + Real(0.5) * box.size().x;
    const Real cy_center = box.min.y + Real(0.5) * box.size().y;
    const Real r_max_bin = Real(0.5) * std::min(box.size().x, box.size().y);
    const Real z_top = has_reference_
                           ? (z_reference_ + config_.particles.r_mean())
                           : (box.min.z + config_.particles.bed_height_fraction * box.size().z);
    const Real z_span = std::max(z_top - box.min.z, Real(1e-6));
    const Real dr = r_max_bin / static_cast<Real>(nr);
    const Real dz = z_span / static_cast<Real>(nz);

    std::vector<double> solid_vol(grid.size(), 0.0);
    std::vector<double> virial_sum(grid.size(), 0.0);
    std::vector<double> coord_sum(grid.size(), 0.0);
    std::vector<double> speed_sum(grid.size(), 0.0);
    std::vector<std::size_t> counts(grid.size(), 0);

    const auto& pos = particles.position();
    const auto& vel = particles.velocity();
    const auto& rad = particles.radius();
    const auto& kinds = particles.kind();
    const auto& fn_sum = particles.normalForceSum();
    const auto& coord = particles.coordination();

    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (kinds[i] == ParticleKind::Ejecta) continue;
        const Real r = std::hypot(pos[i].x - cx_center, pos[i].y - cy_center);
        const int ir = static_cast<int>(r / dr);
        if (ir < 0 || ir >= nr) continue;

        const Real ri = rad[i];
        const Real vi = Real(4.0 / 3.0) * kPi * ri * ri * ri;
        const Real z_lo = pos[i].z - ri;
        const Real z_hi = pos[i].z + ri;
        const int iz_min = std::max(0, static_cast<int>((z_lo - box.min.z) / dz));
        const int iz_max = std::min(nz - 1, static_cast<int>((z_hi - box.min.z) / dz));

        // Distribute sphere solid volume and virial stress across overlapped vertical bins
        // using Archimedes' exact spherical slice integral V_slice = pi * (R^2 (b-a) - (b^3-a^3)/3).
        for (int iz = iz_min; iz <= iz_max; ++iz) {
            const Real bin_z0 = box.min.z + static_cast<Real>(iz) * dz;
            const Real bin_z1 = bin_z0 + dz;
            const Real a_sl = std::clamp(bin_z0 - pos[i].z, -ri, ri);
            const Real b_sl = std::clamp(bin_z1 - pos[i].z, -ri, ri);
            if (b_sl <= a_sl) continue;
            const Real v_slice = kPi * (ri * ri * (b_sl - a_sl)
                                        - (b_sl * b_sl * b_sl - a_sl * a_sl * a_sl) / Real(3.0));
            const Real weight = vi > 0 ? (v_slice / vi) : Real(0);

            const std::size_t cell = static_cast<std::size_t>(iz) * static_cast<std::size_t>(nr)
                                     + static_cast<std::size_t>(ir);
            solid_vol[cell] += static_cast<double>(v_slice);
            virial_sum[cell] += static_cast<double>(weight * ri * fn_sum[i]);
            coord_sum[cell] += static_cast<double>(coord[i]);
            speed_sum[cell] += static_cast<double>(vel[i].norm());
            ++counts[cell];
        }
    }

    for (int iz = 0; iz < nz; ++iz) {
        for (int ir = 0; ir < nr; ++ir) {
            const std::size_t c = static_cast<std::size_t>(iz) * static_cast<std::size_t>(nr)
                                  + static_cast<std::size_t>(ir);
            const Real r0 = static_cast<Real>(ir) * dr;
            const Real r1 = static_cast<Real>(ir + 1) * dr;
            const Real ring_vol = kPi * (r1 * r1 - r0 * r0) * dz;
            grid[c].r_m = Real(0.5) * (r0 + r1);
            grid[c].z_m = box.min.z + (static_cast<Real>(iz) + Real(0.5)) * dz;
            grid[c].phi_local = ring_vol > 0
                                    ? std::min(Real(0.64), static_cast<Real>(solid_vol[c]) / ring_vol)
                                    : Real(0);
            // Isotropic compressive stress sigma = (1 / (3 V_ring)) * sum_i r_i |F_n,i| [kPa]
            grid[c].stress_kPa = ring_vol > 0
                                     ? static_cast<Real>(virial_sum[c] / (3.0 * ring_vol * 1000.0))
                                     : Real(0);
            if (counts[c] > 0) {
                grid[c].coordination = static_cast<Real>(coord_sum[c] / static_cast<double>(counts[c]));
                grid[c].speed_m_s = static_cast<Real>(speed_sum[c] / static_cast<double>(counts[c]));
            }
        }
    }
    return grid;
}

CraterObservables CraterAnalyzer::analyze(const core::ParticleSystem& particles,
                                          Real impact_energy) const {
    CraterObservables out;
    out.impact_energy = impact_energy;
    out.z_reference = z_reference_;
    out.rotational_energy = particles.rotationalKineticEnergy();
    out.impact_pressure_kPa = Real(0.5) * config_.projectile.density
                              * config_.projectile.speed * config_.projectile.speed / Real(1000.0);
    if (!has_reference_) return out;

    const core::Box box = config_.domain.box();
    const Real cell_x = box.size().x / grid_;
    const Real cell_y = box.size().y / grid_;
    const Real cell_area = cell_x * cell_y;
    const Real margin = Real(0.55) * margin_factor_ * config_.particles.r_max;
    const bool large_box = std::min(box.size().x, box.size().y) >= Real(10.0) * config_.projectile.radius;
    const Real wall_pad = large_box ? Real(2.5) * config_.particles.r_max : Real(0);
    const auto insideUsableBox = [&](int cx, int cy) {
        if (wall_pad <= 0) return true;
        const Real x = box.min.x + (static_cast<Real>(cx) + Real(0.5)) * cell_x;
        const Real y = box.min.y + (static_cast<Real>(cy) + Real(0.5)) * cell_y;
        return x >= box.min.x + wall_pad && x <= box.max.x - wall_pad
            && y >= box.min.y + wall_pad && y <= box.max.y - wall_pad;
    };

    // Smoothing radius: covers ~2 grains, the scale at which the surface stops
    // being a list of grains and becomes a surface again.
    const int smooth_radius = std::max(
        1, static_cast<int>(std::lround(Real(2.0) * config_.particles.r_max
                                        / std::min(cell_x, cell_y))));
    const auto current = smoothMap(surfaceMap(particles, /*bed_only=*/true), grid_, smooth_radius);
    const auto current_total = smoothMap(surfaceMap(particles, /*bed_only=*/false), grid_, smooth_radius);

    // Occupied bed height includes top grain radius above the top layer centers (z_reference_).
    out.bed_height = std::max<Real>(0, z_reference_ - box.min.z);
    const Real occupied_height = std::max<Real>(
        config_.particles.r_mean(), (z_reference_ + config_.particles.r_mean()) - box.min.z);
    {
        const Real volume = box.size().x * box.size().y * occupied_height;
        out.phi_bed = particles.solidsVolumeOf(ParticleKind::Bed) / volume;
    }

    // ---------- 0) Diagnostic: bed level change --------------------------------
    {
        std::vector<Real> deltas;
        deltas.reserve(reference_map_.size());
        for (std::size_t cell = 0; cell < reference_map_.size(); ++cell) {
            if (std::isnan(reference_map_[cell]) || std::isnan(current[cell])) continue;
            deltas.push_back(reference_map_[cell] - current[cell]);
        }
        if (!deltas.empty()) {
            const auto mid = deltas.begin() + static_cast<std::ptrdiff_t>(deltas.size() / 2);
            std::nth_element(deltas.begin(), mid, deltas.end());
            out.bed_shift = *mid;
        }
    }

    // ---------- 1) Bed drift correction -----------------------------------------
    const Vec3 impact{box.min.x + Real(0.5) * box.size().x, box.min.y + Real(0.5) * box.size().y, 0};
    const Real exclude = Real(4.0) * config_.projectile.radius;
    std::vector<Real> raw_delta(reference_map_.size(), std::numeric_limits<Real>::quiet_NaN());
    for (std::size_t cell = 0; cell < reference_map_.size(); ++cell) {
        if (std::isnan(reference_map_[cell]) || std::isnan(current[cell])) continue;
        raw_delta[cell] = reference_map_[cell] - current[cell];
    }

    double outer_sum = 0;
    std::size_t outer_n = 0;
    for (int cy = 0; cy < grid_; ++cy) {
        for (int cx = 0; cx < grid_; ++cx) {
            if (!insideUsableBox(cx, cy)) continue;
            const Real d = raw_delta[static_cast<std::size_t>(cy) * static_cast<std::size_t>(grid_)
                                     + static_cast<std::size_t>(cx)];
            if (std::isnan(d)) continue;
            const Real x = box.min.x + (cx + Real(0.5)) * cell_x;
            const Real y = box.min.y + (cy + Real(0.5)) * cell_y;
            const Real dx = x - impact.x;
            const Real dy = y - impact.y;
            if (std::sqrt(dx * dx + dy * dy) < exclude) continue;
            outer_sum += static_cast<double>(d);
            ++outer_n;
        }
    }
    if (outer_n < 10) {
        out.classification = "Sand Mound / no crater";
        return out;
    }
    const Real drift = static_cast<Real>(outer_sum / static_cast<double>(outer_n));

    // ---------- 2) Per-cell deformation corrected for the drift ------------------
    std::vector<Real> delta(reference_map_.size(), std::numeric_limits<Real>::quiet_NaN());
    Real d_exc_local = 0;
    for (int cy = 0; cy < grid_; ++cy) {
        for (int cx = 0; cx < grid_; ++cx) {
            if (!insideUsableBox(cx, cy)) continue;
            const auto cell = static_cast<std::size_t>(cy) * static_cast<std::size_t>(grid_)
                              + static_cast<std::size_t>(cx);
            const Real z = raw_delta[cell];      // deformation: ref - current
            if (std::isnan(z)) continue;
            // Deformation corrected for the drift: + excavated, - raised.
            delta[cell] = z - drift;
            d_exc_local = std::max(d_exc_local, delta[cell]);
        }
    }

    // ---------- 2b) Subsurface packing beneath the crater floor & Virial Stress -
    {
        const Real cx_center = box.min.x + Real(0.5) * box.size().x;
        const Real cy_center = box.min.y + Real(0.5) * box.size().y;
        const Real r_col = std::max(config_.projectile.radius, Real(4.0) * config_.particles.r_max);
        // Top of the compressed plug lies strictly below the current crater floor
        // (z_reference_ - d_exc_local), avoiding the empty air of the excavated cavity.
        const Real z_floor = std::max(box.min.z + Real(2.5) * config_.particles.r_max,
                                      z_reference_ - d_exc_local);
        const Real z_col_top = std::max(box.min.z + Real(2.0) * config_.particles.r_max,
                                        z_floor - config_.particles.r_max);
        const Real r_eff = r_col + Real(0.5) * config_.particles.r_mean();
        const Real h_eff = std::max(z_col_top - box.min.z, Real(1e-9)) + config_.particles.r_mean();
        const Real col_vol = kPi * r_eff * r_eff * h_eff;

        Real sub_solid_vol = 0;
        Real sub_virial = 0;
        Real max_grain_stress_kPa = 0;
        Real horiz_momentum = 0;

        const auto& pos = particles.position();
        const auto& vel = particles.velocity();
        const auto& rad = particles.radius();
        const auto& mass = particles.mass();
        const auto& kinds = particles.kind();
        const auto& fn_sum = particles.normalForceSum();

        for (std::size_t i = 0; i < particles.size(); ++i) {
            const Real dx = pos[i].x - cx_center;
            const Real dy = pos[i].y - cy_center;
            const Real r_xy = std::hypot(dx, dy);

            if (kinds[i] == ParticleKind::Bed) {
                const Real v_horiz = std::hypot(vel[i].x, vel[i].y);
                horiz_momentum += mass[i] * v_horiz;
            }

            if (r_xy <= r_col && pos[i].z >= box.min.z && pos[i].z <= z_col_top) {
                const Real ri = rad[i];
                const Real vi = Real(4.0 / 3.0) * kPi * ri * ri * ri;
                sub_solid_vol += vi;
                sub_virial += ri * fn_sum[i];
                if (vi > 0) {
                    const Real sigma_i_kPa = (ri * fn_sum[i]) / (Real(3.0) * vi * Real(1000.0));
                    max_grain_stress_kPa = std::max(max_grain_stress_kPa, sigma_i_kPa);
                }
            }
        }
        const Real phi_plug = col_vol > 0 ? (sub_solid_vol / col_vol) : out.phi_bed;
        // Uniaxial column compaction estimate from crater floor compression:
        const Real h0_col = std::max(z_reference_ - box.min.z, Real(1e-6));
        const Real h_comp = std::max(Real(0.35) * h0_col, h0_col - Real(0.55) * d_exc_local);
        const Real phi_uniaxial = out.phi_bed * (h0_col / h_comp);
        out.phi_subsurface = std::clamp(std::max(phi_plug, phi_uniaxial), out.phi_bed, Real(0.64));
        out.mean_virial_stress_kPa = col_vol > 0 ? sub_virial / (Real(3.0) * col_vol * Real(1000.0)) : 0;
        out.max_virial_stress_kPa = max_grain_stress_kPa;

        const Real m_proj_est = Real(4.0 / 3.0) * kPi * std::pow(config_.projectile.radius, 3)
                                * config_.projectile.density;
        const Real p0_proj = m_proj_est * std::max(config_.projectile.speed, Real(1e-6));
        out.momentum_transfer_vh = p0_proj > 0 ? horiz_momentum / p0_proj : 0;
    }

    // ---------- 3) Connected component of the excavation -------------------------
    const auto idx = [&](int cx, int cy) {
        return static_cast<std::size_t>(cy) * static_cast<std::size_t>(grid_)
               + static_cast<std::size_t>(cx);
    };
    std::size_t deepest = delta.size();
    Real deepest_delta = -std::numeric_limits<Real>::infinity();
    for (int cy = 0; cy < grid_; ++cy) {
        for (int cx = 0; cx < grid_; ++cx) {
            if (!insideUsableBox(cx, cy)) continue;
            const Real d = delta[idx(cx, cy)];
            if (!std::isnan(d) && d > margin && d > deepest_delta) {
                deepest_delta = d;
                deepest = idx(cx, cy);
            }
        }
    }

    std::size_t count = 0;
    Real sum_x = 0, sum_y = 0;
    Real v_in_integral = 0;
    std::vector<Real> xs, ys;
    Real x_min = 0, x_max = 0, y_min = 0, y_max = 0;
    std::vector<char> keep(delta.size(), 0);
    if (deepest < delta.size()) {
        std::vector<std::size_t> stack{deepest};
        keep[deepest] = 1;
        while (!stack.empty()) {
            const std::size_t cell = stack.back();
            stack.pop_back();
            const int cx = static_cast<int>(cell % static_cast<std::size_t>(grid_));
            const int cy = static_cast<int>(cell / static_cast<std::size_t>(grid_));
            const int neighbours[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
            for (const auto& nb : neighbours) {
                const int nx = cx + nb[0], ny = cy + nb[1];
                if (nx < 0 || ny < 0 || nx >= grid_ || ny >= grid_) continue;
                if (!insideUsableBox(nx, ny)) continue;
                const std::size_t nc = idx(nx, ny);
                if (keep[nc]) continue;
                const Real d = delta[nc];
                if (std::isnan(d) || d <= margin) continue;
                keep[nc] = 1;
                stack.push_back(nc);
            }
        }
        for (int cy = 0; cy < grid_; ++cy) {
            for (int cx = 0; cx < grid_; ++cx) {
                const std::size_t cell = idx(cx, cy);
                if (!keep[cell]) continue;
                d_exc_local = std::max(d_exc_local, delta[cell]);
                v_in_integral += delta[cell] * cell_area;
                const Real x = box.min.x + (cx + Real(0.5)) * cell_x;
                const Real y = box.min.y + (cy + Real(0.5)) * cell_y;
                if (count == 0) { x_min = x; x_max = x; y_min = y; y_max = y; }
                x_min = std::min(x_min, x);
                x_max = std::max(x_max, x);
                y_min = std::min(y_min, y);
                y_max = std::max(y_max, y);
                xs.push_back(x);
                ys.push_back(y);
                sum_x += x;
                sum_y += y;
                ++count;
            }
        }
    }

    out.crater_cells = count;

    // Ejecta (always reported, whether or not there is a cavity).
    std::size_t n_ejecta = 0;
    Real speed_sum = 0;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (particles.kind()[i] != ParticleKind::Ejecta) continue;
        out.ejecta_mass += particles.mass()[i];
        speed_sum += particles.velocity()[i].norm();
        ++n_ejecta;
    }
    out.ejecta_speed_mean = n_ejecta > 0 ? speed_sum / static_cast<Real>(n_ejecta) : 0;

    // Analytical models can always be evaluated from bed and impact parameters.
    const Real bed_density = config_.particles.rho_grain * std::max(out.phi_bed, Real(1e-3));
    const Real dyn_pressure_Pa = out.impact_pressure_kPa * Real(1000.0);
    const Real h_column = std::min(out.bed_height, Real(2.5) * config_.projectile.radius);
    out.d_log_model = models::diameterLogModel(out.impact_energy);
    out.d_uehara_model = models::diameterUeharaModel(out.impact_energy, bed_density,
                                                     config_.particles.rho_grain,
                                                     2 * config_.projectile.radius,
                                                     std::max(config_.particles.friction, Real(0.1)),
                                                     std::max(std::abs(config_.gravity.gz), Real(0.1)));
    out.d_heckel_model = models::depthHeckelUniaxialModel(h_column, out.phi_bed, dyn_pressure_Pa);
    out.h_rim_model = models::rimHeightLogLogModel(out.impact_energy);

    if (count < 3) {
        out.classification = "Sand Mound / no crater";
        return out;
    }

    out.d_exc = d_exc_local;
    out.v_in_num = v_in_integral;

    // ---------- 4) Cavity ellipse by moments ------------------------------------
    const Real mean_x = sum_x / static_cast<Real>(count);
    const Real mean_y = sum_y / static_cast<Real>(count);
    Real var_xx = 0, var_yy = 0, cov_xy = 0;
    for (std::size_t k = 0; k < xs.size(); ++k) {
        const Real dx = xs[k] - mean_x;
        const Real dy = ys[k] - mean_y;
        var_xx += dx * dx;
        var_yy += dy * dy;
        cov_xy += dx * dy;
    }
    var_xx /= static_cast<Real>(count);
    var_yy /= static_cast<Real>(count);
    cov_xy /= static_cast<Real>(count);

    const Real trace = var_xx + var_yy;
    const Real diff = var_xx - var_yy;
    const Real root = std::sqrt(std::max<Real>(0, diff * diff / 4 + cov_xy * cov_xy));
    const Real lambda_major = trace / 2 + root;
    const Real lambda_minor = std::max<Real>(0, trace / 2 - root);
    const Real a = 2 * std::sqrt(lambda_major);
    const Real b_axis = 2 * std::sqrt(lambda_minor);
    out.D = 2 * a;
    out.d_minor = 2 * b_axis;
    out.theta_rad = Real(0.5) * std::atan2(2 * cov_xy, diff);
    out.epsilon = a > 0 ? std::sqrt(std::max<Real>(0, 1 - (b_axis * b_axis) / (a * a))) : 0;

    // ---------- 5) Rim height, raised volume V_ex, inner slope, and central peak -
    Real h_rim_local = 0;
    Real v_ex_integral = 0;
    Real slope_sum = 0;
    std::size_t slope_count = 0;
    const Real rim_margin = std::min(Real(0.0012), Real(0.30) * margin);
    const Real max_rim_radius = std::max(Real(1.85) * a, Real(2.2) * config_.projectile.radius);
    {
        const Real pad = std::max(Real(3.0) * std::max(cell_x, cell_y), Real(0.5) * a);
        const int cx_lo = static_cast<int>((x_min - pad - box.min.x) / cell_x);
        const int cx_hi = static_cast<int>((x_max + pad - box.min.x) / cell_x);
        const int cy_lo = static_cast<int>((y_min - pad - box.min.y) / cell_y);
        const int cy_hi = static_cast<int>((y_max + pad - box.min.y) / cell_y);
        for (int cy = std::max(0, cy_lo); cy <= std::min(grid_ - 1, cy_hi); ++cy) {
            for (int cx = std::max(0, cx_lo); cx <= std::min(grid_ - 1, cx_hi); ++cx) {
                if (!insideUsableBox(cx, cy)) continue;
                const std::size_t c = idx(cx, cy);
                const Real d = delta[c];
                if (std::isnan(d)) continue;
                const Real x = box.min.x + (cx + Real(0.5)) * cell_x;
                const Real y = box.min.y + (cy + Real(0.5)) * cell_y;
                const Real r_c = std::hypot(x - mean_x, y - mean_y);
                if (-d > rim_margin && r_c <= max_rim_radius) {
                    h_rim_local = std::max(h_rim_local, -d);
                    v_ex_integral += (-d) * cell_area;
                }
                // Inner crater wall slope |grad z| inside the outer half of the cavity
                if (keep[c] && cx > 0 && cx + 1 < grid_ && cy > 0 && cy + 1 < grid_) {
                    if (r_c >= Real(0.45) * a && r_c <= Real(1.05) * a) {
                        const Real dx_val = delta[idx(cx + 1, cy)] - delta[idx(cx - 1, cy)];
                        const Real dy_val = delta[idx(cx, cy + 1)] - delta[idx(cx, cy - 1)];
                        if (!std::isnan(dx_val) && !std::isnan(dy_val)) {
                            const Real gx = dx_val / (Real(2.0) * cell_x);
                            const Real gy = dy_val / (Real(2.0) * cell_y);
                            slope_sum += std::hypot(gx, gy);
                            ++slope_count;
                        }
                    }
                }
            }
        }
    }

    // Central peak ("porous piston" of projectile grains resting in the cavity center)
    Real h_cp_local = 0;
    Real v_cp_local = 0;
    for (int cy = 0; cy < grid_; ++cy) {
        for (int cx = 0; cx < grid_; ++cx) {
            const std::size_t c = idx(cx, cy);
            if (std::isnan(current[c]) || std::isnan(current_total[c])) continue;
            const Real x = box.min.x + (cx + Real(0.5)) * cell_x;
            const Real y = box.min.y + (cy + Real(0.5)) * cell_y;
            if (std::hypot(x - mean_x, y - mean_y) <= Real(0.55) * a) {
                const Real dh = current_total[c] - current[c];
                if (dh > margin) {
                    h_cp_local = std::max(h_cp_local, dh);
                    v_cp_local += dh * cell_area;
                }
            }
        }
    }

    const Real max_physical_rim = std::max(Real(6.0) * config_.particles.r_max,
                                           Real(0.4) * config_.projectile.radius);
    out.h_rim = std::min(h_rim_local, max_physical_rim);
    out.d_max = out.h_rim + out.d_exc;
    out.rim_slope = slope_count > 0 ? slope_sum / static_cast<Real>(slope_count) : 0;
    out.aspect_zd = out.D > 0 ? out.d_exc / out.D : 0;
    out.v_ex = v_ex_integral;
    out.v_exc = out.v_in_num - out.v_ex;
    out.h_cp = h_cp_local;
    out.v_cp = v_cp_local;

    // ---------- 6) Cavity volume ------------------------------------------------
    out.v_in = (core::kPi / 3) * a * a * out.d_exc;

    // ---------- 7) Morphological classification ---------------------------------
    if (out.h_cp > Real(2.0) * config_.particles.r_max && out.aspect_zd > 0.08) {
        out.classification = "Complex (central peak)";
    } else if (out.aspect_zd > 0.20) {
        out.classification = "Simple (deep)";
    } else if (out.aspect_zd > 0.05) {
        out.classification = "Simple";
    } else {
        out.classification = "Sand Mound";
    }
    return out;
}

}  // namespace granimpact::analysis
