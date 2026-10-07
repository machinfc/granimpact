#include "granimpact/analysis/CraterAnalyzer.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

namespace granimpact::analysis {

using core::ParticleKind;
using core::Real;
using core::Vec3;

namespace models {

Real diameterLogModel(Real energy_joules) {
    if (energy_joules <= 1.0) return 0;   // the log law was calibrated for E > 1 J
    return 0.089 * std::log10(energy_joules);
}

Real diameterUeharaModel(Real energy_joules, Real bed_density, Real grain_density,
                         Real projectile_diameter) {
    if (energy_joules <= 0 || grain_density <= 0) return 0;
    const Real ratio = bed_density / grain_density;
    return 1.84 * std::pow(ratio, Real(0.25)) * std::pow(projectile_diameter, Real(0.75))
           * std::pow(energy_joules, Real(0.25));
}

Real depthHeckelModel(Real h, Real phi0, Real phi_pressure) {
    if (phi_pressure <= 0) return 0;
    const Real ratio = std::clamp(phi0 / phi_pressure, Real(0), Real(1));
    return h * (1.0 - ratio);
}

}  // namespace models

namespace {

 // Spatial smoothing of a surface map. With coarse grains, one grid cell holds
 // one or two grains, so the maximum z per cell jumps a whole diameter whenever a
 // grain moves. Averaging over a neighbourhood of a few grains removes that
 // granularity noise and leaves the crater shape.
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

    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (bed_only && particles.kind()[i] != ParticleKind::Bed) continue;
        const Vec3 p = particles.position()[i];
        const auto cx = static_cast<int>((p.x - box.min.x) / box.size().x * grid_);
        const auto cy = static_cast<int>((p.y - box.min.y) / box.size().y * grid_);
        if (cx < 0 || cy < 0 || cx >= grid_ || cy >= grid_) continue;
        Real& cell = map[static_cast<std::size_t>(cy) * static_cast<std::size_t>(grid_)
                         + static_cast<std::size_t>(cx)];
        // Maximum z: that is the surface a depth camera would see.
        cell = std::isnan(cell) ? p.z : std::max(cell, p.z);
    }
    return map;
}

void CraterAnalyzer::recordReferenceSurface(const core::ParticleSystem& particles) {
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

CraterObservables CraterAnalyzer::analyze(const core::ParticleSystem& particles,
                                          Real impact_energy) const {
    CraterObservables out;
    out.impact_energy = impact_energy;
    out.z_reference = z_reference_;
    if (!has_reference_) return out;

    const core::Box box = config_.domain.box();
    const Real cell_x = box.size().x / grid_;
    const Real cell_y = box.size().y / grid_;
    const Real margin = margin_factor_ * config_.particles.r_max;

    // Smoothing radius: covers ~2 grains, the scale at which the surface stops
    // being a list of grains and becomes a surface again.
    const int smooth_radius = std::max(
        1, static_cast<int>(std::lround(Real(2.0) * config_.particles.r_max
                                        / std::min(cell_x, cell_y))));
    const auto current = smoothMap(surfaceMap(particles, /*bed_only=*/true), grid_, smooth_radius);

    // ---------- 0) Diagnostic: bed level change --------------------------------
    // Median deformation over the whole domain: how much the bed sank (or rose)
    // during the impact phase. Reported as is, so that it is visible whether the case
    // has a large drift (and therefore whether the step-1 correction matters).
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
    // The deformation (reference - current) mixes two things: the crater and the slow
    // drift of the bed, which keeps settling a few tenths of a grain during the impact
    // phase. The drift is measured on the outer ring (outside the impact zone) and
    // subtracted.
    const Vec3 impact{box.min.x + Real(0.5) * box.size().x, box.min.y + Real(0.5) * box.size().y, 0};
    const Real exclude = Real(4.0) * config_.projectile.radius;
    // Raw per-cell deformation: reference (settled bed) minus current state.
    std::vector<Real> raw_delta(reference_map_.size(), std::numeric_limits<Real>::quiet_NaN());
    for (std::size_t cell = 0; cell < reference_map_.size(); ++cell) {
        if (std::isnan(reference_map_[cell]) || std::isnan(current[cell])) continue;
        raw_delta[cell] = reference_map_[cell] - current[cell];
    }

    // Drift correction: mean deformation on the outer ring, that is, far from the
    // impact. Only the LEVEL is corrected (not the shape): fitting a plane to the
    // ring and evaluating it at the centre would extrapolate beyond the data, and the
    // settled surface is not flat.
    double outer_sum = 0;
    std::size_t outer_n = 0;
    for (int cy = 0; cy < grid_; ++cy) {
        for (int cx = 0; cx < grid_; ++cx) {
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
        out.bed_height = std::max<Real>(0, z_reference_ - box.min.z);
        const Real volume = box.size().x * box.size().y * std::max(out.bed_height, Real(1e-12));
        out.phi_bed = particles.solidsVolume() / volume;
        return out;
    }
    const Real drift = static_cast<Real>(outer_sum / static_cast<double>(outer_n));

    // ---------- 2) Per-cell deformation corrected for the drift ------------------
    std::vector<Real> delta(reference_map_.size(), std::numeric_limits<Real>::quiet_NaN());
    Real d_exc_local = 0;
    for (int cy = 0; cy < grid_; ++cy) {
        for (int cx = 0; cx < grid_; ++cx) {
            const auto cell = static_cast<std::size_t>(cy) * static_cast<std::size_t>(grid_)
                              + static_cast<std::size_t>(cx);
            const Real z = raw_delta[cell];      // deformation: ref - current
            if (std::isnan(z)) continue;
            // Deformation corrected for the drift: + excavated, - raised.
            delta[cell] = z - drift;
            d_exc_local = std::max(d_exc_local, delta[cell]);
        }
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
            const Real d = delta[idx(cx, cy)];
            if (!std::isnan(d) && d > margin && d > deepest_delta) {
                deepest_delta = d;
                deepest = idx(cx, cy);
            }
        }
    }

    std::size_t count = 0;
    Real sum_x = 0, sum_y = 0;
    std::vector<Real> xs, ys;
    Real x_min = 0, x_max = 0, y_min = 0, y_max = 0;
    if (deepest < delta.size()) {
        std::vector<char> keep(delta.size(), 0);
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
    out.bed_height = std::max<Real>(0, z_reference_ - box.min.z);
    {
        const Real volume = box.size().x * box.size().y * std::max(out.bed_height, Real(1e-12));
        out.phi_bed = particles.solidsVolume() / volume;
    }

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

    if (count < 3) {
        out.classification = "Sand Mound / no crater";
        return out;
    }

    out.d_exc = d_exc_local;

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

    // ---------- 5) Rim: maximum elevation above the plane in the crater box -----
    Real h_rim_local = 0;
    {
        const Real pad = 2.0 * std::max(cell_x, cell_y);
        const int cx_lo = static_cast<int>((x_min - pad - box.min.x) / cell_x);
        const int cx_hi = static_cast<int>((x_max + pad - box.min.x) / cell_x);
        const int cy_lo = static_cast<int>((y_min - pad - box.min.y) / cell_y);
        const int cy_hi = static_cast<int>((y_max + pad - box.min.y) / cell_y);
        for (int cy = std::max(0, cy_lo); cy <= std::min(grid_ - 1, cy_hi); ++cy) {
            for (int cx = std::max(0, cx_lo); cx <= std::min(grid_ - 1, cx_hi); ++cx) {
                const Real d = delta[idx(cx, cy)];
                if (std::isnan(d)) continue;
                h_rim_local = std::max(h_rim_local, -d);
            }
        }
    }
    // A real rim cannot exceed a couple of grain layers: higher elevations
    // are grains in flight, not relief.
    out.h_rim = std::min(h_rim_local, Real(2.0) * config_.particles.r_max);
    out.d_max = out.h_rim + out.d_exc;
    out.aspect_zd = out.D > 0 ? out.d_exc / out.D : 0;

    // ---------- 6) Cavity volume and models -------------------------------------
    out.v_in = (core::kPi / 3) * a * a * out.d_exc;
    const Real bed_density = config_.particles.rho_grain * std::max(out.phi_bed, Real(1e-3));
    out.d_log_model = models::diameterLogModel(out.impact_energy);
    out.d_uehara_model = models::diameterUeharaModel(out.impact_energy, bed_density,
                                                     config_.particles.rho_grain,
                                                     2 * config_.projectile.radius);
    out.d_heckel_model = models::depthHeckelModel(
        out.bed_height, out.phi_bed, std::min<Real>(Real(0.64), out.phi_bed * Real(1.15)));

    // ---------- 7) Morphological classification ---------------------------------
    if (out.aspect_zd > 0.20) {
        out.classification = "Simple (deep)";
    } else if (out.aspect_zd > 0.05) {
        out.classification = "Simple";
    } else {
        out.classification = "Sand Mound";
    }
    return out;
}

}  // namespace granimpact::analysis
