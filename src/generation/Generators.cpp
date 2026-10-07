#include "granimpact/generation/Generators.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <array>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace granimpact::generation {

using core::kPi;
using core::ParticleKind;
using core::Real;
using core::Vec3;

namespace {

 // Hash grid for O(1) neighbour queries: used when placing each grain to
 // guarantee that none overlaps the ones already placed.
class NeighborHash {
public:
    explicit NeighborHash(Real cell) : cell_(cell) {}

    void insert(const Vec3& p, Real radius) {
        const auto key = keyOf(p);
        buckets_[key].push_back({p, radius});
    }

    // Maximum radius allowed at `p` without overlapping any neighbour (clearance aside).
    [[nodiscard]] Real maxAllowedRadius(const Vec3& p, Real clearance) const {
        Real allowed = std::numeric_limits<Real>::max();
        for (int dz = -1; dz <= 1; ++dz) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dx = -1; dx <= 1; ++dx) {
                    const auto key = keyOf(p + Vec3{dx * cell_, dy * cell_, dz * cell_});
                    const auto it = buckets_.find(key);
                    if (it == buckets_.end()) continue;
                    for (const auto& [q, rq] : it->second) {
                        const Real dist = (q - p).norm();
                        allowed = std::min(allowed, dist - rq - clearance);
                    }
                }
            }
        }
        return allowed;
    }

private:
    [[nodiscard]] std::int64_t keyOf(const Vec3& p) const {
        const auto i = static_cast<std::int64_t>(std::floor(p.x / cell_));
        const auto j = static_cast<std::int64_t>(std::floor(p.y / cell_));
        const auto k = static_cast<std::int64_t>(std::floor(p.z / cell_));
        return ((i * 73856093) ^ (j * 19349663) ^ (k * 83492791));
    }

    Real cell_;
    std::unordered_map<std::int64_t, std::vector<std::pair<Vec3, Real>>> buckets_;
};

}  // namespace

Real measuredPacking(const core::ParticleSystem& particles, Real volume) {
    if (volume <= 0) return 0;
    return particles.solidsVolume() / volume;
}

BedInfo generateBed(core::ParticleSystem& particles, const core::SimConfig& config,
                    std::mt19937_64& rng) {
    const auto& pc = config.particles;
    const core::Box box = config.domain.box();

    BedInfo info;
    info.target_phi = pc.packing_fraction;
    info.region_height = pc.bed_height_fraction * box.size().z;

    // =========================== Bed method ============================
    // A loose granular bed CANNOT be seeded as an empty lattice: without contacts it is
    // not in equilibrium and collapses as soon as the simulation starts (the
    // subsidence would then be read as a false crater). The method that works:
    // 1) seed a DENSE nested lattice (layers shifted by half a cell) with grains
    // in contact -> a jammed configuration, disordered by the jitter;
    // 2) EXPAND it isotropically to the target density -> the grains end up almost
    // touching and the contact topology is preserved;
    // 3) let settling readjust it locally (see Simulation).
    // Result: stable bed, no overlaps, at the requested density and reproducible.
    // The real final density is MEASURED and reported (never assumed).
    // ==========================================================================
    const Real s = Real(2.01) * pc.r_max;            // horizontal spacing (no overlaps)
    const Real dz_dense = Real(1.42) * pc.r_max;     // anidado: contacto capa a capa

    // Mean grain volume: (4/3) pi E[r^3], with r uniform in [r_min, r_max].
    const Real e_r3 = (std::pow(pc.r_max, 4) - std::pow(pc.r_min, 4))
                      / (4 * (pc.r_max - pc.r_min));
    const Real v_mean = Real(4.0 / 3.0) * kPi * e_r3;

    const Real phi_dense = v_mean / (s * s * dz_dense);          // ~0.6-0.75
    // If the dense packing the lattice allows is already LESS dense than the target
    // (which happens with wide radius distributions), it is not compressed: compressing
    // would create overlaps. It is seeded as is and the density obtained is reported.
    const Real expansion = std::max(Real(1.0), std::cbrt(phi_dense / pc.packing_fraction));
    info.lattice_spacing = s * expansion;
    info.layer_spacing = dz_dense * expansion;

    // Region where the dense lattice is built: after expanding by `expansion` it must
    // fit exactly in the usable box.
    const Real dense_x = (box.size().x - 2 * pc.r_max) / expansion;
    const Real dense_y = (box.size().y - 2 * pc.r_max) / expansion;
    const Real dense_z = (info.region_height - 2 * pc.r_max) / expansion;

    const int nx = std::max(1, static_cast<int>(std::floor(dense_x / s)));
    const int ny = std::max(1, static_cast<int>(std::floor(dense_y / s)));
    const int nz = std::max(1, static_cast<int>(std::floor(dense_z / dz_dense)));

    const Vec3 center{box.min.x + Real(0.5) * box.size().x, box.min.y + Real(0.5) * box.size().y,
                      box.min.z};
    const auto expand = [&](const Vec3& p) {
        return Vec3{center.x + (p.x - center.x) * expansion,
                    center.y + (p.y - center.y) * expansion,
                    center.z + (p.z - center.z) * expansion};
    };

    NeighborHash hash(Real(2.05) * pc.r_max * expansion);
    std::uniform_real_distribution<Real> uni(Real(0), Real(1));
    const Real jitter = Real(0.08) * s;
    const Real clearance = Real(0.01) * pc.r_mean();

    const auto place = [&](const Vec3& raw_position, Real radius) {
        const Vec3 p = expand(raw_position);
        const auto index = particles.add(p, radius, pc.rho_grain, ParticleKind::Bed);
        particles.setMaterialOf(index, pc.young, pc.poisson, pc.friction, pc.restitution);
        hash.insert(p, radius);
    };

    for (int iz = 0; iz < nz; ++iz) {
        const Real off = (iz % 2 == 0) ? Real(0) : Real(0.5) * s;
        for (int iy = 0; iy < ny; ++iy) {
            for (int ix = 0; ix < nx; ++ix) {
                const Vec3 raw{center.x - Real(0.5) * dense_x + (ix + Real(0.5)) * s + off
                                   + jitter * (uni(rng) - Real(0.5)),
                               center.y - Real(0.5) * dense_y + (iy + Real(0.5)) * s + off
                                   + jitter * (uni(rng) - Real(0.5)),
                               box.min.z + pc.r_max + (iz + Real(0.5)) * dz_dense
                                   + Real(0.3) * jitter * (uni(rng) - Real(0.5))};
                Real radius = pc.r_min + uni(rng) * (pc.r_max - pc.r_min);
                // The grain is placed already expanded, so the overlap
                // check is done at the final position.
                radius = std::min(radius, hash.maxAllowedRadius(expand(raw), clearance));
                if (radius < pc.r_min) continue;
                place(raw, radius);
                ++info.count;
            }
        }
    }

    // Seeded density: measured over the LATTICE VOLUME (nx*s x ny*s x nz*dz), not over
    // the full region. The lattice is a block with an edge of radii that is not
    // filled; averaging that with the region would give a falsely low density. After
    // settling the bed spreads out and then it is measured over the real occupied volume.
    const Real envelope = static_cast<Real>(nx) * info.lattice_spacing
                          * (static_cast<Real>(ny) * info.lattice_spacing)
                          * (static_cast<Real>(nz) * info.layer_spacing);
    info.seeded_phi = measuredPacking(particles, envelope);
    return info;
}

std::size_t generateProjectile(core::ParticleSystem& particles, const core::SimConfig& config,
                               std::mt19937_64& rng) {
    const auto& pj = config.projectile;
    const auto& pc = config.particles;
    const core::Box box = config.domain.box();

    const Vec3 center{box.min.x + Real(0.5) * box.size().x, box.min.y + Real(0.5) * box.size().y,
                      pj.start_z};
    const Real theta = pj.angle_deg * kPi / 180.0;
    const Vec3 velocity{pj.speed * std::cos(theta), 0.0, -pj.speed * std::sin(theta)};

    const Real s = Real(2.02) * pc.r_max;
    const Real dz = Real(0.85) * s;
    const int steps = static_cast<int>(std::ceil(pj.radius / s)) + 1;
    const Real jitter = Real(0.08) * s;
    const Real clearance = Real(0.02) * pc.r_mean();

    NeighborHash hash(Real(2.05) * pc.r_max);
    std::uniform_real_distribution<Real> uni(Real(0), Real(1));

    std::size_t placed = 0;
    for (int iz = -steps; iz <= steps; ++iz) {
        const Real off = (iz % 2 == 0) ? Real(0) : Real(0.5) * s;
        for (int iy = -steps; iy <= steps; ++iy) {
            for (int ix = -steps; ix <= steps; ++ix) {
                const Vec3 offset{(ix + Real(0.5)) * s + off + jitter * (uni(rng) - Real(0.5)),
                                  (iy + Real(0.5)) * s + off + jitter * (uni(rng) - Real(0.5)),
                                  (iz + Real(0.5)) * dz + jitter * (uni(rng) - Real(0.5))};
                // The grain must fit inside the projectile sphere.
                const Real dist = offset.norm();
                if (dist > pj.radius - pc.r_max) continue;

                const Vec3 position = center + offset;
                Real radius = pc.r_min + uni(rng) * (pc.r_max - pc.r_min);
                radius = std::min(radius, pj.radius - dist);
                radius = std::min(radius, hash.maxAllowedRadius(position, clearance));
                if (radius < pc.r_min) continue;

                const auto index = particles.add(position, radius, pj.density,
                                                 ParticleKind::Projectile);
                particles.setMaterialOf(index, pj.young, pj.poisson, pc.friction,
                                        pc.restitution);
                particles.velocity()[static_cast<std::size_t>(index)] = velocity;
                hash.insert(position, radius);
                ++placed;
            }
        }
    }
    return placed;
}

}  // namespace granimpact::generation
