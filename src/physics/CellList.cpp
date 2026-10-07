#include "granimpact/physics/CellList.hpp"

#include <algorithm>
#include <cmath>

namespace granimpact::physics {

using core::Real;

const std::vector<CellList::Offset>& CellList::halfOffsets() {
    static const std::vector<Offset> offsets = {
        {0, 0, 0},
        {1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {1, 1, 0},
        {-1, -1, 1}, {0, -1, 1}, {1, -1, 1},
        {-1, 0, 1}, {0, 0, 1}, {1, 0, 1},
        {-1, 1, 1}, {0, 1, 1}, {1, 1, 1},
    };
    return offsets;
}

Real CellList::suggestCellSize(Real max_radius, Real safety) {
    // Classic rule: the cell must cover the maximum diameter, with margin.
    return 2.0 * max_radius * safety;
}

void CellList::build(const core::ParticleSystem& particles, const core::Box& box, Real cell_size) {
    if (cell_size <= 0) throw std::invalid_argument("CellList: cell_size must be > 0");

    const core::Vec3 size = box.max - box.min;
    origin_ = box.min;
    inv_cell_ = {1.0 / cell_size, 1.0 / cell_size, 1.0 / cell_size};
    dims_ = {std::max(1, static_cast<int>(std::floor(size.x / cell_size))),
             std::max(1, static_cast<int>(std::floor(size.y / cell_size))),
             std::max(1, static_cast<int>(std::floor(size.z / cell_size)))};

    const std::size_t n = particles.size();
    ordered_.resize(n);
    counts_.assign(cellCount() + 1, 0);

    // 1) count per cell
    std::vector<std::uint32_t> cell_of(n);
    for (std::size_t i = 0; i < n; ++i) {
        const core::Vec3 p = particles.position()[i];
        const int cx = std::clamp(static_cast<int>((p.x - origin_.x) * inv_cell_.x), 0, dims_[0] - 1);
        const int cy = std::clamp(static_cast<int>((p.y - origin_.y) * inv_cell_.y), 0, dims_[1] - 1);
        const int cz = std::clamp(static_cast<int>((p.z - origin_.z) * inv_cell_.z), 0, dims_[2] - 1);
        const std::uint32_t c = index(cx, cy, cz);
        cell_of[i] = c;
        ++counts_[c + 1];
    }

    // 2) prefijos (CSR)
    for (std::size_t c = 1; c < counts_.size(); ++c) counts_[c] += counts_[c - 1];
    cell_start_.assign(counts_.begin(), counts_.end());

    // 3) place each particle in its cell
    std::vector<std::uint32_t> cursor(cell_start_.begin(), cell_start_.end() - 1);
    for (std::size_t i = 0; i < n; ++i) {
        ordered_[cursor[cell_of[i]]++] = static_cast<core::Index>(i);
    }
}

}  // namespace granimpact::physics
