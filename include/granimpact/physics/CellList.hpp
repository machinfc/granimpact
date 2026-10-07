#pragma once
// Neighbour search with a linked cell list: O(N) instead of O(N^2).
// Implemented with flat arrays (CSR) instead of pointer-linked lists,
// because that layout is the one that maps 1:1 to GPU memory later.
#include <array>
#include <cstdint>
#include <vector>

#include "granimpact/core/ParticleSystem.hpp"
#include "granimpact/core/Types.hpp"

namespace granimpact::physics {

// Local aliases: they keep signatures readable without losing the core's single type.
using core::Index;
using core::Real;
using core::Vec3;

class CellList {
public:
    using Offset = std::array<int, 3>;

    // Safe cell size: it must be >= 2*r_max so that no contact is missed.
    static Real suggestCellSize(Real max_radius, Real safety = 1.05);

    void build(const core::ParticleSystem& particles, const core::Box& box, Real cell_size);

    // Walks every candidate pair (i<j) exactly once.
    template <typename Fn>
    void forEachPair(Fn&& fn) const {
        for (int cz = 0; cz < dims_[2]; ++cz) {
            for (int cy = 0; cy < dims_[1]; ++cy) {
                for (int cx = 0; cx < dims_[0]; ++cx) {
                    const std::uint32_t c = index(cx, cy, cz);
                    for (const Offset& o : halfOffsets()) {
                        const int nx = cx + o[0], ny = cy + o[1], nz = cz + o[2];
                        if (!inside(nx, ny, nz)) continue;
                        const std::uint32_t n = index(nx, ny, nz);
                        if (n < c) continue;                 // only half the neighbourhood
                        for (std::uint32_t a = begin(c); a < begin(c + 1); ++a) {
                            const std::uint32_t b0 = (n == c) ? a + 1 : begin(n);
                            for (std::uint32_t b = b0; b < begin(n + 1); ++b) {
                                // The pair is always emitted as (lower, higher): the
                                // contract is explicit for whoever consumes pairs.
                                const core::Index pi = ordered_[a];
                                const core::Index pj = ordered_[b];
                                if (pi < pj) fn(pi, pj); else fn(pj, pi);
                            }
                        }
                    }
                }
            }
        }
    }

    [[nodiscard]] std::uint32_t cellCount() const {
        return static_cast<std::uint32_t>(dims_[0] * dims_[1] * dims_[2]);
    }
    [[nodiscard]] const std::array<int, 3>& dims() const { return dims_; }
    [[nodiscard]] const std::vector<core::Index>& ordered() const { return ordered_; }

private:
    // 13 offsets + the cell itself = the 14 needed (half the neighbourhood).
    static const std::vector<Offset>& halfOffsets();

    [[nodiscard]] std::uint32_t index(int cx, int cy, int cz) const {
        return static_cast<std::uint32_t>((cz * dims_[1] + cy) * dims_[0] + cx);
    }
    [[nodiscard]] bool inside(int cx, int cy, int cz) const {
        return cx >= 0 && cy >= 0 && cz >= 0 && cx < dims_[0] && cy < dims_[1] && cz < dims_[2];
    }
    [[nodiscard]] std::uint32_t begin(std::uint32_t cell) const {
        return cell_start_[cell];
    }

    core::Vec3 origin_{};
    core::Vec3 inv_cell_{};
    std::array<int, 3> dims_{1, 1, 1};
    std::vector<std::uint32_t> cell_start_;   // CSR: start of each cell
    std::vector<core::Index> ordered_;        // indices grouped by cell
    std::vector<std::uint32_t> counts_;
};

}  // namespace granimpact::physics
