// The cell list is the part of the code with the highest risk of a silent bug:
// skip one pair and the simulation 'works' but the physics is wrong.
#include <gtest/gtest.h>

#include <random>
#include <set>
#include <vector>

#include "granimpact/core/ParticleSystem.hpp"
#include "granimpact/physics/CellList.hpp"

using namespace granimpact;

namespace {

core::ParticleSystem randomCloud(std::size_t n, core::Real box_side, core::Real radius,
                                std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<core::Real> uni(radius, box_side - radius);
    core::ParticleSystem ps(n);
    for (std::size_t i = 0; i < n; ++i) {
        ps.add({uni(rng), uni(rng), uni(rng)}, radius, 2650.0, core::ParticleKind::Bed);
    }
    ps.setMaterial(70e9, 0.3, 0.4, 0.5);
    return ps;
}

}  // namespace

TEST(CellList, CadaParSeVisitaUnaSolaVez) {
    const core::Real side = 0.05;
    auto ps = randomCloud(200, side, 2e-3, 7);
    const core::Box box{{0, 0, 0}, {side, side, side}};

    physics::CellList cells;
    cells.build(ps, box, physics::CellList::suggestCellSize(2e-3));

    std::set<std::pair<int, int>> seen;
    std::size_t count = 0;
    cells.forEachPair([&](core::Index i, core::Index j) {
        ASSERT_LT(i, j) << "pairs must come ordered and without repetition";
        const auto inserted = seen.insert({i, j}).second;
        EXPECT_TRUE(inserted) << "par repetido: " << i << "," << j;
        ++count;
    });
    EXPECT_EQ(count, seen.size());
}

TEST(CellList, EncuentraTodosLosParesDentroDelRadioDeCorte) {
    const core::Real side = 0.06;
    const core::Real radius = 1.5e-3;
    auto ps = randomCloud(600, side, radius, 4242);
    const core::Box box{{0, 0, 0}, {side, side, side}};

    physics::CellList cells;
    cells.build(ps, box, physics::CellList::suggestCellSize(radius));

    std::set<std::pair<int, int>> found;
    cells.forEachPair([&](core::Index i, core::Index j) { found.insert({i, j}); });

    std::size_t expected = 0;
    for (std::size_t i = 0; i < ps.size(); ++i) {
        for (std::size_t j = i + 1; j < ps.size(); ++j) {
            const core::Real dist = (ps.position()[i] - ps.position()[j]).norm();
            if (dist < 2 * radius) {
                ++expected;
                EXPECT_TRUE(found.count({static_cast<int>(i), static_cast<int>(j)}) == 1)
                    << "contact pair not detected: " << i << "," << j;
            }
        }
    }
    EXPECT_GT(expected, 0u);
}

TEST(CellList, ElTamanoDeCeldaCubreElDiametroMaximo) {
    EXPECT_NEAR(physics::CellList::suggestCellSize(1e-3), 2e-3 * 1.05, 1e-15);
    EXPECT_GT(physics::CellList::suggestCellSize(1e-3), 2e-3);
}

TEST(CellList, ReconstruyeLaGeometriaDeLasCeldas) {
    const core::Real side = 0.1;
    auto ps = randomCloud(50, side, 2e-3, 11);
    const core::Box box{{0, 0, 0}, {side, side, side}};

    physics::CellList cells;
    cells.build(ps, box, 0.02);   // 5 x 5 x 5 celdas
    EXPECT_EQ(cells.dims()[0], 5);
    EXPECT_EQ(cells.dims()[1], 5);
    EXPECT_EQ(cells.dims()[2], 5);
    EXPECT_EQ(cells.cellCount(), 125u);
    EXPECT_EQ(cells.ordered().size(), ps.size());
}
