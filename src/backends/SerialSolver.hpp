#pragma once
 // Reference solver: one thread, no dependencies. It is the yardstick against
 // which the accelerated stages are validated (same result, different time).
#include "granimpact/sim/ContactSolvers.hpp"

namespace granimpact::backends {

class CpuContactSolver final : public ContactSolver {
public:
    explicit CpuContactSolver(const core::SimConfig& config) : model_(config) {}
    [[nodiscard]] std::string name() const override { return "serial"; }
    physics::ContactStats computeForces(core::ParticleSystem& particles,
                                        const physics::CellList& cells, core::Real dt) override {
        return model_.computeForces(particles, cells, dt);
    }

private:
    physics::ContactModel model_;
};

}  // namespace granimpact::backends
