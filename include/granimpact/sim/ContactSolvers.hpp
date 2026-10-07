#pragma once
 // Contact solver interface: the ONLY piece of the simulator that gets accelerated.
//
 // The core (generation, integration, analysis, outputs) does not know which backend
 // runs underneath; it just asks to 'compute the contact forces' and receives stats.
 // That is why stages 4 (OpenMP) and 5 (CUDA) are added without touching anything
// docs/architecture.md.
#include <memory>
#include <string>
#include <vector>

#include "granimpact/core/ParticleSystem.hpp"
#include "granimpact/core/SimConfig.hpp"
#include "granimpact/physics/CellList.hpp"
#include "granimpact/physics/ContactModel.hpp"

namespace granimpact::backends {

 // Interface satisfied by every contact solver (CPU and GPU).
class ContactSolver {
public:
    virtual ~ContactSolver() = default;
    [[nodiscard]] virtual std::string name() const = 0;
    virtual physics::ContactStats computeForces(core::ParticleSystem& particles,
                                                const physics::CellList& cells, core::Real dt) = 0;
};

struct BackendStatus {
    std::string name;
    std::string description;
    bool available{false};
};

#ifdef GRANIMPACT_ENABLE_OPENMP
// Defined in src/backends/OpenMpSolver.cpp (stage 4).
std::unique_ptr<ContactSolver> makeOpenMpContactSolver(const core::SimConfig& config);
#endif
#ifdef GRANIMPACT_ENABLE_CUDA
// Defined in src/backends/CudaSolver.cu (stage 5).
std::unique_ptr<ContactSolver> makeCudaContactSolver(const core::SimConfig& config);
#endif

 // State of each backend in this build: `granimpact_serial --backends`.
[[nodiscard]] std::vector<BackendStatus> backendStatus();

 // Factory. Throws std::invalid_argument with an actionable message if the backend
 // does not exist in this build (it never pretends to have a stage it does not have).
[[nodiscard]] std::unique_ptr<ContactSolver> makeContactSolver(const std::string& name,
                                                             const core::SimConfig& config);

}  // namespace granimpact::backends
