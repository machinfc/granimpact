 // Reference (serial) implementation and contact solver factory.
//
 // Adding an accelerated stage means adding a file right here:
 // src/backends/Backend.cpp (this factory)
//   src/backends/OpenMpBackend.cpp     (etapa 4)
//   src/backends/CudaBackend.cu        (etapa 5)
 // Nothing else has to be touched: not the core, not the analysis, not the outputs.
#include "SerialSolver.hpp"

#include <stdexcept>

namespace granimpact::backends {

std::vector<BackendStatus> backendStatus() {
    std::vector<BackendStatus> list;
    list.push_back({"serial", "single-thread CPU, correctness reference", true});
#ifdef GRANIMPACT_ENABLE_OPENMP
    list.push_back({"openmp", "CPU multihilo (OpenMP)", false});
#else
    list.push_back({"openmp", "multithreaded CPU; build with -DGRANIMPACT_ENABLE_OPENMP=ON", false});
#endif
#ifdef GRANIMPACT_ENABLE_CUDA
    list.push_back({"cuda", "GPU (CUDA); kernels en src/backends/CudaBackend.cu", false});
#else
    list.push_back({"cuda", "GPU; build with -DGRANIMPACT_ENABLE_CUDA=ON and nvcc", false});
#endif
    return list;
}

std::unique_ptr<ContactSolver> makeContactSolver(const std::string& name,
                                                const core::SimConfig& config) {
    if (name.empty() || name == "serial") {
        return std::make_unique<CpuContactSolver>(config);
    }
#ifdef GRANIMPACT_ENABLE_OPENMP
    if (name == "openmp") {
        return makeOpenMpContactSolver(config);   // src/backends/OpenMpSolver.cpp
    }
#endif
#ifdef GRANIMPACT_ENABLE_CUDA
    if (name == "cuda") {
        return makeCudaContactSolver(config);     // src/backends/CudaSolver.cu
    }
#endif
    // Accelerated backends are registered here once they exist (stages 4 and 5).
    // Until then, the message says exactly what is missing: it never fakes their existence.
    if (name == "openmp" || name == "cuda") {
        throw std::invalid_argument(
            "backend '" + name + "' is not implemented in this build yet. "
            "Use --backend serial. See docs/architecture.md (stages 4-5).");
    }
    throw std::invalid_argument("unknown backend: '" + name + "'");
}

}  // namespace granimpact::backends
