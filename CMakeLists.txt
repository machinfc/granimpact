# ============================================================================
#  GranImpact · DEM simulator for granular impact
#  ---------------------------------------------------------------------------
#  Build philosophy: ONE library with all the physics (granimpact_core) and thin
#  executables on top. OpenMP and CUDA are enabled by options and add code
#  (src/backends/*) without touching the core -> "just keep adding code".
#
#  Presets available: cmake --preset debug | release | asan | omp | cuda
# ============================================================================
cmake_minimum_required(VERSION 3.20)

project(granimpact
    VERSION 0.1.0
    DESCRIPTION "Discrete Element Method simulator for granular impact craters"
    HOMEPAGE_URL "https://github.com/machinfc/granimpact"
    LANGUAGES CXX)

# ------------------------------------------------------------------ options --
option(GRANIMPACT_BUILD_TESTS      "Build and register the unit tests" ON)
option(GRANIMPACT_BUILD_TOOLS      "Build auxiliary tools"           ON)
option(GRANIMPACT_ENABLE_OPENMP    "Enable the OpenMP backend"                OFF)
option(GRANIMPACT_ENABLE_CUDA      "Enable the CUDA backend"                  OFF)
option(GRANIMPACT_ENABLE_SANITIZERS "Build with ASan + UBSan (development)"  OFF)
option(GRANIMPACT_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(GRANIMPACT_USE_SYSTEM_DEPS  "Use system nlohmann_json/googletest (no download)" OFF)

# ------------------------------------------------------------------ standard --
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)   # for clangd / clang-tidy
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ------------------------------------------------------------------ warnings --
set(GRANIMPACT_WARNINGS
    -Wall -Wextra -Wpedantic -Wshadow -Wnon-virtual-dtor -Wcast-align
    -Wunused -Woverloaded-virtual -Wconversion -Wsign-conversion -Wdouble-promotion)
if(GRANIMPACT_WARNINGS_AS_ERRORS)
    list(APPEND GRANIMPACT_WARNINGS -Werror)
endif()

# ---------------------------------------------------------- dependencies ------
include(FetchContent)
set(FETCHCONTENT_QUIET OFF)

if(GRANIMPACT_USE_SYSTEM_DEPS)
    find_package(nlohmann_json 3.10 REQUIRED)
else()
    FetchContent_Declare(nlohmann_json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.11.3
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(nlohmann_json)
endif()

# ------------------------------------------------------------- library --------
file(GLOB_RECURSE GRANIMPACT_SOURCES CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/core/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/physics/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/generation/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/analysis/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/io/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/sim/*.cpp)

if(GRANIMPACT_ENABLE_CUDA)
    enable_language(CUDA)
    file(GLOB_RECURSE GRANIMPACT_CUDA_SOURCES CONFIGURE_DEPENDS
        ${CMAKE_CURRENT_SOURCE_DIR}/src/backends/cuda/*.cu)
    set(GRANIMPACT_SOURCES ${GRANIMPACT_SOURCES} ${GRANIMPACT_CUDA_SOURCES})
endif()

# CPU backends: serial is always there; OpenMP/CUDA only if enabled.
file(GLOB_RECURSE GRANIMPACT_BACKEND_SOURCES CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/backends/*.cpp)
list(APPEND GRANIMPACT_SOURCES ${GRANIMPACT_BACKEND_SOURCES})
if(GRANIMPACT_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
endif()

add_library(granimpact_core STATIC ${GRANIMPACT_SOURCES})
add_library(granimpact::core ALIAS granimpact_core)

target_include_directories(granimpact_core
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
           $<INSTALL_INTERFACE:include>)

target_link_libraries(granimpact_core PUBLIC nlohmann_json::nlohmann_json)
target_compile_options(granimpact_core PRIVATE ${GRANIMPACT_WARNINGS})

if(GRANIMPACT_ENABLE_OPENMP)
    target_link_libraries(granimpact_core PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(granimpact_core PUBLIC GRANIMPACT_ENABLE_OPENMP=1)
endif()

if(GRANIMPACT_ENABLE_CUDA)
    target_compile_definitions(granimpact_core PUBLIC GRANIMPACT_ENABLE_CUDA=1)
    target_link_libraries(granimpact_core PUBLIC CUDA::cudart)
endif()

# --------------------------------------------------------------- executables --
add_executable(granimpact_serial src/apps/main_serial.cpp)
target_link_libraries(granimpact_serial PRIVATE granimpact_core)
target_compile_options(granimpact_serial PRIVATE ${GRANIMPACT_WARNINGS})

# The accelerated executables are added here once they exist (stages 4-5):
# if(GRANIMPACT_ENABLE_OPENMP) add_executable(granimpact_openmp src/apps/main_openmp.cpp) ... endif()

# ------------------------------------------------------------------ tests -----
if(GRANIMPACT_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# -------------------------------------------------------------- sanitizers ----
if(GRANIMPACT_ENABLE_SANITIZERS)
    add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer -g)
    add_link_options(-fsanitize=address,undefined)
endif()

# -------------------------------------------------------------------- install -
include(GNUInstallDirs)
install(TARGETS granimpact_core granimpact_serial
        EXPORT granimpactTargets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY configs/ DESTINATION ${CMAKE_INSTALL_DATADIR}/granimpact/configs)
install(DIRECTORY docs/ DESTINATION ${CMAKE_INSTALL_DATADIR}/granimpact/docs)
install(DIRECTORY python/granimpact_post
        DESTINATION ${CMAKE_INSTALL_DATADIR}/granimpact/python)
install(DIRECTORY scripts/ DESTINATION ${CMAKE_INSTALL_DATADIR}/granimpact/scripts
        USE_SOURCE_PERMISSIONS)

# `granimpact` launcher: picks the best compiled backend (cuda > openmp > serial).
# Once stages 4-5 exist, building them is enough; the launcher detects them.
install(PROGRAMS scripts/granimpact DESTINATION ${CMAKE_INSTALL_BINDIR})

message(STATUS "GranImpact ${PROJECT_VERSION}")
message(STATUS "  build type ....... ${CMAKE_BUILD_TYPE}")
message(STATUS "  tests ............ ${GRANIMPACT_BUILD_TESTS}")
message(STATUS "  backend OpenMP ... ${GRANIMPACT_ENABLE_OPENMP}")
message(STATUS "  backend CUDA ..... ${GRANIMPACT_ENABLE_CUDA}")
message(STATUS "  sanitizers ....... ${GRANIMPACT_ENABLE_SANITIZERS}")
