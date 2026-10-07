# One executable per test file: that is what allows running only what matters
# while iterating on the physics (ctest -R contact).
include(FetchContent)

if(GRANIMPACT_USE_SYSTEM_DEPS)
    find_package(GTest REQUIRED)
    set(GRANIMPACT_GTEST_TARGET GTest::gtest_main)
else()
    FetchContent_Declare(googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0
        GIT_SHALLOW TRUE)
    # Not installed: we only want the test target.
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
    set(GRANIMPACT_GTEST_TARGET gtest_main)
endif()

function(granimpact_add_test name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE granimpact_core ${GRANIMPACT_GTEST_TARGET})
    target_compile_options(${name} PRIVATE ${GRANIMPACT_WARNINGS})
    add_test(NAME ${name} COMMAND ${name})
    set_tests_properties(${name} PROPERTIES LABELS "unit")
endfunction()

granimpact_add_test(test_integrator)
granimpact_add_test(test_contact_model)
granimpact_add_test(test_cell_list)
granimpact_add_test(test_generation)
granimpact_add_test(test_crater_analysis)
granimpact_add_test(test_simulation_smoke)

# ------------------------------------------------------------------- python ----
# Post-processing is tested too: these are the same CSVs the simulator produces.
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    add_test(NAME test_postprocess
             COMMAND ${Python3_EXECUTABLE} -m unittest discover -s python/tests -v
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..)
endif()
