# Polyfem
# License: MIT

if(TARGET polyfem::polyfem)
    return()
endif()

get_filename_component(POLYFEMPY_REPO_ROOT "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
set(POLYFEMPY_POLYFEM_ROOT "${POLYFEMPY_REPO_ROOT}/polyfem")

message(STATUS "Third-party: using 'polyfem' submodule")

if(NOT EXISTS "${POLYFEMPY_POLYFEM_ROOT}/CMakeLists.txt")
    message(
        FATAL_ERROR
        "polyfem submodule is missing at: ${POLYFEMPY_POLYFEM_ROOT}\n"
        "Run: git submodule update --init --recursive\n"
        "Or clone with: git clone --recurse-submodules <repo-url>"
    )
endif()

option(POLYFEM_WITH_PYTHON "Enable Python expressions in PolyFEM" ON)

function(polyfempy_add_backend)
    # The parent supplies extension linking instead of embedding an interpreter.
    # Keep this override local so the user's cached option remains unchanged.
    set(POLYFEM_WITH_PYTHON OFF)
    add_subdirectory("${POLYFEMPY_POLYFEM_ROOT}" "${CMAKE_BINARY_DIR}/polyfem")
endfunction()

polyfempy_add_backend()

if(POLYFEM_WITH_PYTHON)
    find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)
    include("${POLYFEMPY_POLYFEM_ROOT}/cmake/recipes/pybind11.cmake")
    target_link_libraries(polyfem PUBLIC pybind11::module)
    target_compile_definitions(polyfem PUBLIC
        POLYFEM_WITH_PYTHON
        POLYFEM_PYTHON_EXECUTABLE="${Python3_EXECUTABLE}")
endif()
