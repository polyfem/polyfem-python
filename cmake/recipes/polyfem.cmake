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

add_subdirectory("${POLYFEMPY_POLYFEM_ROOT}" "${CMAKE_BINARY_DIR}/polyfem")
