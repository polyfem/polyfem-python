################################################################################
# Check required CMake version
set(REQUIRED_CMAKE_VERSION "3.25.0")
cmake_minimum_required(VERSION ${REQUIRED_CMAKE_VERSION})

if(INPUT_POLYFEMPY_DATA_ROOT)
    if(IS_ABSOLUTE "${INPUT_POLYFEMPY_DATA_ROOT}")
        set(POLYFEMPY_DATA_ROOT "${INPUT_POLYFEMPY_DATA_ROOT}")
    else()
        set(POLYFEMPY_DATA_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/${INPUT_POLYFEMPY_DATA_ROOT}")
    endif()
else()
    set(POLYFEMPY_DATA_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/polyfem-data")
endif()

project(polyfempy)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

################################################################################

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/recipes/")

include(polyfem_cpm_cache)

# Color output
include(UseColors)

# Prepend function
include(PrependCurrentPath)

# Extra warnings
include(Warnings)

# Enable compiler features required by dependencies.
include(CXXFeatures)

# Sort projects inside the solution
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Generate position independent code by default
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################
# Polyfem library
################################################################################

# Use nanobind for Python bindings (pybind11 support removed)
include(CPM)
# Force robin-map to 1.4.1 (or >=1.3.0) BEFORE including polyfem
# This ensures CPM uses the correct version, overriding any indirect dependencies
# (e.g., polyfem's ipc-toolkit may try to use robin-map 1.2.1, which is too old)
# nanobind requires: >= 1.3.0, < 2.0.0
    CPMAddPackage(
        NAME robin-map
        GITHUB_REPOSITORY Tessil/robin-map
        GIT_TAG v1.4.1  # Use 1.4.1 (recommended) or at least >= 1.3.0, < 2.0.0
        OPTIONS
            "ROBIN_MAP_BUILD_TESTS OFF"
            "ROBIN_MAP_BUILD_EXAMPLES OFF"
    )
    if(TARGET tsl::robin_map)
        message(STATUS "Pre-declared robin-map 1.4.1 for nanobind compatibility")
    else()
        message(WARNING "Failed to pre-declare robin-map 1.4.1 - may cause version conflicts")
endif()

# dependencies
set(POLYFEM_PYTHON_MODULE ON CACHE BOOL "Use Python support inside a Python extension" FORCE)
set(POLYFEM_WITH_TESTS OFF CACHE BOOL "Build PolyFEM tests")
include(polyfem)

# This executable embeds schemas at build time. Only runtime binaries belong in the wheel.
if(TARGET jse_embed_spec_tool)
    set_target_properties(jse_embed_spec_tool PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tools/$<CONFIG>"
        RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_CURRENT_BINARY_DIR}/tools/Release")
endif()

    # Use a shared nanobind runtime alongside the extension.
    set(NB_STATIC_LIB OFF CACHE BOOL "Disable nanobind static library" FORCE)
    
    message(STATUS "Using nanobind for Python bindings")
    include(nanobind)

# Mesh data is only needed for examples and numerical tests.
if(INPUT_POLYFEMPY_DATA_ROOT)
    include(polyfem_data)
endif()

################################################################################
# Subdirectories
################################################################################

    # Place the shared runtime and the extension in the same package.
    nanobind_add_module(polyfempy NB_SHARED)
    
    # Verify no static library target was created (safety check)
    if(TARGET nanobind-static)
        message(FATAL_ERROR 
            "Unexpected static nanobind target.\n"
            "Keep NB_STATIC_LIB disabled and use NB_SHARED for this extension."
        )
    endif()
    
    add_subdirectory(src)
    target_include_directories(polyfempy PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_link_libraries(polyfempy PRIVATE polyfem::polyfem)
    # Define macro for source files
    target_compile_definitions(polyfempy PRIVATE USE_NANOBIND)

# Ensure the produced module has the correct Python extension suffix.
# Prefer FindPython's canonical variable `Python_EXTENSION_MODULE_SUFFIX`.
# Example on Windows: `.cp311-win_amd64.pyd`
if(DEFINED Python_EXTENSION_MODULE_SUFFIX AND NOT "${Python_EXTENSION_MODULE_SUFFIX}" STREQUAL "")
    set_target_properties(polyfempy PROPERTIES PREFIX "" SUFFIX "${Python_EXTENSION_MODULE_SUFFIX}")
elseif(DEFINED Python_EXTENSION_SUFFIX AND NOT "${Python_EXTENSION_SUFFIX}" STREQUAL "")
    # Backwards/alt name in some CMake setups.
    set_target_properties(polyfempy PROPERTIES PREFIX "" SUFFIX "${Python_EXTENSION_SUFFIX}")
elseif(WIN32)
    # Last resort: ensure it is importable as a Python extension on Windows.
    set_target_properties(polyfempy PROPERTIES PREFIX "" SUFFIX ".pyd")
else()
    message(WARNING "Python extension module suffix not defined; leaving default target suffix.")
endif()
