# nanobind (https://github.com/pywrap/nanobind)
# License: BSD-style
if(TARGET nanobind::nanobind)
    return()
endif()

message(STATUS "Third-party: creating target 'nanobind::nanobind'")

if (POLICY CMP0094)  # https://cmake.org/cmake/help/latest/policy/CMP0094.html
    cmake_policy(SET CMP0094 NEW)  # FindPython should return the first matching Python
endif ()

# needed on GitHub Actions CI: actions/setup-python does not touch registry/frameworks on Windows/macOS
# this mirrors PythonInterp behavior which did not consult registry/frameworks first
if (NOT DEFINED Python_FIND_REGISTRY)
    set(Python_FIND_REGISTRY "LAST")
endif ()
if (NOT DEFINED Python_FIND_FRAMEWORK)
    set(Python_FIND_FRAMEWORK "LAST")
endif ()

# Find Python for nanobind
if(DEFINED PYTHON_EXECUTABLE AND NOT DEFINED Python_EXECUTABLE)
    set(Python_EXECUTABLE "${PYTHON_EXECUTABLE}" CACHE FILEPATH "Python executable for FindPython" FORCE)
endif()
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
set(PYTHON_EXECUTABLE ${Python_EXECUTABLE})

# Try to find nanobind from installed package first
execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE nanobind_DIR
    RESULT_VARIABLE nanobind_found_result
    ERROR_QUIET
)

if(nanobind_found_result EQUAL 0 AND EXISTS "${nanobind_DIR}")
    message(STATUS "Found nanobind via pip: ${nanobind_DIR}")

    # Check if nanobind has bundled submodule (ext/robin-map)
    # pip-installed nanobind may or may not include submodules
    set(NANOBIND_EXT_DIR "${nanobind_DIR}/../ext")
    set(NANOBIND_ROBIN_MAP_DIR "${NANOBIND_EXT_DIR}/robin-map")
    
    if(EXISTS "${NANOBIND_ROBIN_MAP_DIR}/include/tsl/robin_map.h")
        message(STATUS "Found nanobind bundled robin-map submodule at: ${NANOBIND_ROBIN_MAP_DIR}")
        # Use bundled submodule - this ensures correct 1.x version
        set(NB_USE_SUBMODULE_DEPS ON CACHE BOOL "Use nanobind bundled submodule dependencies" FORCE)
    else()
        message(STATUS "nanobind bundled robin-map not found, will use CPM or existing declaration")
        # No bundled submodule - we need to provide compatible robin-map 1.x
        # Check if robin-map is already declared (e.g., by top-level CMakeLists.txt)
        set(NB_USE_SUBMODULE_DEPS OFF CACHE BOOL "Use external tsl::robin_map" FORCE)
        
        # Check if tsl::robin_map target already exists (from top-level CMakeLists.txt)
        if(NOT TARGET tsl::robin_map)
            message(STATUS "robin-map not pre-declared, fetching compatible version via CPM")
            # Fetch compatible robin-map 1.x via CPM before finding nanobind
            # This ensures find_dependency(tsl-robin-map) finds the correct version
            # nanobind requires: >= 1.3.0, < 2.0.0
            include(CPM)
            CPMAddPackage(
                NAME robin-map
                GITHUB_REPOSITORY Tessil/robin-map
                GIT_TAG v1.4.1  # Use 1.4.1 (recommended) or at least >= 1.3.0, < 2.0.0
                OPTIONS
                    "ROBIN_MAP_BUILD_TESTS OFF"
                    "ROBIN_MAP_BUILD_EXAMPLES OFF"
            )
            
            if(TARGET tsl::robin_map)
                message(STATUS "Successfully fetched robin-map 1.x via CPM")
            else()
                message(FATAL_ERROR "Failed to fetch compatible robin-map 1.x via CPM")
            endif()
        else()
            message(STATUS "Using pre-declared robin-map (should be >= 1.3.0, < 2.0.0)")
        endif()
    endif()
    
    # Explicitly disable static library mode (use shared library)
    set(NB_STATIC_LIB OFF CACHE BOOL "Disable nanobind static library" FORCE)
    
    # Make sure CMake searches this directory for nanobindConfig
    list(PREPEND CMAKE_PREFIX_PATH "${nanobind_DIR}")

    # Find nanobind package
    find_package(nanobind CONFIG REQUIRED)
    
    # Diagnostic: Verify robin_map version after find_package
    # Check if any robin_map header is found and verify version
    # Note: We don't use NO_DEFAULT_PATH to allow searching CMAKE_PREFIX_PATH (including conda paths)
    find_path(ROBIN_MAP_HEADER_PATH
        NAMES tsl/robin_map.h
        PATHS
            "${NANOBIND_ROBIN_MAP_DIR}/include"
        PATH_SUFFIXES
            include
    )
    
    if(ROBIN_MAP_HEADER_PATH)
        message(STATUS "Using robin_map header at: ${ROBIN_MAP_HEADER_PATH}")
        # Try to detect version by reading the header
        file(READ "${ROBIN_MAP_HEADER_PATH}/tsl/robin_map.h" ROBIN_MAP_HEADER_CONTENT LIMIT 1000)
        if(ROBIN_MAP_HEADER_CONTENT MATCHES "TSL_ROBIN_MAP_VERSION_MAJOR[ \t]+([0-9]+)")
            set(ROBIN_MAP_VERSION_MAJOR ${CMAKE_MATCH_1})
            message(STATUS "Detected robin_map major version: ${ROBIN_MAP_VERSION_MAJOR}")
            if(NOT ROBIN_MAP_VERSION_MAJOR EQUAL 1)
                message(FATAL_ERROR 
                    "Incompatible tsl::robin_map version detected: ${ROBIN_MAP_VERSION_MAJOR}.x\n"
                    "nanobind requires version >= 1.3.0, < 2.0.0.\n"
                    "Found header at: ${ROBIN_MAP_HEADER_PATH}\n"
                    "This will cause compilation to fail. Please ensure a compatible 1.x version is used."
                )
            endif()
        else()
            # Try alternative version detection
            if(ROBIN_MAP_HEADER_CONTENT MATCHES "version[ \t]+([0-9]+)\\.([0-9]+)")
                set(ROBIN_MAP_VERSION_MAJOR ${CMAKE_MATCH_1})
                if(NOT ROBIN_MAP_VERSION_MAJOR EQUAL 1)
                    message(FATAL_ERROR "Incompatible robin_map version: ${CMAKE_MATCH_1}.${CMAKE_MATCH_2}")
                endif()
            else()
                message(WARNING "Could not detect robin_map version from header, proceeding with caution")
            endif()
        endif()
    else()
        message(WARNING "Could not locate robin_map header for version check - may cause build failure")
    endif()
    
    # After find_package, ensure static library is not enabled
    if(TARGET nanobind-static)
        message(WARNING "nanobind-static target found, but should use shared library mode")
    endif()

else()
    message(FATAL_ERROR 
        "nanobind not found via pip. Please install it with:\n"
        "  python -m pip install nanobind\n"
        "Static nanobind build is not supported (requires tsl::robin_map dependency)."
    )
endif()

# Verify nanobind_add_module is available
if(COMMAND nanobind_add_module)
    message(STATUS "nanobind_add_module() function is available")
else()
    message(WARNING "nanobind_add_module() function not found, may need manual setup")
endif()
