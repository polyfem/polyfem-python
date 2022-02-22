################################################################################
# Check required CMake version
set(REQUIRED_CMAKE_VERSION "3.14.0")
if(POLYFEM_TOPLEVEL_PROJECT)
    cmake_minimum_required(VERSION ${REQUIRED_CMAKE_VERSION})
else()
    # Don't use cmake_minimum_required here to avoid implicitly overriding parent policies
    if(${CMAKE_VERSION} VERSION_LESS ${REQUIRED_CMAKE_VERSION})
        message(FATAL_ERROR "CMake required version to build PolyFEM is ${REQUIRED_CMAKE_VERSION}")
    endif()
endif()

if(INPUT_POLYFEMPY_DATA_ROOT)
    set(POLYFEMPY_DATA_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/${INPUT_POLYFEMPY_DATA_ROOT}/")
else()
    set(POLYFEMPY_DATA_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/data/")
endif()
if(NOT EXISTS ${POLYFEMPY_DATA_ROOT})
    file(MAKE_DIRECTORY ${POLYFEMPY_DATA_ROOT})
endif()

project(polyfempy)
################################################################################

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/recipes/")

# Color output
include(UseColors)

# Prepend function
include(PrependCurrentPath)

# Extra warnings
include(Warnings)

# Use C++11/14
include(CXXFeatures)

# Sort projects inside the solution
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Generate position independent code by default
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################
# Polyfem library
################################################################################

# dependencies
include(polyfem)
include(pybind11)
include(pybind11_json)

#for testing purpose
include(polyfem_data)

################################################################################
# Subdirectories
################################################################################
add_library(polyfempy MODULE src/binding.cpp src/raster.cpp)
target_link_libraries(polyfempy PRIVATE polyfem::polyfem pybind11::module pybind11::json)

set_target_properties(polyfempy PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}" SUFFIX "${PYTHON_MODULE_EXTENSION}")




