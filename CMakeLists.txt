################################################################################
# Check required CMake version
set(REQUIRED_CMAKE_VERSION "3.25.0")
cmake_minimum_required(VERSION ${REQUIRED_CMAKE_VERSION})

if(INPUT_POLYFEMPY_DATA_ROOT)
    set(POLYFEMPY_DATA_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/${INPUT_POLYFEMPY_DATA_ROOT}/")
else()
    set(POLYFEMPY_DATA_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/data/")
endif()
if(NOT EXISTS ${POLYFEMPY_DATA_ROOT})
    file(MAKE_DIRECTORY ${POLYFEMPY_DATA_ROOT})
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
pybind11_add_module(polyfempy)
add_subdirectory(src)
target_include_directories(polyfempy PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(polyfempy PRIVATE polyfem::polyfem pybind11::module pybind11::json)

set_target_properties(polyfempy PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}" SUFFIX "${PYTHON_MODULE_EXTENSION}")
