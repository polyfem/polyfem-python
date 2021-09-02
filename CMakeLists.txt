################################################################################
cmake_minimum_required(VERSION 3.1)
project(polyfempy)
################################################################################

if(INPUT_THIRD_PARTY_DIR)
	set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${INPUT_THIRD_PARTY_DIR}/)
else()
	set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/)
endif()
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

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


# Setup dependencies
include(PolyfemPythonDependencies)

################################################################################
# Polyfem library
################################################################################

set(BUILD_TESTS		OFF CACHE BOOL "" FORCE)
set(DOWNLOAD_GTEST	OFF CACHE BOOL "" FORCE)



# polyfem
polyfem_python_download_polyfem()
add_subdirectory(${THIRD_PARTY_DIR}/polyfem)
add_library(nlohmann_json ALIAS json)

# pybind11
polyfem_python_download_pybind11()
add_subdirectory(${THIRD_PARTY_DIR}/pybind11)

# pybind11_json
polyfem_python_download_pybind11_json()
add_subdirectory(${THIRD_PARTY_DIR}/pybind11_json)




#for testing purpose
polyfem_python_download_data()


################################################################################
# Subdirectories
################################################################################
add_library(polyfempy MODULE src/binding.cpp src/raster.cpp)
target_link_libraries(polyfempy PRIVATE polyfem pybind11_json)
target_link_libraries(polyfempy PRIVATE pybind11::module)
set_target_properties(polyfempy PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}" SUFFIX "${PYTHON_MODULE_EXTENSION}")




