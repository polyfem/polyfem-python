# Prepare dependencies
#
# For each third-party library, if the appropriate target doesn't exist yet,
# download it via external project, and add_subdirectory to build it alongside
# this project.

### Configuration
set(POLYFEM_PYTHON_ROOT     "${CMAKE_CURRENT_LIST_DIR}/..")
set(POLYFEM_PYTHON_EXTERNAL "${POLYFEM_PYTHON_ROOT}/3rdparty")

# Download and update 3rdparty libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
include(PolyfemPythonDownloadExternal)

################################################################################
# Required libraries
################################################################################

