################################################################################
include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(POLYFEM_PYTHON_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(POLYFEM_PYTHON_EXTRA_OPTIONS "")
endif()

# Shortcut function
function(polyfem_python_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${POLYFEM_PYTHON_EXTERNAL}/${name}
        DOWNLOAD_DIR ${POLYFEM_PYTHON_EXTERNAL}/.cache/${name}
        QUIET
        ${POLYFEM_PYTHON_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################


function(polyfem_python_download_polyfem)
    polyfem_python_download_project(polyfem
        GIT_REPOSITORY https://github.com/polyfem/polyfem.git
        GIT_TAG        c62a62e83311689f05927f2196e47e6d59bd29c4
    )
endfunction()

function(polyfem_python_download_pybind11)
    polyfem_python_download_project(pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG        4c36fb7b1236fce25e00b63f357ccc36dc006662
    )
endfunction()

## data
function(polyfem_python_download_data)
    polyfem_python_download_project(data
        GIT_REPOSITORY https://github.com/polyfem/polyfem-data.git
        GIT_TAG        11767f7ce80f8932b5da4410d45ddb37b5be681d
    )
endfunction()