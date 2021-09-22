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
        GIT_TAG        75298361e234570cb8fee2df21438ccb6b826a4f
    )
endfunction()

function(polyfem_python_download_pybind11)
    polyfem_python_download_project(pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG        v2.7.1
    )
endfunction()

# BSD 3
function(polyfem_python_download_pybind11_json)
    polyfem_python_download_project(pybind11_json
        GIT_REPOSITORY https://github.com/pybind/pybind11_json.git
        GIT_TAG        0.2.11
    )
endfunction()

## data
function(polyfem_python_download_data)
    polyfem_python_download_project(data
        GIT_REPOSITORY https://github.com/polyfem/polyfem-data.git
        GIT_TAG        09ba9e7db4c1a87cdd16f53e0bf54e19d599a465
    )
endfunction()
