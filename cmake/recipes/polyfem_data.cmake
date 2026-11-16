# data
# License: MIT

message(STATUS "Third-party: using 'polyfem data' submodule")

if(NOT EXISTS "${POLYFEMPY_DATA_ROOT}/README.md")
    message(
        FATAL_ERROR
        "polyfem-data submodule is missing at: ${POLYFEMPY_DATA_ROOT}\n"
        "Run: git submodule update --init --recursive\n"
        "Or configure with: -DINPUT_POLYFEMPY_DATA_ROOT=path/to/polyfem-data"
    )
endif()

message(STATUS "Using polyfem data directory: ${POLYFEMPY_DATA_ROOT}")
set(polyfem_data_POPULATED TRUE)
