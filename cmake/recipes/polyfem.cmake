# Polyfem
# License: MIT

if(TARGET polyfem::polyfem)
    return()
endif()

message(STATUS "Third-party: creating target 'polyfem::polyfem'")

if (POLICY CMP0079)  # https://cmake.org/cmake/help/latest/policy/CMP0079.html
    cmake_policy(SET CMP0079 NEW)  # FindPython should return the first matching Python
endif ()


include(FetchContent)
FetchContent_Declare(
    polyfem
    GIT_REPOSITORY https://github.com/polyfem/polyfem.git
    GIT_TAG c0b5cd24021b26690095779707e2b0e1656073f5
    GIT_SHALLOW FALSE
)
FetchContent_MakeAvailable(polyfem)


