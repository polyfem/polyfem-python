# Polyfem
# License: MIT

if(TARGET polyfem::polyfem)
    return()
endif()

message(STATUS "Third-party: creating target 'polyfem::polyfem'")

include(FetchContent)
FetchContent_Declare(
    polyfem
    GIT_REPOSITORY https://github.com/polyfem/polyfem.git
    GIT_TAG fcab1f5339e1f25eb2f2d214d17b7ec68a1bc0a2
    GIT_SHALLOW FALSE
)
FetchContent_MakeAvailable(polyfem)


