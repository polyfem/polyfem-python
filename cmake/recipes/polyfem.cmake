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
    GIT_TAG 76afcb9bd6d2dcd24a344478f1d7025540b05395
    GIT_SHALLOW FALSE
)
FetchContent_MakeAvailable(polyfem)


