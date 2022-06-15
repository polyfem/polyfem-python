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
    GIT_TAG ea0e5cc25f2ffba64ea47225cf106460935cda5f
    GIT_SHALLOW FALSE
)
FetchContent_MakeAvailable(polyfem)


