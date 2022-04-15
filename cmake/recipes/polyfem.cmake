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
    GIT_TAG ffb1d130e70fc913ec4546c8683d6d22bd58a700
    GIT_SHALLOW FALSE
)
FetchContent_MakeAvailable(polyfem)


