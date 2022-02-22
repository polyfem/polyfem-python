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
    GIT_TAG 709a609aa4e3221929ff1c83bb0a88dfcad112b6
    GIT_SHALLOW FALSE
)
FetchContent_MakeAvailable(polyfem)


