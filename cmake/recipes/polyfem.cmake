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
    GIT_TAG 12ac634833f91a3946cff26db01972fdb2ec3214
    GIT_SHALLOW FALSE
)
FetchContent_MakeAvailable(polyfem)


