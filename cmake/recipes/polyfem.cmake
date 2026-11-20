include_guard(GLOBAL)
include(CPM)
string(COMPARE EQUAL "${SKBUILD_STATE}" "sdist" download_only)
CPMAddPackage(
    NAME polyfem
    GITHUB_REPOSITORY polyfem/polyfem
    GIT_TAG 76205bfd1247570abbb9e552347b8b44cc62b0b2
    DOWNLOAD_ONLY "${download_only}")
