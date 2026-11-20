include_guard(GLOBAL)
include(CPM)
find_package(Git REQUIRED)
string(COMPARE EQUAL "${SKBUILD_STATE}" "sdist" download_only)
CPMAddPackage(
    NAME polyfem
    URL https://github.com/polyfem/polyfem/archive/836c923cba93772bfc72daeaa42c7e0fcd6ca3e2.tar.gz
    DOWNLOAD_ONLY "${download_only}"
    PATCH_COMMAND "${GIT_EXECUTABLE}" apply
        "${CMAKE_CURRENT_LIST_DIR}/../patches/polyfem.patch")
