# libigl (https://github.com/libigl/libigl)
# License: MPL

#
# Copyright 2020 Adobe. All rights reserved.
# This file is licensed to you under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License. You may obtain a copy
# of the License at http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under
# the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR REPRESENTATIONS
# OF ANY KIND, either express or implied. See the License for the specific language
# governing permissions and limitations under the License.
#

if(TARGET igl::core)
    return()
endif()

# Keep dependency patches in the submodule and build integration in the parent.
message(STATUS "Third-party: creating target 'igl::core'")

set(LIBIGL_PREDICATES ON CACHE BOOL "Use exact predicates" FORCE)

include(CMakeDependentOption)
cmake_dependent_option(LIBIGL_RESTRICTED_TRIANGLE "Build target igl_restricted::triangle" ON "POLYFEM_WITH_TRIANGLE" ON)

find_package(Git REQUIRED)

file(GLOB_RECURSE patches_for_libigl CONFIGURE_DEPENDS
    "${POLYFEMPY_POLYFEM_ROOT}/cmake/patches/igl_*.patch"
)

include(CPM)
CPMAddPackage(
    NAME libigl
    GITHUB_REPOSITORY "libigl/libigl"
    GIT_TAG "aeeea9b416a5b5474cb9d12c07ff1238a4c83ba1"
    PATCH_COMMAND "${GIT_EXECUTABLE}" apply --ignore-whitespace --whitespace=nowarn ${patches_for_libigl}
)
