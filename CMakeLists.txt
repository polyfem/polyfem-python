cmake_minimum_required(VERSION 3.25)
project(polyfempy LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/recipes")

find_package(Python 3.10 COMPONENTS Interpreter Development.Module REQUIRED)
if(NOT SKBUILD_STATE STREQUAL "sdist")
    include(nanobind)
endif()
include(polyfem)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/python-from-jse/tools/generate_with_overrides.py")
    execute_process(
        COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tools/generate_polyfem_api.py"
                --polyfem-source-dir "${polyfem_SOURCE_DIR}"
        COMMAND_ERROR_IS_FATAL ANY)
elseif(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/python/polyfempy/generated_api/generated_api.py")
    message(FATAL_ERROR "Initialize the generator submodule: git submodule update --init --recursive")
endif()
if(SKBUILD_STATE STREQUAL "sdist")
    return()
endif()

nanobind_add_module(polyfempy NB_STATIC)
add_subdirectory(src)
target_include_directories(polyfempy PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(polyfempy PRIVATE polyfem::polyfem)
install(TARGETS polyfempy LIBRARY DESTINATION polyfempy COMPONENT python
                        RUNTIME DESTINATION polyfempy COMPONENT python)
install(FILES
    python-from-jse/generator/__init__.py
    python-from-jse/generator/model_builder.py
    python-from-jse/generator/id_relationships.py
    python-from-jse/generator/selection_refs.py
    DESTINATION polyfempy/generated_api/generator COMPONENT python)
install(FILES "${polyfem_SOURCE_DIR}/LICENSE" DESTINATION polyfempy/licenses
    RENAME PolyFEM.LICENSE COMPONENT python)
install(FILES "${nanobind_SOURCE_DIR}/LICENSE" DESTINATION polyfempy/licenses
    RENAME nanobind.LICENSE COMPONENT python)
