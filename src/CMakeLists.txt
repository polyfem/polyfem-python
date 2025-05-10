set(SOURCES
  binding.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(polyfempy PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

add_subdirectory(differentiable)
add_subdirectory(mesh)
add_subdirectory(state)
add_subdirectory(solver)
