set(SOURCES
  adjoint.cpp
  diff_cache.cpp
  objective.cpp
  utils.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(polyfempy PRIVATE ${SOURCES})
