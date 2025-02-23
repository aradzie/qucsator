include_directories(
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/components
  ${CMAKE_SOURCE_DIR}/src/components/helpers
  ${CMAKE_SOURCE_DIR}/src/math
)

set(
  SOURCES
  cross.cpp
  ground.cpp
  itrafo.cpp
  open.cpp
  short.cpp
  tee.cpp
)

add_library(coreCompHelpers OBJECT ${SOURCES})
