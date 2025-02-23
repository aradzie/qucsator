include_directories(
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/components
  ${CMAKE_SOURCE_DIR}/src/components/probes
  ${CMAKE_SOURCE_DIR}/src/math
)

set(
  SOURCES
  iprobe.cpp
  vprobe.cpp
  wprobe.cpp
)

add_library(coreCompProbes OBJECT ${SOURCES})
