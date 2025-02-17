include_directories(
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/components
  ${CMAKE_SOURCE_DIR}/src/components/devices
  ${CMAKE_SOURCE_DIR}/src/components/microstrip
  ${CMAKE_SOURCE_DIR}/src/math
)

set(
  SOURCES
  bondwire.cpp
  circularloop.cpp
  cpwgap.cpp
  cpwline.cpp
  cpwopen.cpp
  cpwshort.cpp
  cpwstep.cpp
  mscorner.cpp
  mscoupled.cpp
  mscross.cpp
  msgap.cpp
  mslange.cpp
  msline.cpp
  msmbend.cpp
  msopen.cpp
  msrstub.cpp
  msstep.cpp
  mstee.cpp
  msvia.cpp
  spiralinductor.cpp
  substrate.cpp
)

add_library(coreMicrostrip OBJECT ${SOURCES})
