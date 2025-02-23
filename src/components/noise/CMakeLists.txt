include_directories(
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/components
  ${CMAKE_SOURCE_DIR}/src/components/noise
  ${CMAKE_SOURCE_DIR}/src/math
)

set(
  SOURCES
  iinoise.cpp
  inoise.cpp
  ivnoise.cpp
  ivnoise.cpp
  vnoise.cpp
  vvnoise.cpp
)

add_library(coreCompNoise OBJECT ${SOURCES})
