include_directories(
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/math
)

set(
  SOURCES
  complex.cpp
  fspecial.cpp
  matrix.cpp
  real.cpp
)

add_library(coreMath OBJECT ${SOURCES})
