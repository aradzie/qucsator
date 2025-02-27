include_directories(
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/analyses
  ${CMAKE_SOURCE_DIR}/src/math
)

set(
  SOURCES
  acsolver.cpp
  analysis.cpp
  dcsolver.cpp
  hbsolver.cpp
  nasolver.h
  parasweep.cpp
  spsolver.cpp
  trsolver.cpp
  trsolver_dc.cpp
)

add_library(coreAnalyses OBJECT ${SOURCES})
