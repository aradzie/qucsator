include_directories(
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/components
  ${CMAKE_SOURCE_DIR}/src/math
)

set(
  SOURCES
  circline.cpp
  circulator.cpp
  coaxline.cpp
  ctline.cpp
  ecvs.cpp
  pac.cpp
  rectline.cpp
  rfedd.cpp
  rlcg.cpp
  spdeembed.cpp
  spembed.cpp
  spfile.cpp
  strafo.cpp
  taperedline.cpp
  tline.cpp
  tline4p.cpp
  trafo.cpp
  twistedpair.cpp
)

add_library(coreComponents OBJECT ${SOURCES})
