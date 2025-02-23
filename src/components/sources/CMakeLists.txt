include_directories(
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/components
  ${CMAKE_SOURCE_DIR}/src/components/sources
  ${CMAKE_SOURCE_DIR}/src/math
)

set(
  SOURCES
  cccs.cpp
  ccvs.cpp
  iac.cpp
  idc.cpp
  iexp.cpp
  ifile.cpp
  ipulse.cpp
  irect.cpp
  vac.cpp
  vam.cpp
  vccs.cpp
  vcvs.cpp
  vdc.cpp
  vexp.cpp
  vfile.cpp
  vpm.cpp
  vpulse.cpp
  vrect.cpp
)

add_library(coreCompSources OBJECT ${SOURCES})
