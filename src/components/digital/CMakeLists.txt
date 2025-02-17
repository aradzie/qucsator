include_directories(
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/components
  ${CMAKE_SOURCE_DIR}/src/components/digital
  ${CMAKE_SOURCE_DIR}/src/math
)

set(
  SOURCES
  and.cpp
  buffer.cpp
  digisource.cpp
  digital.cpp
  inverter.cpp
  nand.cpp
  nor.cpp
  or.cpp
  xnor.cpp
  xor.cpp
)

add_library(coreDigital OBJECT ${SOURCES})
