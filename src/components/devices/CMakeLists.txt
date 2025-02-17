include_directories(
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/components
  ${CMAKE_SOURCE_DIR}/src/components/devices
  ${CMAKE_SOURCE_DIR}/src/math
)

set(
  SOURCES
  bjt.cpp
  device.cpp
  diac.cpp
  diode.cpp
  eqndefined.cpp
  jfet.cpp
  mosfet.cpp
  thyristor.cpp
  triac.cpp
  tunneldiode.cpp
)

add_library(coreDevices OBJECT ${SOURCES})
