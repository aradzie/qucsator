include_directories(
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/components
  ${CMAKE_SOURCE_DIR}/src/components/linear
  ${CMAKE_SOURCE_DIR}/src/math
)

set(
  SOURCES
  amplifier.cpp
  attenuator.cpp
  biastee.cpp
  capacitor.cpp
  capq.cpp
  coupler.cpp
  dcblock.cpp
  dcfeed.cpp
  gyrator.cpp
  hybrid.cpp
  indq.cpp
  inductor.cpp
  isolator.cpp
  mutual.cpp
  mutual2.cpp
  mutualx.cpp
  phaseshifter.cpp
  relais.cpp
  resistor.cpp
  tswitch.cpp
)

add_library(coreCompLinear OBJECT ${SOURCES})
