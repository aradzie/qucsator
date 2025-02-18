set(
  TEMPLATE_SOURCES
  eqnsys.h
  hash.h
  ptrlist.h
  states.h
  tmatrix.h
  tridiag.h
  tvector.h
  tvector.h
  valuelist.h
)

set(
  SOURCES
  ${TEMPLATE_SOURCES}
  circuit.cpp
  dataset.cpp
  devstates.cpp
  differentiate.cpp
  environment.cpp
  equation.cpp
  evaluate.cpp
  exception.cpp
  exceptionstack.cpp
  fourier.cpp
  history.cpp
  input.cpp
  integrator.cpp
  interpolator.cpp
  logging.c
  matvec.cpp
  module.cpp
  net.cpp
  nodelist.cpp
  nodeset.cpp
  object.cpp
  parasweep.cpp
  property.cpp
  range.cpp
  receiver.cpp
  spline.cpp
  strlist.cpp
  sweep.cpp
  transient.cpp
  variable.cpp
  vector.cpp
)

include_directories(
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/analyses
  ${CMAKE_SOURCE_DIR}/src/components
  ${CMAKE_SOURCE_DIR}/src/math
  ${CMAKE_SOURCE_DIR}/src/parsers
)

add_subdirectory(analyses)
add_subdirectory(components)
add_subdirectory(components/devices)
add_subdirectory(components/digital)
add_subdirectory(components/microstrip)
add_subdirectory(math)
add_subdirectory(parsers)

add_executable(qucsator main.cpp ${SOURCES})

target_link_libraries(
  qucsator PRIVATE
  coreAnalyses
  coreComponents
  coreDevices
  coreDigital
  coreMath
  coreMicrostrip
  coreParsers
)
