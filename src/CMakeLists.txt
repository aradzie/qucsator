set(
  TEMPLATE_SOURCES
  eqnsys.h
  hash.h
  ptrlist.h
  states.h
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
  logging.c
  matvec.cpp
  module.cpp
  net.cpp
  nodelist.cpp
  nodeset.cpp
  object.cpp
  property.cpp
  range.cpp
  receiver.cpp
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
add_subdirectory(components/digital)
add_subdirectory(components/helpers)
add_subdirectory(components/linear)
add_subdirectory(components/microstrip)
add_subdirectory(components/noise)
add_subdirectory(components/nonlinear)
add_subdirectory(components/probes)
add_subdirectory(components/sources)
add_subdirectory(math)
add_subdirectory(parsers)

add_executable(qucsator main.cpp ${SOURCES})

target_link_libraries(
  qucsator PRIVATE
  coreAnalyses
  coreCompDigtal
  coreCompHelpers
  coreCompLinear
  coreCompMicrostrip
  coreCompNoise
  coreCompNonlinear
  coreCompProbes
  coreCompSources
  coreComponents
  coreMath
  coreParsers
)

add_executable(example example.cpp ${SOURCES})

target_link_libraries(
  example PRIVATE
  coreAnalyses
  coreCompDigtal
  coreCompHelpers
  coreCompLinear
  coreCompMicrostrip
  coreCompNoise
  coreCompNonlinear
  coreCompProbes
  coreCompSources
  coreComponents
  coreMath
  coreParsers
)
