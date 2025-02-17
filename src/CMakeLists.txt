find_package(FLEX REQUIRED)
find_package(BISON REQUIRED)

set(GENERATED_SOURCES)

foreach(type csv citi dataset mdl netlist touchstone zvr)
  BISON_TARGET("Parser${type}"
          "${CMAKE_SOURCE_DIR}/src/parse_${type}.ypp"
          "${CMAKE_SOURCE_DIR}/src/parse_${type}.cpp"
          DEFINES_FILE "${CMAKE_SOURCE_DIR}/src/parse_${type}.h")
  FLEX_TARGET("Lexer${type}"
          "${CMAKE_SOURCE_DIR}/src/scan_${type}.lpp"
          "${CMAKE_SOURCE_DIR}/src/scan_${type}.cpp")
  ADD_FLEX_BISON_DEPENDENCY("Lexer${type}" "Parser${type}")

  list(APPEND GENERATED_SOURCES "${CMAKE_SOURCE_DIR}/src/parse_${type}.cpp")
  list(APPEND GENERATED_SOURCES "${CMAKE_SOURCE_DIR}/src/scan_${type}.cpp")
endforeach()

set(
  TEMPLATE_SOURCES
  eqnsys.h
  hash.h
  nasolution.h
  nasolver.h
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
  ${GENERATED_SOURCES}
  ${TEMPLATE_SOURCES}
  acsolver.cpp
  analysis.cpp
  check_citi.cpp
  check_csv.cpp
  check_dataset.cpp
  check_mdl.cpp
  check_netlist.cpp
  check_touchstone.cpp
  check_zvr.cpp
  circuit.cpp
  dataset.cpp
  dcsolver.cpp
  devstates.cpp
  differentiate.cpp
  environment.cpp
  equation.cpp
  evaluate.cpp
  exception.cpp
  exceptionstack.cpp
  fourier.cpp
  hbsolver.cpp
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
  spsolver.cpp
  strlist.cpp
  sweep.cpp
  transient.cpp
  trsolver.cpp
  variable.cpp
  vector.cpp
)

include_directories(
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/components
  ${CMAKE_SOURCE_DIR}/src/math
)

add_subdirectory(components)
add_subdirectory(components/devices)
add_subdirectory(components/digital)
add_subdirectory(components/microstrip)
add_subdirectory(math)

add_executable(qucsator main.cpp ${SOURCES})

target_link_libraries(qucsator PRIVATE coreComponents coreDevices coreDigital coreMicrostrip coreMath)
