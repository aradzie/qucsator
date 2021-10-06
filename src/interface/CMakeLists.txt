include_directories(
  ${qucs-core_SOURCE_DIR}
  ${qucs-core_CURRENT_SOURCE_DIR}
  ${qucs-core_SOURCE_DIR}/src/math
  ${qucs-core_SOURCE_DIR}/src/components
  ${qucs-core_BINARY_DIR}/src/components)

set(INTERFACE_SRC qucs_interface.cpp e_trsolver.cpp)

set(HEADERS qucs_interface.h)

add_library(coreInterface OBJECT ${INTERFACE_SRC})

add_dependencies(coreInterface coreVerilog)

install(FILES ${HEADERS} DESTINATION include/qucs-core)
