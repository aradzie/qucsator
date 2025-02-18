find_package(FLEX REQUIRED)
find_package(BISON REQUIRED)

set(GENERATED_SOURCES)

foreach(type csv citi dataset mdl netlist touchstone zvr)
  BISON_TARGET("Parser${type}"
          "${CMAKE_SOURCE_DIR}/src/parsers/parse_${type}.ypp"
          "${CMAKE_SOURCE_DIR}/src/parsers/parse_${type}.cpp"
          DEFINES_FILE "${CMAKE_SOURCE_DIR}/src/parsers/parse_${type}.h")
  FLEX_TARGET("Lexer${type}"
          "${CMAKE_SOURCE_DIR}/src/parsers/scan_${type}.lpp"
          "${CMAKE_SOURCE_DIR}/src/parsers/scan_${type}.cpp")
  ADD_FLEX_BISON_DEPENDENCY("Lexer${type}" "Parser${type}")

  list(APPEND GENERATED_SOURCES "${CMAKE_SOURCE_DIR}/src/parsers/parse_${type}.cpp")
  list(APPEND GENERATED_SOURCES "${CMAKE_SOURCE_DIR}/src/parsers/scan_${type}.cpp")
endforeach()

set(SOURCES
  ${GENERATED_SOURCES}
  check_citi.cpp
  check_csv.cpp
  check_dataset.cpp
  check_mdl.cpp
  check_netlist.cpp
  check_touchstone.cpp
  check_zvr.cpp
)

include_directories(
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/parsers
)

add_library(coreParsers OBJECT ${SOURCES})
