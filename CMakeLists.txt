cmake_minimum_required(VERSION 3.30)

project(qucsator CXX C)

set(CMAKE_VERBOSE_MAKEFILE OFF)

set(CMAKE_CXX_STANDARD 17)

add_compile_definitions(DEBUG)

add_subdirectory(src)
