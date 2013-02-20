# set up project
project("dune-grid" C CXX)

# general stuff
cmake_minimum_required(VERSION 2.8.8)
message("CMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}")
# make sure our own modules are found
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

add_subdirectory("m4")
add_subdirectory("dune")
add_subdirectory("doc")
add_subdirectory("lib")

export(TARGETS dunegrid ${DUNE_ALBERTA_LIBS} FILE dune-grid-targets.cmake)
# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
