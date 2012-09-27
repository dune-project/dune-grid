# set up project
project("dune-grid" C CXX)

# general stuff
cmake_minimum_required(VERSION 2.8)

#find dune-common and set the module path
find_package(dune-common)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
  "${CMAKE_SOURCE_DIR}/cmake/modules")

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

add_subdirectory("m4")
add_subdirectory("dune")
add_subdirectory("doc")
add_subdirectory("lib")
add_subdirectory("cmake/modules")

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
