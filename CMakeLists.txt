# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

cmake_minimum_required(VERSION 3.13)
project("dune-grid" C CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
        "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
  "${PROJECT_SOURCE_DIR}/cmake/modules")

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

# create the dune library for dune-grid
dune_add_library(dunegrid)
target_link_libraries(dunegrid PUBLIC ${DUNE_LIBS})

# create a dune library for the albertagrid
if(Alberta_FOUND)
  foreach(_dim ${ALBERTA_WORLD_DIMS})
    dune_add_library(dunealbertagrid${_dim}d NO_MODULE_LIBRARY)
    target_link_libraries(dunealbertagrid${_dim}d PUBLIC
      dunegrid Alberta::AlbertaGrid_${_dim}d)
    target_compile_definitions(dunealbertagrid${_dim}d PUBLIC
      ALBERTA_DIM=${_dim} ENABLE_ALBERTA=1)
  endforeach(_dim)
endif()

# create alias HAVE_UG
set(HAVE_UG ${HAVE_DUNE_UGGRID})

add_subdirectory("dune")
add_subdirectory("doc")
add_subdirectory("src")
add_subdirectory("cmake/modules")

# if Python bindings are enabled, include necessary sub directories.
if(DUNE_ENABLE_PYTHONBINDINGS)
  add_subdirectory("python")
endif()

# Use a  custom section to export path to the grid examples.
set(DUNE_CUSTOM_PKG_CONFIG_SECTION "#Export the directory with the grid example for downstream modules
if(${ProjectName}_INSTALLED)
  set(DUNE_GRID_EXAMPLE_GRIDS_PATH \"\${PACKAGE_PREFIX_DIR}/${CMAKE_INSTALL_DOCDIR}/grids/\")
else()
  set(DUNE_GRID_EXAMPLE_GRIDS_PATH \"\${PACKAGE_PREFIX_DIR}/doc/grids/\")
endif()
")
# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
