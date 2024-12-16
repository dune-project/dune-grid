# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

# Link all test targets in this directory against Dune::Grid
link_libraries(Dune::Grid)

dune_add_test(NAME test-dgf-yasp
              SOURCES test-dgf-yasp.cc
              COMPILE_DEFINITIONS DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\"
             )

dune_add_test(NAME test-dgf-yasp-offset
              SOURCES test-dgf-yasp.cc
              COMPILE_DEFINITIONS DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\"
                                  TESTCOORDINATES
              CMD_ARGS ${PROJECT_SOURCE_DIR}/doc/grids/dgf/test2d_offset.dgf)

dune_add_test(NAME test-dgf-oned
              SOURCES test-dgf-oned.cc
              COMPILE_DEFINITIONS DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\"
             )

if(Alberta_FOUND)
  add_executable(test-dgf-alberta test-dgf-alberta.cc)
  add_dune_alberta_flags(test-dgf-alberta WORLDDIM 2)
  target_compile_definitions(test-dgf-alberta PUBLIC
    DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
  dune_add_test(TARGET test-dgf-alberta)
endif()

dune_add_test(NAME test-dgf-ug
              SOURCES test-dgf-ug.cc
              COMPILE_DEFINITIONS DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\"
              CMAKE_GUARD dune-uggrid_FOUND)

dune_add_test(NAME test-dgf-gmsh-ug
              SOURCES test-dgf-ug.cc
              COMPILE_DEFINITIONS DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\" DGFTEST_USE_GMSH=1
              CMAKE_GUARD dune-uggrid_FOUND)
