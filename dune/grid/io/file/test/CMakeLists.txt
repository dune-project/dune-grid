# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

# TODO when porting this file to dune_add_test, the GMSH test needed to be adapted.
# It previously relied on the flags added to the executable to switch different grid types.
# I have adapted it quickly by adding a preprocesor variable. This could be and should
# be done more nicely.
#

dune_add_test(SOURCES conformvolumevtktest.cc)

dune_add_test(SOURCES gnuplottest.cc)

dune_add_test(SOURCES nonconformboundaryvtktest.cc)

dune_add_test(SOURCES printgridtest.cc)

dune_add_test(NAME subsamplingvtktest
              SOURCES subsamplingvtktest.cc test-linking.cc
              TIMEOUT 600)

dune_add_test(SOURCES vtktest.cc
              MPI_RANKS 1 2
              TIMEOUT 1200)

dune_add_test(SOURCES vtksequencetest.cc)

dune_add_test(SOURCES starcdreadertest.cc
              LINK_LIBRARIES dunegrid
              COMPILE_DEFINITIONS DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\"
              CMAKE_GUARD dune-uggrid_FOUND)

# the gmsh tests
dune_add_test(NAME gmshtest-onedgrid
              SOURCES gmshtest.cc
              LINK_LIBRARIES dunegrid
              COMPILE_DEFINITIONS GMSH_ONEDGRID
                                  DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")

dune_add_test(NAME gmshtest-uggrid
              SOURCES gmshtest.cc
              LINK_LIBRARIES dunegrid
              COMPILE_DEFINITIONS GMSH_UGGRID
                                  DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\"
              CMAKE_GUARD dune-uggrid_FOUND)

if(Alberta_FOUND)
  add_executable(gmshtest-alberta2d gmshtest.cc)
  target_link_libraries(gmshtest-alberta2d PRIVATE dunegrid)
  add_dune_alberta_flags(gmshtest-alberta2d WORLDDIM 2)
  target_compile_definitions(gmshtest-alberta2d PUBLIC GMSH_ALBERTAGRID DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
  dune_add_test(TARGET gmshtest-alberta2d)

  add_executable(gmshtest-alberta3d gmshtest.cc)
  target_link_libraries(gmshtest-alberta3d PRIVATE dunegrid)
  add_dune_alberta_flags(gmshtest-alberta3d WORLDDIM 3)
  target_compile_definitions(gmshtest-alberta3d PUBLIC GMSH_ALBERTAGRID DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
  dune_add_test(TARGET gmshtest-alberta3d)
endif()
