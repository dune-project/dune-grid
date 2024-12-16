# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

# TODO when porting this file to dune_add_test, the GMSH test needed to be adapted.
# It previously relied on the flags added to the executable to switch different grid types.
# I have adapted it quickly by adding a preprocessor variable. This could be and should
# be done more nicely.
#

# Link all test targets in this directory against Dune::Grid
link_libraries(Dune::Grid)

dune_add_test(SOURCES conformvolumevtktest.cc)

dune_add_test(SOURCES gnuplottest.cc)

dune_add_test(SOURCES nonconformboundaryvtktest.cc)

dune_add_test(SOURCES printgridtest.cc)

if(TARGET Python3::Python)
  dune_add_test(NAME subsamplingvtktest
                LINK_LIBRARIES ${DUNE_LIBS} Python3::Python
                SOURCES subsamplingvtktest.cc test-linking.cc
                TIMEOUT 600
                COMPILE_DEFINITIONS PYTHON_INTERPRETER=L"${DUNE_PYTHON_VIRTUALENV_EXECUTABLE}")

  dune_add_test(SOURCES vtktest.cc
                LINK_LIBRARIES ${DUNE_LIBS} Python3::Python
                MPI_RANKS 1 2
                TIMEOUT 1200
                COMPILE_DEFINITIONS PYTHON_INTERPRETER=L"${DUNE_PYTHON_VIRTUALENV_EXECUTABLE}")

  dune_add_test(SOURCES vtksequencetest.cc
                LINK_LIBRARIES ${DUNE_LIBS} Python3::Python
                COMPILE_DEFINITIONS PYTHON_INTERPRETER=L"${DUNE_PYTHON_VIRTUALENV_EXECUTABLE}")
endif()

dune_add_test(SOURCES starcdreadertest.cc
              COMPILE_DEFINITIONS DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\"
              CMAKE_GUARD dune-uggrid_FOUND)

# the gmsh tests
dune_add_test(NAME gmshtest-onedgrid
              SOURCES gmshtest.cc
              COMPILE_DEFINITIONS GMSH_ONEDGRID
                                  DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")

dune_add_test(NAME gmshtest-uggrid
              SOURCES gmshtest.cc
              COMPILE_DEFINITIONS GMSH_UGGRID
                                  DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\"
              CMAKE_GUARD dune-uggrid_FOUND)

if(Alberta_FOUND)
  add_executable(gmshtest-alberta2d gmshtest.cc)
  add_dune_alberta_flags(gmshtest-alberta2d WORLDDIM 2)
  target_compile_definitions(gmshtest-alberta2d PUBLIC GMSH_ALBERTAGRID DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
  dune_add_test(TARGET gmshtest-alberta2d)

  add_executable(gmshtest-alberta3d gmshtest.cc)
  add_dune_alberta_flags(gmshtest-alberta3d WORLDDIM 3)
  target_compile_definitions(gmshtest-alberta3d PUBLIC GMSH_ALBERTAGRID DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
  dune_add_test(TARGET gmshtest-alberta3d)
endif()
