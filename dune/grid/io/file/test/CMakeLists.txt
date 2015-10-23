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
              SOURCES subsamplingvtktest.cc test-linking.cc)

dune_add_test(SOURCES vtktest.cc
              MPI_RANKS 1 2
              SKIP_ON_77)

dune_add_test(SOURCES vtksequencetest.cc)

if(AMIRAMESH_FOUND)
  dune_add_test(SOURCES amirameshtest.cc)
endif()

if(UG_FOUND)
  dune_add_test(SOURCES starcdreadertest.cc
                COMPILE_DEFINITIONS DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
endif()

# the gmsh tests
dune_add_test(NAME gmshtest-onedgrid
              SOURCES gmshtest.cc
              COMPILE_DEFINITIONS GMSH_ONEDGRID
                                  DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")

if(HAVE_UG)
  dune_add_test(NAME gmshtest-uggrid
                SOURCES gmshtest.cc
                COMPILE_DEFINITIONS GMSH_UGGRID
                                    DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
endif()

if(HAVE_ALUGRID)
  dune_add_test(NAME gmshtest-alugrid
                SOURCES gmshtest.cc
                COMPILE_DEFINITIONS GMSH_ALUGRID
                                    DUNE_GRID_CHECK_USE_DEPRECATED_ENTITY_AND_INTERSECTION_INTERFACE
                                    DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
endif()

if(HAVE_ALBERTA)
  add_executable(gmshtest-alberta2d gmshtest.cc)
  add_dune_alberta_flags(gmshtest-alberta2d WORLDDIM 2)
  target_compile_definitions(gmshtest-alberta2d PUBLIC GMSH_ALBERTAGRID DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
  dune_add_test(TARGET gmshtest-alberta2d)

  add_executable(gmshtest-alberta3d gmshtest.cc)
  add_dune_alberta_flags(gmshtest-alberta3d WORLDDIM 3)
  target_compile_definitions(gmshtest-alberta3d PUBLIC GMSH_ALBERTAGRID DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
  dune_add_test(TARGET gmshtest-alberta3d)
endif()
