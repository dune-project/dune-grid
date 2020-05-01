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

dune_add_test(SOURCES amirameshtest.cc
              CMAKE_GUARD AMIRAMESH_FOUND)

dune_add_test(SOURCES starcdreadertest.cc
              COMPILE_DEFINITIONS DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\"
              CMAKE_GUARD dune-uggrid_FOUND)

# the gmsh tests
foreach(FLAVOUR IN ITEMS dynamicgmsh gmsh)
  dune_add_test(NAME ${FLAVOUR}test-onedgrid
                SOURCES ${FLAVOUR}test.cc
                LINK_LIBRARIES dunegrid
                COMPILE_DEFINITIONS GMSH_ONEDGRID
                                    DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")

  dune_add_test(NAME ${FLAVOUR}test-uggrid
                SOURCES ${FLAVOUR}test.cc
                COMPILE_DEFINITIONS GMSH_UGGRID
                                    DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\"
                CMAKE_GUARD dune-uggrid_FOUND)

  if(HAVE_ALBERTA)
    foreach(DIM IN ITEMS 2 3)
      add_executable(${FLAVOUR}test-alberta${DIM}d ${FLAVOUR}test.cc)
      add_dune_alberta_flags(${FLAVOUR}test-alberta${DIM}d WORLDDIM ${DIM})
      target_compile_definitions(${FLAVOUR}test-alberta${DIM}d PUBLIC GMSH_ALBERTAGRID DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
      dune_add_test(TARGET ${FLAVOUR}test-alberta${DIM}d)
    endforeach(DIM)
  endif()
endforeach(FLAVOUR)
