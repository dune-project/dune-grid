if(ALBERTA_FOUND)
  set(ALBERTA_PROGRAMS
    test_alberta_1_1
    test_alberta_1_2
    test_alberta_2_2
    test_alberta_1_3
    test_alberta_2_3
    test_alberta_3_3
    test_alberta_generic)
  set(ALBERTA_EXTRA_PRGRAMS test_alberta)
  message(AUTHOR_WARNING "TODO: Check whether test_alberta is treated correctly.")
 endif(ALBERTA_FOUND)

if(ALUGRID_FOUND)
  set(ALUGRID_PROGRAMS test_alugrid)
endif(ALUGRID_FOUND)

if(UG_FOUND)
  set(UG_PROGRAMS
    test-ug
    test-parallel-ug)
endif(UG_FOUND)

if(UG_FOUND AND ALUGRID_FOUND)
  set(DGFALUGRID_UG_PROGRAMS test-dgfalu-uggrid-combination)
endif(UG_FOUND AND ALUGRID_FOUND)
set(GRIDDIM 2)
set(WORLDDIM 2)
set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" "WORLDDIM=${GRIDDIM}" "YASPGRID")

#
## Some parameters for the geogrid test
#
set(COORDFUNCTION "IdenticalCoordFunction<double,${GRIDDIM}>")
# set(COORDFUNCTION "Helix")
# set(COORDFUNCTION  "ThickHelix")
# set(COORDFUNCTION "DGFCoordFunction\\<${GRIDDIM},${GRIDDIM}\\>")
# set(COORDFUNCTION "DGFCoordFunction\\<2,3\\>")
set(CACHECOORDFUNCTION 0)

#
## define the lists of tests to build and run
#

# tests where program to build and program to run are equal
set(TESTS test_sgrid test_oned test_yaspgrid test_geogrid
  ${ALBERTA_PROGRAMS} ${UG_PROGRAMS} ${ALUGRID_PROGAMS}
  ${DGFALUGRID_UG_PROGRAMS})

set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS "DUNE_DEVEL_MODE")
set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS "DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${CMAKE_SOURCE_DIR}/doc/grids/\"")

get_property(props DIRECTORY PROPERTY COMPILE_DEFINITIONS)
message("props=${props}")
add_executable(test_sgrid test-sgrid.cc)
add_executable(test_oned test-oned.cc)
add_executable(test_yaspgrid test-yaspgrid.cc)
add_dune_mpi_flags(test_yaspgrid)

if(ALBERTA_FOUND)
  add_executable(test_alberta test-alberta.cc)
  add_dune_alberta_flags(test_alberta WORLDDIM ${GRIDDIM})
  add_executable(test_alberta_1_1 test-alberta.cc)
  add_dune_alberta_flags(test_alberta_1_1 GRIDDIM 1 WORLDDIM 1)
  add_executable(test_alberta_1_2 test-alberta.cc)
  add_dune_alberta_flags(test_alberta_1_2 GRIDDIM 1 WORLDDIM 2)
  add_executable(test_alberta_2_2 test-alberta.cc)
  add_dune_alberta_flags(test_alberta_2_2 GRIDDIM 2 WORLDDIM 2)
  add_executable(test_alberta_1_3 test-alberta.cc)
  add_dune_alberta_flags(test_alberta_1_3 GRIDDIM 1 WORLDDIM 3)
  add_executable(test_alberta_2_3 test-alberta.cc)
  add_dune_alberta_flags(test_alberta_2_3 GRIDDIM 2 WORLDDIM 3)
  add_executable(test_alberta_3_3 test-alberta.cc)
  add_dune_alberta_flags(test_alberta_3_3 GRIDDIM 3 WORLDDIM 3)
  add_executable(test_alberta_generic test-alberta.cc)
  add_dune_alberta_flags(test_alberta_generic USE_GENERIC)
endif(ALBERTA_FOUND)

if(ALUGRID_FOUND)
  add_executable(test_alugrid test-alugrid.cc)
  add_dune_alugrid_flags(test_alugrid)
endif(ALUGRID_FOUND)

add_executable(test_geogrid test-geogrid.cc)
set_property(TARGET test_geogrid APPEND PROPERTY COMPILE_DEFINITIONS
  COORDFUNCTION=${COORDFUNCTION} CACHECOORDFUNCTION=${CACHECOORDFUNCTION})

if(UG_FOUND)
  add_executable(test_ug test-ug.cc)
  add_dune_ug_flags(test_ug)
  add_executable(test_parallel_ug test-parallel-ug.cc)
  add_dune_ug_flags(test_parallel_ug)
endif(UG_FOUND)

if(ALUGRID_FOUND AND UG_FOUND)
  add_executable(test_dgfalu_uggrid_combination test-dgfalu-uggrid-combination.cc)
  add_dune_ug_flags(test_dgfalu_uggrid_combination)
  add_dune_alugrid_flags(test_dgfalu_uggrid_combination)
endif(ALUGRID_FOUND AND UG_FOUND)

foreach(_exe ${TESTS})
  target_link_libraries(${_exe} "dunegrid" ${DUNE_LIBS})
  add_test(${_exe} ${_exe})
endforeach(_exe ${TESTS})


set(SOURCES gridcheck.cc staticcheck.hh checkindexset.cc checkgeometryinfather.cc checkintersectionit.cc checkcommunicate.cc checkiterators.cc checktwists.cc check-albertareader.cc basicunitcube.hh checkadaptation.cc checkgeometry.cc checkpartition.cc checkentityseed.cc)
message(AUTHOR_WARNING "Make sure that ${SOURCES} are shipped")

set(CLEANFILES *.gcda *.gcno semantic.cache simplex-testgrid*.dgf.* cube-testgrid*.dgf.* dgfparser.log)
message(AUTHOR_WARNING "TODO: Make sure that ${CLEANFILES} are cleaned during make clean")
