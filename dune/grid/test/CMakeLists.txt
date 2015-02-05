if(ALBERTA_FOUND)
  set(ALBERTA_PROGRAMS
    test-alberta-1-1
    test-alberta-1-2
    test-alberta-2-2
    test-alberta-1-3
    test-alberta-2-3
    test-alberta-3-3
    test-alberta-generic)
endif(ALBERTA_FOUND)

if(ALUGRID_FOUND)
  set(ALUGRID_PROGRAMS test-alugrid)
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
set(TESTS
  test-geogrid test-identitygrid test-oned test-sgrid test-yaspgrid
  ${ALBERTA_PROGRAMS} ${ALUGRID_PROGRAMS} ${UG_PROGRAMS}
  ${DGFALUGRID_UG_PROGRAMS} test-mcmg-geogrid
  testiteratorranges)

set_property(DIRECTORY APPEND PROPERTY
  COMPILE_DEFINITIONS "DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\"")

add_executable(test-geogrid EXCLUDE_FROM_ALL test-geogrid.cc)
add_executable(test-identitygrid EXCLUDE_FROM_ALL test-identitygrid.cc)
add_executable(test-oned EXCLUDE_FROM_ALL test-oned.cc)
add_executable(test-mcmg-geogrid EXCLUDE_FROM_ALL test-mcmg-geogrid.cc)
add_executable(test-sgrid EXCLUDE_FROM_ALL test-sgrid.cc)
add_executable(test-yaspgrid EXCLUDE_FROM_ALL test-yaspgrid.cc)
add_dune_mpi_flags(test-yaspgrid)

set_property(TARGET test-sgrid test-oned test-yaspgrid test-geogrid APPEND
  PROPERTY COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" "WORLDDIM=${GRIDDIM}" "YASPGRID")

set_property(TARGET test-geogrid APPEND PROPERTY COMPILE_DEFINITIONS
  COORDFUNCTION=${COORDFUNCTION} CACHECOORDFUNCTION=${CACHECOORDFUNCTION})

if(ALBERTA_FOUND)
  add_executable(test-alberta EXCLUDE_FROM_ALL test-alberta.cc)
  add_dune_alberta_flags(test-alberta WORLDDIM ${GRIDDIM})
  add_executable(test-alberta-1-1 EXCLUDE_FROM_ALL test-alberta.cc)
  add_dune_alberta_flags(test-alberta-1-1 GRIDDIM 1 WORLDDIM 1)
  add_executable(test-alberta-1-2 EXCLUDE_FROM_ALL test-alberta.cc)
  add_dune_alberta_flags(test-alberta-1-2 GRIDDIM 1 WORLDDIM 2)
  add_executable(test-alberta-2-2 EXCLUDE_FROM_ALL test-alberta.cc)
  add_dune_alberta_flags(test-alberta-2-2 GRIDDIM 2 WORLDDIM 2)
  add_executable(test-alberta-1-3 EXCLUDE_FROM_ALL test-alberta.cc)
  add_dune_alberta_flags(test-alberta-1-3 GRIDDIM 1 WORLDDIM 3)
  add_executable(test-alberta-2-3 EXCLUDE_FROM_ALL test-alberta.cc)
  add_dune_alberta_flags(test-alberta-2-3 GRIDDIM 2 WORLDDIM 3)
  add_executable(test-alberta-3-3 EXCLUDE_FROM_ALL test-alberta.cc)
  add_dune_alberta_flags(test-alberta-3-3 GRIDDIM 3 WORLDDIM 3)
  add_executable(test-alberta-generic EXCLUDE_FROM_ALL test-alberta.cc)
  add_dune_alberta_flags(test-alberta-generic USE_GENERIC GRIDDIM 2)
endif(ALBERTA_FOUND)

if(ALUGRID_FOUND)
  add_executable(test-alugrid EXCLUDE_FROM_ALL test-alugrid.cc)
  add_dune_alugrid_flags(test-alugrid)
endif(ALUGRID_FOUND)

if(UG_FOUND)
  add_executable(test-ug EXCLUDE_FROM_ALL test-ug.cc)
  add_dune_ug_flags(test-ug)
  add_executable(test-parallel-ug EXCLUDE_FROM_ALL test-parallel-ug.cc)
  add_dune_ug_flags(test-parallel-ug)
  # test-loadbalancing
  add_executable(test-loadbalancing EXCLUDE_FROM_ALL test-loadbalancing.cc)
  add_dune_ug_flags(test-loadbalancing)
  add_dune_parmetis_flags(test-loadbalancing)
endif(UG_FOUND)

if(ALUGRID_FOUND AND UG_FOUND)
  add_executable(test-dgfalu-uggrid-combination EXCLUDE_FROM_ALL test-dgfalu-uggrid-combination.cc)
  add_dune_ug_flags(test-dgfalu-uggrid-combination)
  add_dune_alugrid_flags(test-dgfalu-uggrid-combination)
endif(ALUGRID_FOUND AND UG_FOUND)

add_executable(testiteratorranges EXCLUDE_FROM_ALL testiteratorranges.cc)
add_dune_mpi_flags(testiteratorranges)


# tell gridcheck to use old entity and intersections interface for unported grids
set_property(TARGET
  test-geogrid
  test-identitygrid
  test-mcmg-geogrid
  ${ALBERTA_PROGRAMS}
  ${ALUGRID_PROGRAMS}
  APPEND PROPERTY COMPILE_DEFINITIONS "DUNE_GRID_CHECK_USE_DEPRECATED_ENTITY_AND_INTERSECTION_INTERFACE"
  )


foreach(_exe ${TESTS})
  target_link_libraries(${_exe} "dunegrid" ${DUNE_LIBS})
  add_test(${_exe} ${_exe})
endforeach(_exe ${TESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS})

# define HAVE_DUNE_GRID for the dgfparser
set_property(TARGET ${TESTS} APPEND PROPERTY
  COMPILE_DEFINITIONS HAVE_DUNE_GRID=1)

set(SOURCES
  basicunitcube.hh
  check-albertareader.cc
  checkadaptation.cc
  checkadaptation.hh
  checkcommunicate.cc
  checkcommunicate.hh
  checkentityseed.cc
  checkentityseed.hh
  checkgeometry.cc
  checkgeometryinfather.cc
  checkindexset.cc
  checkintersectionit.cc
  checkiterators.cc
  checkpartition.cc
  checktwists.cc
  functions.hh
  gridcheck.cc
  staticcheck.hh)

# install the tests as we want to support testing 3rdparty grids with installed dune-grid
install(FILES ${SOURCES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/test)

# clean up additional files during make clean. Unfortunately it is not possible with
# cmake to use search patterns here. Automake wanted to clean the following files:
# *.gcda *.gcno semantic.cache simplex-testgrid*.dgf.* cube-testgrid*.dgf.* dgfparser.log
# add files to the list as they appear.
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "semantic.cache;dgfparser.log")
