dune_add_test(SOURCES test-geogrid.cc
              COMPILE_DEFINITIONS COORDFUNCTION=IdenticalCoordFunction<double,2>
                                  CACHECOORDFUNCTION=0
                                  DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\"
                                  GRIDDIM=2
                                  WORLDDIM=2
                                  YASPGRID)

dune_add_test(SOURCES test-identitygrid.cc)

# TODO: WTF? I need to set GRIDDIM to 2 to compile OneDGrid????????
dune_add_test(SOURCES test-oned.cc
              COMPILE_DEFINITIONS GRIDDIM=2 WORLDDIM=2)

dune_add_test(SOURCES test-mcmg-geogrid.cc)

dune_add_test(SOURCES test-yaspgrid.cc
              COMPILE_DEFINITIONS GRIDDIM=2 WORLDDIM=2)

dune_add_test(SOURCES testiteratorranges.cc)

dune_add_test(SOURCES test-hierarchicsearch.cc)


if(ALUGRID_FOUND)
  dune_add_test(SOURCES test-alugrid.cc
                COMPILE_DEFINITIONS DUNE_GRID_CHECK_USE_DEPRECATED_ENTITY_AND_INTERSECTION_INTERFACE)
endif(ALUGRID_FOUND)

if(UG_FOUND)
  dune_add_test(SOURCES test-ug.cc)

  dune_add_test(SOURCES test-parallel-ug.cc)

  dune_add_test(SOURCES test-loadbalancing.cc)
endif(UG_FOUND)

if(ALUGRID_FOUND AND UG_FOUND)
  dune_add_test(SOURCES test-dgfalu-uggrid-combination.cc)
endif(ALUGRID_FOUND AND UG_FOUND)

# The alberta tests are only alibi-ported, until the  grid and world dimension
# are configuretime parameters and we can treat alberta just as any other grid manager
# - buildsystemwise. PLEASE DONT LOOK AT THIS IF YOU WANT TO KNOW HOW TO WRITE TESTS.
if(ALBERTA_FOUND)
  add_executable(test-alberta EXCLUDE_FROM_ALL test-alberta.cc)
  add_dune_alberta_flags(test-alberta WORLDDIM ${GRIDDIM})
  target_compile_definitions(test-alberta PUBLIC DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
  dune_add_test(TARGET test-alberta)

  add_executable(test-alberta-1-1 EXCLUDE_FROM_ALL test-alberta.cc)
  add_dune_alberta_flags(test-alberta-1-1 GRIDDIM 1 WORLDDIM 1)
  target_compile_definitions(test-alberta-1-1 PUBLIC DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
  dune_add_test(TARGET test-alberta-1-1)

  add_executable(test-alberta-1-1-no-deprecated EXCLUDE_FROM_ALL test-alberta.cc)
  add_dune_alberta_flags(test-alberta-1-1-no-deprecated GRIDDIM 1 WORLDDIM 1)
  target_compile_definitions(test-alberta-1-1-no-deprecated PUBLIC
    DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\"
    DISABLE_DEPRECATED_METHOD_CHECK=1)
  dune_add_test(TARGET test-alberta-1-1-no-deprecated)

  add_executable(test-alberta-1-2 EXCLUDE_FROM_ALL test-alberta.cc)
  add_dune_alberta_flags(test-alberta-1-2 WORLDDIM 2)
  target_compile_definitions(test-alberta-1-2 PUBLIC DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
  dune_add_test(TARGET test-alberta-1-2)

  add_executable(test-alberta-2-2 EXCLUDE_FROM_ALL test-alberta.cc)
  add_dune_alberta_flags(test-alberta-2-2 WORLDDIM 2)
  target_compile_definitions(test-alberta-2-2 PUBLIC DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
  dune_add_test(TARGET test-alberta-2-2)

  add_executable(test-alberta-1-3 EXCLUDE_FROM_ALL test-alberta.cc)
  add_dune_alberta_flags(test-alberta-1-3 WORLDDIM 3)
  target_compile_definitions(test-alberta-1-3 PUBLIC DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
  dune_add_test(TARGET test-alberta-1-3)

  add_executable(test-alberta-2-3 EXCLUDE_FROM_ALL test-alberta.cc)
  add_dune_alberta_flags(test-alberta-2-3 WORLDDIM 3)
  target_compile_definitions(test-alberta-2-3 PUBLIC DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
  dune_add_test(TARGET test-alberta-2-3)

  add_executable(test-alberta-3-3 EXCLUDE_FROM_ALL test-alberta.cc)
  add_dune_alberta_flags(test-alberta-3-3 WORLDDIM 3)
  target_compile_definitions(test-alberta-3-3 PUBLIC DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
  dune_add_test(TARGET test-alberta-3-3)

  add_executable(test-alberta-generic EXCLUDE_FROM_ALL test-alberta.cc)
  add_dune_alberta_flags(test-alberta-generic USE_GENERIC WORLDDIM 2)
  target_compile_definitions(test-alberta-generic PUBLIC DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
  dune_add_test(TARGET test-alberta-generic)
endif(ALBERTA_FOUND)

# install the test tools as we want to support testing 3rdparty grids with installed dune-grid
install(FILES basicunitcube.hh
              check-albertareader.hh
              checkadaptation.hh
              checkcommunicate.hh
              checkentitylifetime.hh
              checkentityseed.hh
              checkgeometry.hh
              checkgeometryinfather.hh
              checkindexset.hh
              checkintersectionit.hh
              checkintersectionlifetime.hh
              checkiterators.hh
              checkjacobians.hh
              checkpartition.hh
              checktwists.hh
              functions.hh
              gridcheck.hh
              staticcheck.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/test)

# clean up additional files during make clean. Unfortunately it is not possible with
# cmake to use search patterns here. Automake wanted to clean the following files:
# *.gcda *.gcno semantic.cache simplex-testgrid*.dgf.* cube-testgrid*.dgf.* dgfparser.log
# add files to the list as they appear.
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "semantic.cache;dgfparser.log")
