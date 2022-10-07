# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_subdirectory(yasp)

dune_add_test(SOURCES geometrygrid-coordfunction-copyconstructor.cc)

dune_add_test(NAME test-geogrid-yaspgrid
              SOURCES test-geogrid.cc
              LINK_LIBRARIES dunegrid
              COMPILE_DEFINITIONS COORDFUNCTION=IdenticalCoordFunction<double,2>
                                  CACHECOORDFUNCTION=0
                                  DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\"
                                  GRIDTYPE=Dune::YaspGrid<2>)

dune_add_test(NAME test-geogrid-uggrid
              SOURCES test-geogrid.cc
              LINK_LIBRARIES dunegrid
              COMPILE_DEFINITIONS COORDFUNCTION=IdenticalCoordFunction<double,2>
                                  CACHECOORDFUNCTION=0
                                  DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\"
                                  GRIDTYPE=Dune::UGGrid<2>
              CMAKE_GUARD dune-uggrid_FOUND)

dune_add_test(SOURCES test-gridinfo.cc)

dune_add_test(SOURCES test-identitygrid.cc)

dune_add_test(SOURCES test-oned.cc
              LINK_LIBRARIES dunegrid)

dune_add_test(SOURCES test-mcmg-geogrid.cc)

dune_add_test(SOURCES testiteratorranges.cc)

dune_add_test(SOURCES test-hierarchicsearch.cc)

dune_add_test(SOURCES test-ug.cc
              LINK_LIBRARIES dunegrid
              COMPILE_DEFINITIONS DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\"
              CMAKE_GUARD dune-uggrid_FOUND)

dune_add_test(SOURCES test-parallel-ug.cc
              LINK_LIBRARIES dunegrid
              CMAKE_GUARD dune-uggrid_FOUND
              MPI_RANKS 1 2 3 4 8
              TIMEOUT 300)

dune_add_test(SOURCES test-loadbalancing.cc
              LINK_LIBRARIES dunegrid
              CMAKE_GUARD dune-uggrid_FOUND)

dune_add_test(SOURCES issue-53-uggrid-intersections.cc
              LINK_LIBRARIES dunegrid
              CMAKE_GUARD dune-uggrid_FOUND)

# The alberta tests are only alibi-ported, until the  grid and world dimension
# are configuretime parameters and we can treat alberta just as any other grid manager
# - buildsystemwise. PLEASE DON'T LOOK AT THIS IF YOU WANT TO KNOW HOW TO WRITE TESTS.
if(Alberta_FOUND)
  foreach(WORLDDIM RANGE 1 3)
    foreach(GRIDDIM RANGE 1 ${WORLDDIM})
      set(_test test-alberta-${GRIDDIM}-${WORLDDIM})
      add_executable(${_test} EXCLUDE_FROM_ALL test-alberta.cc)
      target_link_libraries(${_test} PRIVATE dunegrid)
      add_dune_alberta_flags(${_test} WORLDDIM ${WORLDDIM})
      target_compile_definitions(${_test} PUBLIC
        GRIDDIM=${GRIDDIM}
        DUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")
      dune_add_test(TARGET ${_test})
      unset(_test)
    endforeach(GRIDDIM)
  endforeach(WORLDDIM)
endif()

# install the test tools as we want to support testing 3rdparty grids with installed dune-grid
install(FILES basicunitcube.hh
              check-albertareader.hh
              checkadaptation.hh
              checkcomcorrectness.hh
              checkcommunicate.hh
              checkentitylifetime.hh
              checkentityseed.hh
              checkgeometry.hh
              checkgeometryinfather.hh
              checkgridfactory.hh
              checkidset.hh
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
