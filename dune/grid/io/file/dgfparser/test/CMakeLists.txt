add_definitions(-DDUNE_GRID_EXAMPLE_GRIDS_PATH=\"${PROJECT_SOURCE_DIR}/doc/grids/\")

set(TESTS test-dgf-yasp test-dgf-oned)

if(ALUGRID_FOUND)
  list(APPEND TESTS test-dgf-alu)
endif(ALUGRID_FOUND)

if(ALBERTA_FOUND)
  list(APPEND TESTS test-dgf-alberta)
endif(ALBERTA_FOUND)

if(UG_FOUND)
  list(APPEND TESTS test-dgf-ug)
endif(UG_FOUND)

foreach(_test ${TESTS})
  add_executable(${_test} test-dgf.cc)
  target_link_libraries(${_test} dunegrid ${DUNE_LIBS})
  add_test(${_test} ${_test})
endforeach(_test ${TESTS})

# test-dgf-yasp
set_property(TARGET test-dgf-yasp APPEND PROPERTY
    COMPILE_DEFINITIONS GRIDDIM=3 YASPGRID HAVE_DUNE_GRID=1)
# test-dgf-yasp-offset
add_executable(test-dgf-yasp-offset test-dgf.cc)
target_link_libraries(test-dgf-yasp-offset dunegrid ${DUNE_LIBS})
add_test(NAME test-dgf-yasp-offset
  COMMAND test-dgf-yasp-offset "${PROJECT_SOURCE_DIR}/doc/grids/dgf/test2d_offset.dgf")
set(gridtypeOffset "YaspGrid<2,EquidistantOffsetCoordinates<double,2>>")
set_property(TARGET test-dgf-yasp-offset
  APPEND PROPERTY
  COMPILE_DEFINITIONS GRIDDIM=3 YASPGRID HAVE_DUNE_GRID=1
    COMPLETE_GRID_TYPE=${gridtypeOffset})
# test-dgf-oned
set_property(TARGET test-dgf-oned APPEND PROPERTY
    COMPILE_DEFINITIONS GRIDDIM=1 ONEDGRID HAVE_DUNE_GRID=1)

if(ALUGRID_FOUND)
  add_dune_alugrid_flags(test-dgf-alu)
  set_property(TARGET test-dgf-alu APPEND PROPERTY
    COMPILE_DEFINITIONS ALUGRID_CUBE GRIDDIM=3 HAVE_DUNE_GRID=1)
endif(ALUGRID_FOUND)

if(ALBERTA_FOUND)
  add_dune_alberta_flags(test-dgf-alberta GRIDDIM 2 WORLDDIM 2)
  set_property(TARGET test-dgf-alberta APPEND PROPERTY
    COMPILE_DEFINITIONS ALBERTAGRID GRIDDIM=2  HAVE_DUNE_GRID=1)
endif(ALBERTA_FOUND)

if(UG_FOUND)
  add_dune_ug_flags(test-dgf-ug)
  set_property(TARGET test-dgf-ug APPEND PROPERTY
    COMPILE_DEFINITIONS UGGRID GRIDDIM=3  HAVE_DUNE_GRID=1)
endif(UG_FOUND)

foreach(_test ${TESTS})
  add_dune_mpi_flags(${_test})
endforeach(_test ${TESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} test-dgf-yasp-offset)
