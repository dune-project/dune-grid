set(GRIDDIM 2)
set(GRIDTYPE YASPGRID)
add_definitions("-DDUNE_GRID_EXAMPLE_GRIDS_PATH=\"${CMAKE_SOURCE_DIR}/doc/grids/\"")

set(TESTS
  gmshtest
  gnuplottest)

set(CONSISTENT_VTK_TESTS
  conformvolumevtktest
  nonconformboundaryvtktest
  vtktest
  vtksequencetest)

set(VTK_TESTS ${CONSISTENT_VTK_TEST} subsamplingvtktest)

if(AMIRAMESH_FOUND AND PSURFACE_FOUND)
  set(AMIRAMESH_TESTS amirameshtest)
endif(AMIRAMESH_FOUND AND PSURFACE_FOUND)

if(ALBERTA_FOUND)
  set(ALBERTA_TESTS gmshtest_alberta2d gmshtest_alberta3d)
endif(ALBERTA_FOUND)

if(ALUGRID_FOUND)
  set(ALUGRID_TESTS gmshtest_alugrid)
endif(ALUGRID_FOUND)

if(UG_FOUND)
  set(UG_TESTS starcdreadertest)
endif(UG_FOUND)

set(BUILD_TESTS ${TESTS} ${UG_TESTS} ${CONSISTENT_VTK_TESTS})
set(ALLTESTS ${BUILD_TESTS}  ${ALUGRID_TESTS} ${ALBERTA_TESTS})

foreach(_test ${BUILD_TESTS})
  add_executable(${_test} ${_test}.cc)
  target_link_libraries(${_test} dunegrid ${DUNE_LIBS})
endforeach(_test ${BUILD_TESTS})

add_executable(subsamplingvtktest subsamplingvtktest.cc test-linking.cc)

foreach(_test ${ALLTESTS})
  add_test(${_test} ${_test})
endforeach(_test ${ALLTESTS})

if(MPI_FOUND)
  add_test(NAME mpivtktest WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND mpirun -np 2 ./vtktest)
endif(MPI_FOUND)

foreach(_test ${AMIRAMESH_TESTS})
  add_dune_amiramesh_flags(${_test})
  add_dune_ug_flags(${_test})
  add_dune_alugrid_flags(${_test})
  add_dune_alberta_flags(${_test} GRIDDIM 2)
  add_dune_psurface_flags(${_test})
endforeach(_test ${AMIRAMESH_TESTS})

add_dune_mpi_flags(${VTK_TESTS})

add_executable(gmshtest_alugrid gmshtest.cc)
add_dune_alugrid_flags(gmshtest_alugrid)

add_executable(gmshtest_alberta2d gmshtest.cc)
add_dune_alberta_flags(gmshtest_alberta2d WORLDDIM 2)
add_dune_grape_flags(gmshtest_alberta2d)
set_property(TARGET gmshtest_alberta2d APPEND_STRING
  PROPERTY COMPILE_FLAGS "-DGRIDDIM=${ALBERTA_DIM}")

add_executable(gmshtest_alberta3d gmshtest.cc)
add_dune_alberta_flags(gmshtest_alberta3d WORLDDIM 3)
add_dune_grape_flags(gmshtest_alberta3d)

foreach(_test ${UG_TESTS})
  add_dune_ug_flags(${_test})
endforeach(_test ${UG_TESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})
