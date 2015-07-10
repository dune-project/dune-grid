
# add test include property manually because the folder is not called "test"
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  PROPERTY TEST_INCLUDE_FILE ${CMAKE_CURRENT_BINARY_DIR}/BuildTests.cmake)

set(TESTS
  test-yaspgrid-backuprestore-equidistant
  test-yaspgrid-backuprestore-equidistantoffset
  test-yaspgrid-backuprestore-tensor
  test-yaspgrid-tensorgridfactory
  test-yaspgrid-yaspfactory-1d
  test-yaspgrid-yaspfactory-2d
  test-yaspgrid-yaspfactory-3d
)

# Add tests to test target
add_directory_test_target(_test_target)

# Add all tests, do not check deprecation
foreach(_test ${TESTS})
  add_executable(${_test} EXCLUDE_FROM_ALL ${_test}.cc)
  add_dune_mpi_flags(${_test})

  target_compile_definitions(${_test}
    PUBLIC
    DISABLE_DEPRECATED_METHOD_CHECK=1
  )
  target_compile_options(${_test}
    PUBLIC
    -Wall
    -Werror
  )

  target_link_libraries(${_test} "dunegrid" ${DUNE_LIBS})
  add_test(${_test} ${_test})
  add_dependencies(${_test_target} ${_test})
endforeach(_test ${TESTS})

# Add all tests again, this time with deprecation check
foreach(_test ${TESTS})
  add_executable(${_test}_dep EXCLUDE_FROM_ALL ${_test}.cc)
  add_dune_mpi_flags(${_test}_dep)

  target_compile_options(${_test}_dep
    PUBLIC
    -Wall
    -Wno-deprecated-declarations
    -Werror
  )

  target_link_libraries(${_test}_dep "dunegrid" ${DUNE_LIBS})
  add_test(${_test}_dep ${_test}_dep)
  add_dependencies(${_test_target} ${_test}_dep)
endforeach(_test ${TESTS})


# If we have MPI, include parallel tests
if(MPI_FOUND)

  set(TESTS_PARALLEL
    test-yaspgrid-backuprestore
    test-yaspgrid-tensorgridfactory
    test-yaspgrid-yaspfactory-1d
    test-yaspgrid-yaspfactory-2d
    test-yaspgrid-yaspfactory-3d
  )

  # Create tests with mpi command and add them to the test target
  foreach(_test ${TESTS_PARALLEL})
    add_test(NAME ${_test}_parallel
      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} "2" ${_test})
    add_dependencies(${_test_target} ${_test}_parallel)
  endforeach(_test ${TESTS_PARALLEL})

endif(MPI_FOUND)

set(SOURCES
  test-yaspgrid-backuprestore.cc
  test-yaspgrid-tensorgridfactory.cc
  test-yaspgrid-yaspfactory-1d.cc
  test-yaspgrid-yaspfactory-2d.cc
  test-yaspgrid-yaspfactory-3d.cc
  test-yaspgrid.hh
)
install(FILES ${SOURCES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/test)
