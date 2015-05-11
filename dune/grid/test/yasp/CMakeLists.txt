
set(TESTS
  test-yaspgrid-backuprestore
  test-yaspgrid-tensorgridfactory
  test-yaspgrid-yaspfactory-1d
  test-yaspgrid-yaspfactory-2d
  test-yaspgrid-yaspfactory-3d
)

# Add all tests
foreach(_test ${TESTS})
  add_executable(${_test} EXCLUDE_FROM_ALL ${_test}.cc)
  add_dune_mpi_flags(${_test})

  target_link_libraries(${_test} "dunegrid" ${DUNE_LIBS})
  add_test(${_test} ${_test})
endforeach(_test ${TESTS})

# add test include property manually because the folder is not called "test"
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  PROPERTY TEST_INCLUDE_FILE ${CMAKE_CURRENT_BINARY_DIR}/BuildTests.cmake)

# Add tests to test target
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS})


set(SOURCES
  test-yaspgrid-backuprestore.cc
  test-yaspgrid-tensorgridfactory.cc
  test-yaspgrid-yaspfactory-1d.cc
  test-yaspgrid-yaspfactory-2d.cc
  test-yaspgrid-yaspfactory-3d.cc
  test-yaspgrid.hh
)
install(FILES ${SOURCES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/test)
