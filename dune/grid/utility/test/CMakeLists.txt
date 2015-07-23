set(TESTS
  globalindexsettest
  tensorgridfactorytest
  structuredgridfactorytest
  vertexordertest
  persistentcontainertest)

foreach(_T ${TESTS})
  add_executable(${_T} ${_T}.cc)
  target_link_libraries(${_T} dunegrid ${DUNE_LIBS})
  add_test(${_T} ${_T})
endforeach(_T ${TESTS})

if(MPI_FOUND)
  add_test(NAME structuredgridfactorytest-parallel
      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} "2" structuredgridfactorytest)
  add_test(NAME tensorgridfactorytest-parallel
      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} "2" tensorgridfactorytest)
endif(MPI_FOUND)

add_dune_ug_flags(${TESTS})
add_dune_mpi_flags(structuredgridfactorytest tensorgridfactorytest)
add_dune_alugrid_flags(vertexordertest persistentcontainertest)

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS})
