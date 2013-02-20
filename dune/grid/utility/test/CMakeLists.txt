set(TESTS
  structuredgridfactorytest
  vertexordertest
  persistentcontainertest)

foreach(_T ${TESTS})
  add_test(${_T} ${_T})
  add_executable(${_T} ${_T}.cc)
endforeach(_T ${TESTS})

add_dune_ug_flags(${TESTS})
add_dune_alugrid_flags(vertexordertest persistentcontainertest)
