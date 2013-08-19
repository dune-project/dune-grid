add_subdirectory(test EXCLUDE_FROM_ALL)

set(HEADERS
  interface.hh
  strided.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/utility/gridviewpartitioning)
