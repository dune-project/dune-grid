add_subdirectory(test EXCLUDE_FROM_ALL)

set(HEADERS
  equidistant.hh
  recursive-equidistant.hh
  recursive-odd.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/utility/partitioner)
