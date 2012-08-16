add_subdirectory(grape)

set(HEADERS
  grapedatadisplay.hh
  grapegriddisplay.hh
  combinedgrapedisplay.hh)
install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/io/visual)
