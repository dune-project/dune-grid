set(HEADERS
  amirameshreader.cc
  amirameshwriter.cc
  psurfaceboundary.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/io/file/amiramesh)
