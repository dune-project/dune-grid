set(HEADERS
  generic2dune.hh
  numbering.cc
  numbering.hh
  sgrid.cc)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/sgrid/)