set(HEADERS
  generic2dune.hh
  numbering.cc
  numbering.hh
  sgrid.cc)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/sgrid/)

set(EXTRA_DIST sgridclasses.fig sgridclasses.eps sgridclasses.png)
message(AUTHOR_WARNING "TODO: Make sure that ${EXTRA_DIST} get distributed")
