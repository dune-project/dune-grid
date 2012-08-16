set(HEADERS
  grids.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/yaspgrid/)

set(EXTRA_DIST grid.fig grid.eps grid.png subgrid.fig subgrid.eps subgrid.png)
message(AUTHOR_WARNING "TODO: Make sure that ${EXTRA_DIST} get distributed")
