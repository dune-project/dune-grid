set(HEADERS
  grids.hh
  yaspgridentity.hh
  yaspgridentitypointer.hh
  yaspgridentityseed.hh
  yaspgridgeometry.hh
  yaspgridhierarchiciterator.hh
  yaspgridindexsets.hh
  yaspgridintersection.hh
  yaspgridintersectioniterator.hh
  yaspgrididset.hh
  yaspgridleveliterator.hh)

exclude_all_but_from_headercheck(grids.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/yaspgrid/)

set(EXTRA_DIST grid.fig grid.eps grid.png subgrid.fig subgrid.eps subgrid.png)
message(AUTHOR_WARNING "TODO: Make sure that ${EXTRA_DIST} get distributed")
