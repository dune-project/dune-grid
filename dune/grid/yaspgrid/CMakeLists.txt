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