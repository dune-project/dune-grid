set(HEADERS
  coordinates.hh
  factory.hh
  sequences.hh
  torus.hh
  yaspgridentity.hh
  yaspgridentitypointer.hh
  yaspgridentityseed.hh
  yaspgridgeometry.hh
  yaspgridhierarchiciterator.hh
  yaspgridindexsets.hh
  yaspgridintersection.hh
  yaspgridintersectioniterator.hh
  yaspgrididset.hh
  yaspgridleveliterator.hh
  ygrid.hh)

exclude_all_but_from_headercheck(torus.hh coordinates.hh ygrid.hh factory.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/yaspgrid/)
