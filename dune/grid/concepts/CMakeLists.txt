set(HEADERS
  entity.hh
  entityiterator.hh
  geometry.hh
  grid.hh
  gridview.hh
  indexset.hh
  intersection.hh
  intersectioniterator.hh)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/concepts)