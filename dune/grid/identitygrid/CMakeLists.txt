set(HEADERS
  identitygridentity.hh
  identitygridentityseed.hh
  identitygridhierarchiciterator.hh
  identitygridgeometry.hh
  identitygridintersectioniterator.hh
  identitygridintersections.hh
  identitygridindexsets.hh
  identitygridleafiterator.hh
  identitygridleveliterator.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/identitygrid)
