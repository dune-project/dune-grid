set(HEADERS
  adaptcallback.hh
  backuprestore.hh
  boundaryprojection.hh
  boundarysegment.hh
  capabilities.hh
  datahandleif.hh
  defaultgridview.hh
  entity.hh
  entityiterator.hh
  entitypointer.hh
  entityrangegenerators.hh
  entityseed.hh
  exceptions.hh
  geometry.hh
  grid.hh
  gridenums.hh
  gridfactory.hh
  gridinfo.hh
  gridview.hh
  indexidset.hh
  intersection.hh
  intersectioniterator.hh
  mcmgmapper.hh
  mapper.hh
  partitionset.hh
  sizecache.hh
  scsgmapper.hh
  universalmapper.hh)

add_subdirectory(test EXCLUDE_FROM_ALL)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/common)
