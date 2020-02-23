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
  rangegenerators.hh
  sizecache.hh
  scsgmapper.hh)

add_subdirectory(test)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/common)
