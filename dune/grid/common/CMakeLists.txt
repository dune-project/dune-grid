
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
  exceptions.hh
  geometry.hh
  grid.hh
  gridenums.hh
  gridfactory.hh
  gridinfo.hh
  gridview.hh
  hierarchiciterator.hh
  indexidset.hh
  intersection.hh
  intersectioniterator.hh
  leafiterator.hh
  leveliterator.hh
  mcmgmapper.hh
  mapper.hh
  refinement.hh
  sizecache.hh
  scsgmapper.hh
  universalmapper.hh
  virtualrefinement.hh
  virtualrefinement.cc)

add_subdirectory(refinement)
add_subdirectory(test EXCLUDE_FROM_ALL)

set(NOINST_DATA cube-to-tet-6.png)

message(AUTHOR_WARNING "Check wether ${NOINST_DATA} ends up in the tarball")

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/common)
