add_subdirectory(partitioner)
add_subdirectory(partitioning)
add_subdirectory(test)
set(HEADERS
  entitycommhelper.hh
  entityrange.hh
  entityset.hh
  globalindexset.hh
  gridinfo-gmsh-main.hh
  gridinfo.hh
  gridtype.hh
  hierarchicsearch.hh
  hostgridaccess.hh
  iterableentityset.hh
  iteratoradapters.hh
  multiindex.hh
  parmetisgridpartitioner.hh
  persistentcontainer.hh
  persistentcontainerinterface.hh
  persistentcontainermap.hh
  persistentcontainervector.hh
  persistentcontainerwrapper.hh
  structuredgridfactory.hh
  tensorgridfactory.hh
  vertexorderfactory.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/utility)
