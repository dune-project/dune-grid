add_subdirectory(partitioner)
add_subdirectory(partitioning)
add_subdirectory(test EXCLUDE_FROM_ALL)
set(HEADERS
  entitycommhelper.hh
  entityrange.hh
  entityset.hh
  grapedataioformattypes.hh
  gridinfo-gmsh-main.hh
  gridinfo.hh
  gridtype.hh
  hierarchicsearch.hh
  hostgridaccess.hh
  iterableentityset.hh
  iteratoradapters.hh
  persistentcontainer.hh
  persistentcontainerinterface.hh
  persistentcontainermap.hh
  persistentcontainervector.hh
  persistentcontainerwrapper.hh
  structuredgridfactory.hh
  tripartit.hh
  vertexorderfactory.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/utility)
