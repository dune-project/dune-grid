add_subdirectory(test EXCLUDE_FROM_ALL)
set(HEADERS
  grapedataioformattypes.hh
  gridinfo-gmsh-main.hh
  gridinfo.hh
  gridtype.hh
  hierarchicsearch.hh
  hostgridaccess.hh
  persistentcontainer.hh
  persistentcontainerwrapper.hh
  structuredgridfactory.hh
  vertexorderfactory.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/utility})
