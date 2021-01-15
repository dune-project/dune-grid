add_subdirectory(test)

if(ALBERTA_FOUND)
  set(ALBERTAGRID_SOURCES
    dgfparser.cc
    elementinfo.cc
    indexsets.cc
    macrodata.cc
    meshpointer.cc)
  # create dimension dependent alberta libraries
  foreach(_dim ${ALBERTA_WORLD_DIMS})
    list(APPEND ALBERTAGRID_LIBS albertagrid_${_dim}d)
    dune_add_library(albertagrid_${_dim}d OBJECT SOURCES ${ALBERTAGRID_SOURCES})
  endforeach(_dim RANGE 1 9)
endif(ALBERTA_FOUND)

#install header
set(HEADERS
  agrid.hh
  albertagrid.cc
  albertaheader.hh
  indexsets.hh
  indexstack.hh
  datahandle.hh
  misc.hh
  macroelement.hh
  elementinfo.hh
  geometrycache.hh
  meshpointer.hh
  macrodata.hh
  dofadmin.hh
  dofvector.hh
  refinement.hh
  coordcache.hh
  level.hh
  undefine-2.0.hh
  undefine-3.0.hh
  entity.hh
  entity.cc
  entityseed.hh
  hierarchiciterator.hh
  algebra.hh
  geometry.hh
  geometry.cc
  projection.hh
  transformation.hh
  leveliterator.hh
  leafiterator.hh
  treeiterator.hh
  intersection.hh
  intersection.cc
  intersectioniterator.hh
  capabilities.hh
  gridfactory.hh
  dgfparser.hh
  albertareader.hh
  gridfamily.hh
  gridview.hh
  persistentcontainer.hh
  backuprestore.hh
  geometryreference.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/albertagrid/)
