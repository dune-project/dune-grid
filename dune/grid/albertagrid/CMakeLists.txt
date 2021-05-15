add_subdirectory(test)

if(Alberta_FOUND)
  set(ALBERTAGRID_SOURCES
    dgfparser.cc
    elementinfo.cc
    indexsets.cc
    macrodata.cc
    meshpointer.cc)

  # Create dimension dependent alberta libraries:
  # Note, all libraries and flags are set PRIVATE only to not mix up
  # flags when linked together.
  foreach(_dim ${ALBERTA_WORLD_DIMS})
    add_library(albertagrid${_dim}d OBJECT ${ALBERTAGRID_SOURCES})
    target_link_libraries(albertagrid${_dim}d PUBLIC Alberta::AlbertaGrid_${_dim}d)
    target_compile_definitions(albertagrid${_dim}d PRIVATE ALBERTA_DIM=${_dim} ENABLE_ALBERTA)
  endforeach(_dim)
endif()

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
