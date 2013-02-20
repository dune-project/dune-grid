add_subdirectory(test EXCLUDE_FROM_ALL)
if(ALBERTA_FOUND)
  set(ALBERTAGRID_SOURCES
    dgfparser.cc
    elementinfo.cc
    indexsets.cc
    macrodata.cc
    meshpointer.cc)
  add_library(albertagrid OBJECT ${ALBERTAGRID_SOURCES})
  add_dune_alberta_flags(albertagrid)
endif(ALBERTA_FOUND)

# create dimension dependant alberta libraries
foreach(_dim RANGE 1 9)
  if(ALBERTA_${_dim}D_FOUND)
    list(APPEND ALBERTAGRID_LIBS albertagrid_${_dim}d)
    add_library(albertagrid_${_dim}d OBJECT ${ALBERTAGRID_SOURCES})
    target_link_libraries(albertagrid_${_dim}d
      ${ALBERTA_BASE_LIBS} alberta_${_dim}d ${DUNE_LIBS})
    set_property(TARGET albertagrid_${_dim}d APPEND_STRING
      COMPILE_DEFINITIONS "-DALBERTA_DIM=${_dim}")
  endif(ALBERTA_${_dim}D_FOUND)
endforeach(_dim RANGE 1 9)

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
  entitypointer.hh
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
  backuprestore.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}//dune/grid/albertagrid/)
