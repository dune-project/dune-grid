if(ALUGRID_FOUND)
  set(SOURCES mappings.cc topology.cc faceutility.cc)
else(ALUGRID_FOUND)
  set(SOURCES topology.cc)
endif(ALUGRID_FOUND)
add_dune_alugrid_flags(SOURCE_ONLY ${SOURCES})
add_library(alu3d OBJECT ${SOURCES})

set(HEADERS
  alu3dinclude.hh
  alu3dgridfactory.hh
  alu3dgridfactory.cc
  datahandle.hh
  grid.hh
  grid_imp.cc
  grid_inline.hh
  mappings.hh
  mappings_imp.cc
  topology.hh
  entity.hh
  entity_imp.cc
  entity_inline.hh
  entityseed.hh
  faceutility.hh
  faceutility_imp.cc
  geometry.hh
  geometry_imp.cc
  indexsets.hh
  iterator.hh
  iterator.cc
  iterator_imp.cc
  alu3diterators.hh
  lbdatahandle.hh
  alugrid.hh
  capabilities.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/alugrid/3d)
