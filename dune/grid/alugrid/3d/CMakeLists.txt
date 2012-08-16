if(ALUGRID_FOUND)
  set(SOURCES mappings.cc topology.cc faceutility.cc)
else(ALUGRID_FOUND)
  set(SOURCES topology.cc)
endif(ALUGRID_FOUND)
add_library(alu3d OBJECT ${SOURCES})
set_property(SOURCE ${SOURCES} APPEND_STRING
  PROPERTY COMPILE_FLAGS ${ALUGRID_COMPILE_FLAGS})

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
