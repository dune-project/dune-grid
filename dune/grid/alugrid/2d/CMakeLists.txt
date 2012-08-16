set(HEADERS
  alu2dinclude.hh
  entity_imp.cc
  geometry_imp.cc
  grid_imp.cc
  intersection.hh
  iterator.hh
  entity.hh
  geometry.hh
  grid.hh
  indexsets.hh
  intersection_imp.cc
  iterator_imp.cc
  datahandle.hh
  alu2dgridfactory.hh
  alugrid.hh
  capabilities.hh
  bndprojection.hh
  intersection.cc
  alu2dgridfactory.cc
  entityseed.hh)

install(FILES ${HEADERS}
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/alugrid/2d)

if(ALUGRID_FOUND)
  add_library(alu2d OBJECT intersection.cc  alu2dgridfactory.cc)
  #target_link_libraries(alu2d ${ALUGRID_LIBS} ${DUNE_LIBS})
  set_property(TARGET alu2d APPEND_STRING
    COMPILE_FLAGS "${ALUGRID_COMPILE_FLAGS}")
endif(ALUGRID_FOUND)
