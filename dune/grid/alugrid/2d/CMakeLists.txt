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
  get_property(prop SOURCE intersection.cc PROPERTY INCLUDE_DIRECTORIES)
  message("--property INCLUDE_DIRECTORIES=${prop} for intersection.cc")
  add_dune_alugrid_flags(SOURCE_ONLY intersection.cc  alu2dgridfactory.cc)
  get_property(prop SOURCE intersection.cc PROPERTY INCLUDE_DIRECTORIES)
  message("--property INCLUDE_DIRECTORIES=${prop} for intersection.cc")
  add_library(alu2d OBJECT intersection.cc  alu2dgridfactory.cc)
endif(ALUGRID_FOUND)
