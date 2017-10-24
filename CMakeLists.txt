set(HEADERS
  capabilities.hh
  entity.hh
  factory.hh
  geometry.hh
  hierarchical.hh
  intersection.hh
  localview.hh
  numpy.hh
  object.hh
  vtk.hh
  function.hh
  gridview.hh
  idset.hh
  indexset.hh
  mapper.hh
  range.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/python/grid)
