set(HEADERS
  capabilities.hh
  commops.hh
  entity.hh
  enums.hh
  factory.hh
  geometry.hh
  hierarchical.hh
  intersection.hh
  localview.hh
  numpy.hh
  numpycommdatahandle.hh
  object.hh
  vtk.hh
  function.hh
  gridview.hh
  idset.hh
  indexset.hh
  mapper.hh
  range.hh
  simplegridfunction.hh
  singleton.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/python/grid)
