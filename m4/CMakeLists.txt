set(ALLM4S
  alberta.m4
  alugrid.m4
  amiramesh.m4
  dune_experimental_grid_extensions.m4
  dune_extra_utils.m4
  dune_grid.m4
  dune_griddim.m4
  dune_gridtype.m4
  grape.m4
  psurface.m4
  ug.m4)

install(FILES ${ALLM4S}
  DESTINATION ${CMAKE_INSTALL_DATADIR}/aclocal)
