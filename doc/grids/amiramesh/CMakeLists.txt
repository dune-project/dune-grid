set(AMIRA_GRIDS
  cube-testgrid-3d.am
  hybrid-testgrid-3d.am
  simplex-testgrid-2d.am
  hybrid-testgrid-2d.am
  simplex-testgrid-3d.am)

install(FILES ${AMIRA_GRIDS}
  DESTINATION ${CMAKE_INSTALL_DOCDIR}/grids/amiramesh)
