set(GRIDS
  grid-1-1.amc
  grid-1-2.amc
  grid-2-2.amc
  grid-1-3.amc
  grid-2-3.amc
  grid-3-3.amc
  grid-2-4.amc
  macro.amc
  periodic-torus.amc)

install(FILES ${GRIDS}
  DESTINATION ${CMAKE_INSTALL_DOCDIR}/grids/amc)
