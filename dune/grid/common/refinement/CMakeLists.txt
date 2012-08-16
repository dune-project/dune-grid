set(HEADERS
  base.cc
  hcube.cc
  simplex.cc
  hcubetriangulation.cc
  prismtriangulation.cc
  pyramidtriangulation.cc)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/common/refinement)
