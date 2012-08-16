set(HEADERS
  geldesc.hh
  ghmesh.cc
  ghmesh.hh
  grapecommon.cc
  grapecommon.hh
  grapedatadisplay.cc
  grapegriddisplay.cc
  grapehmesh.hh
  grapeinclude.hh
  grapewrapper.hh
  partitiondisplay.cc
  combinedgrapedisplay.cc
  grape_undefs.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/io/visual/grape)
