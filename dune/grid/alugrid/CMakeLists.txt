add_subdirectory(common)
add_subdirectory(2d)
add_subdirectory(3d)
install(FILES
  interfaces.hh
  defaultindexsets.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/alugrid)
