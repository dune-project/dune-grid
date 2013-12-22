exclude_from_headercheck(main-volume.hh)

if(ALUGRID_FOUND AND ENABLE_GRID_EXTRA_UTILS)
add_executable(dune-alusimplex-3d-to-vtk-volume alusimplex-3d-volume.cc)
add_executable(dune-alucube-3d-to-vtk-volume alucube-3d-volume.cc)
add_dune_alugrid_flags(dune-alusimplex-3d-to-vtk-volume)
add_dune_alugrid_flags(dune-alucube-3d-to-vtk-volume)
install(TARGETS dune-alusimplex-3d-to-vtk-volume dune-alucube-3d-to-vtk-volume
  DESTINATION ${CMAKE_INSTALL_BINDIR})
message(AUTHOR_WARNING "TODO: Check whether dune-gmsh-to-alusimplex-3d dune-gmsh-to-alucube-3dis handled correctly")
endif(ALUGRID_FOUND AND ENABLE_GRID_EXTRA_UTILS)
