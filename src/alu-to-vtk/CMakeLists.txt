message(AUTHOR_WARNING "TODO: Make sure that header man-volume.hh is in the tarball")

#set(HEADERCHECK_IGNORE man.hh)

if(ALUGRID_FOUND)
add_executable(dune-alusimplex-3d-to-vtk-volume alusimplex-3d-volume.cc)
add_executable(dune-alucube-3d-to-vtk-volume alucube-3d-volume.cc)
add_dune_alugridflags(dune-alusimplex-3d-to-vtk-volume)
add_dune_alugridflags(dune-alucube-3d-to-vtk-volume)
install(PROGRAMS dune-alusimplex-3d-to-vtk-volume dune-alucube-3d-to-vtk-volume
  DESTINATION ${CMAKE_INSTALL_BINDIR})
message(AUTHOR_WARNING "TODO: Check whether dune-gmsh-to-alusimplex-3d dune-gmsh-to-alucube-3dis handled correctly")
endif(ALUGRID_FOUND)
