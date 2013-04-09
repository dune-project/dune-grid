message(AUTHOR_WARNING "TODO: Make sure that header man.hh is in the tarball")

#set(HEADERCHECK_IGNORE man.hh)

if(ALUGRID_FOUND)
add_executable(dune-gmsh-to-alusimplex-3d alusimplex-3d.cc)
add_executable(dune-gmsh-to-alucube-3d alucube-3d.cc)
add_dune_alugridflags(dune-gmsh-to-alusimplex-3d)
add_dune_alugridflags(dune-gmsh-to-alucube-3d)
install(PROGRAMS dune-gmsh-to-alusimplex-3d dune-gmsh-to-alucube-3d
  DESTINATION ${CMAKE_INSTALL_BINDIR})
message(AUTHOR_WARNING "TODO: Check whether dune-gmsh-to-alusimplex-3d dune-gmsh-to-alucube-3dis handled correctly")
endif(ALUGRID_FOUND)
