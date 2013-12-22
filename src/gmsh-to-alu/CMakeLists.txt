exclude_from_headercheck(main.hh)

if(ALUGRID_FOUND AND ENABLE_GRID_EXTRA_UTILS)
add_executable(dune-gmsh-to-alusimplex-3d alusimplex-3d.cc)
add_executable(dune-gmsh-to-alucube-3d alucube-3d.cc)
add_dune_alugrid_flags(dune-gmsh-to-alusimplex-3d)
add_dune_alugrid_flags(dune-gmsh-to-alucube-3d)
install(TARGETS dune-gmsh-to-alusimplex-3d dune-gmsh-to-alucube-3d
  DESTINATION ${CMAKE_INSTALL_BINDIR})
message(AUTHOR_WARNING "TODO: Check whether dune-gmsh-to-alusimplex-3d dune-gmsh-to-alucube-3dis handled correctly")
endif(ALUGRID_FOUND AND ENABLE_GRID_EXTRA_UTILS)
