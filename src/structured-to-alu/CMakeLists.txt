exclude_from_headercheck(main.hh)

if(ALUGRID_FOUND AND ENABLE_GRID_EXTRA_UTILS)
add_executable(dune-structured-to-alusimplex-3d alusimplex-3d.cc)
add_dune_alugrid_flags(dune-structured-to-alusimplex-3d)
install(TARGETS dune-structured-to-alusimplex-3d DESTINATION ${CMAKE_INSTALL_BINDIR})
message(AUTHOR_WARNING "TODO: Check whether dune-structured-to-alusimplex-3d is handled correctly")
endif(ALUGRID_FOUND AND ENABLE_GRID_EXTRA_UTILS)
