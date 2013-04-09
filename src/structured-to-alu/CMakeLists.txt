message(AUTHOR_WARNING "Make sure that header man.hh is in the tarball")

#set(HEADERCHECK_IGNORE man.hh)

if(ALUGRID_FOUND)
add_executable(dune-structured-to-alusimplex-3d alusimplex-3d.cc)
add_dune_alugridflags(dune-structured-to-alusimplex-3d)
install(PROGRAMS dune-structured-to-alusimplex-3d DESTINATION ${CMAKE_INSTALL_BINDIR})
message(AUTHOR_WARNING "TODO: Check whether dune-structured-to-alusimplex-3d is handled correctly")
endif(ALUGRID_FOUND)
