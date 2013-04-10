set(modules
  DuneGridMacros.cmake
  FindAlberta.cmake
  FindAmiraMesh.cmake
  FindGrape.cmake
  FindPsurface.cmake
  GridType.cmake
  UseUG.cmake)
include(GNUInstallDirs)
install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
