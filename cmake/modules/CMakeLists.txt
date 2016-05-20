include(GNUInstallDirs)

install(FILES
  AddAlbertaFlags.cmake
  AddAmiraMeshFlags.cmake
  AddPsurfaceFlags.cmake
  CheckExperimentalGridExtensions.cmake
  DuneGridMacros.cmake
  FindAlberta.cmake
  FindAmiraMesh.cmake
  FindPsurface.cmake
  GridType.cmake
  UseUG.cmake
  DESTINATION ${DUNE_INSTALL_MODULEDIR})
