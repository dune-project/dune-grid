include(GNUInstallDirs)

install(FILES
  AddAlbertaFlags.cmake
  DuneGridMacros.cmake
  FindAlberta.cmake
  GridType.cmake
  UseUG.cmake
  DESTINATION ${DUNE_INSTALL_MODULEDIR})
