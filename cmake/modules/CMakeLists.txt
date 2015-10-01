set(modules
  AddAlbertaFlags.cmake
  AddALUGridFlags.cmake
  AddAmiraMeshFlags.cmake
  AddPsurfaceFlags.cmake
  CheckExperimentalGridExtensions.cmake
  DuneGridMacros.cmake
  FindAlberta.cmake
  FindALUGrid.cmake
  FindAmiraMesh.cmake
  FindPsurface.cmake
  GridType.cmake
  UseUG.cmake)
include(GNUInstallDirs)
install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
