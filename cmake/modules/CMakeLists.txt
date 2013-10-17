set(modules
  AddAlbertaFlags.cmake
  AddALUGridFlags.cmake
  AddAmiraMeshFlags.cmake
  AddGrapeFlags.cmake
  AddPsurfaceFlags.cmake
  CheckExperimentalGridExtensions.cmake
  DuneGridMacros.cmake
  FindAlberta.cmake
  FindALUGrid.cmake
  FindAmiraMesh.cmake
  FindGrape.cmake
  FindPsurface.cmake
  GridType.cmake
  UseUG.cmake)
include(GNUInstallDirs)
install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
