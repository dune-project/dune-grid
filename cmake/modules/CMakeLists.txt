set(modules
  DuneGridMacros.cmake
  FindAlberta.cmake
  FindALUGrid.cmake
  FindAmiraMesh.cmake
  FindGrape.cmake
  FindPsurface.cmake
  GridType.cmake
  UseUG.cmake)
include(GNUInstallDirs)
install(FILES ${modules} DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/modules)
