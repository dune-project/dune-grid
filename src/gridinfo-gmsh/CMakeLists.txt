# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

if(ALBERTA_2D_LIB_FOUND AND ENABLE_GRID_EXTRA_UTILS)
  add_executable(dune-gridinfo-gmsh-alberta-2d gridinfo-gmsh-alberta-2d.cc)
  target_link_libraries(dune-gridinfo-gmsh-alberta-2d PRIVATE dunegrid)
  add_dune_alberta_flags(dune-gridinfo-gmsh-alberta-2d GRIDDIM 2)
  set(PROGRAMS ${PROGRAMS} dune-gridinfo-gmsh-alberta-2d)
endif(ALBERTA_2D_LIB_FOUND AND ENABLE_GRID_EXTRA_UTILS)

if(ALBERTA_3D_LIB_FOUND AND ENABLE_GRID_EXTRA_UTILS)
  add_executable(dune-gridinfo-gmsh-alberta-3d gridinfo-gmsh-alberta-3d.cc)
  target_link_libraries(dune-gridinfo-gmsh-alberta-3d PRIVATE dunegrid)
  add_dune_alberta_flags(dune-gridinfo-gmsh-alberta-3d GRIDDIM 3)
  set(PROGRAMS ${PROGRAMS} dune-gridinfo-gmsh-alberta-3d)
endif(ALBERTA_3D_LIB_FOUND AND ENABLE_GRID_EXTRA_UTILS)

if(dune-uggrid_FOUND AND ENABLE_GRID_EXTRA_UTILS)
  add_executable(dune-gridinfo-gmsh-ug-2d gridinfo-gmsh-ug-2d.cc)
  add_dune_ug_flags(dune-gridinfo-gmsh-ug-2d)
  add_executable(dune-gridinfo-gmsh-ug-3d gridinfo-gmsh-ug-3d.cc)
  add_dune_ug_flags(dune-gridinfo-gmsh-ug-3d)
  set(PROGRAMS ${PROGRAMS} dune-gridinfo-gmsh-ug-2d dune-gridinfo-gmsh-ug-3d)
endif()

if(PROGRAMS)
  install(TARGETS ${PROGRAMS} DESTINATION ${CMAKE_INSTALL_BINDIR})
endif(PROGRAMS)
