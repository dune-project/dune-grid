# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

set(GRIDS
  circle1storder.msh
  circle2ndorder.msh
  circle.geo
  curved2d.geo
  curved2d.msh
  hybrid-testgrid-2d.msh
  hybrid-testgrid-3d.msh
  oned-testgrid.msh
  pyramid1storder.msh
  pyramid2ndorder.msh
  pyramid4.msh
  pyramid.geo
  pyramid.msh
  sphere.msh
  telescope1storder.msh
  telescope2ndorder.msh
  telescope.geo
  telescope.msh
  twotets.geo
  twotets.msh
  unitcube.msh)
dune_symlink_to_source_files(FILES ${GRIDS})
install(FILES ${GRIDS} DESTINATION ${CMAKE_INSTALL_DOCDIR}/grids/gmsh)
