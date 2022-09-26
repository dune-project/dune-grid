# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

set(GRIDS
  basicunitcube-2d.dgf
  cube-2.dgf
  cube_grid.dgf
  cube-testgrid-2-2.dgf
  cube-testgrid-2-3.dgf
  distorted-cube-3.dgf
  examplegrid10a.dgf
  examplegrid10b.dgf
  examplegrid10.dgf
  examplegrid11.dgf
  examplegrid12.dgf
  examplegrid1c.dgf
  examplegrid1gen.dgf
  examplegrid1s.dgf
  examplegrid2a.dgf
  examplegrid2b.dgf
  examplegrid2c.dgf
  examplegrid2d.dgf
  examplegrid2e.dgf
  examplegrid5.dgf
  examplegrid6.dgf
  examplegrid7.dgf
  examplegrid9.dgf
  example-projection.dgf
  grid2Y.dgf
  grid3A.dgf
  grid3Y.dgf
  helix.dgf
  octahedron.dgf
  simplex-testgrid-1-1.dgf
  simplex-testgrid-1-2.dgf
  simplex-testgrid-2-2.dgf
  simplex-testgrid-1-3.dgf
  simplex-testgrid-2-3.dgf
  simplex-testgrid-3-3.dgf
  simplex-testgrid-3-3-large.dgf
  test1d.dgf
  test1d-vertex.dgf
  test2d.dgf
  test2d_offset.dgf
  test2ug.dgf
  test3d.dgf
  torus-2.dgf
  torus-3.dgf
  unstr_cube.dgf)

install(FILES ${GRIDS}
  DESTINATION ${CMAKE_INSTALL_DOCDIR}/grids/dgf)
