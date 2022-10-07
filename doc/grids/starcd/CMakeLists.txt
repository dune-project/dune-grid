# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

set(GRIDS
  star.vrt
  tets.vrt
  withprism.vrt
  withpyramid.vrt
  star.cel
  tets.cel
  withprism.cel
  withpyramid.cel)

install(FILES ${GRIDS}
  DESTINATION ${CMAKE_INSTALL_DOCDIR}/grids/starcd)
