# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

install(FILES
    entity.hh
    entityiterator.hh
    entityset.hh
    geometry.hh
    grid.hh
    gridview.hh
    indexset.hh
    intersection.hh
    intersectioniterator.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/concepts)
