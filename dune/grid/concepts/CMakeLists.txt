# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

install(FILES
    datahandle.hh
    entity.hh
    entityiterator.hh
    geometry.hh
    grid.hh
    gridview.hh
    indexidset.hh
    intersection.hh
    intersectioniterator.hh
    messagebuffer.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/concepts)

add_subdirectory(archetypes)