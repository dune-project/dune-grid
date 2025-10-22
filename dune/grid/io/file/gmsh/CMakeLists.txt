# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

set(HEADERS
  filereader.hh
  gmsh2parser.hh
  gmsh4reader.hh
  gmsh4reader.impl.hh
  gridcreatorinterface.hh
  types.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/io/file/gmsh)

target_sources(dunegrid PRIVATE
  types.cc)

add_subdirectory(gridcreators)
add_subdirectory(utility)
