# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_subdirectory(blocks)
add_subdirectory(test)
add_subdirectory(utils)

target_sources(dunegrid PRIVATE dgfparser.cc dgfug.cc)

install(FILES
  dgfexception.hh
  dgfidentitygrid.hh
  dgfug.hh
  dgfparser.hh
  dgfgeogrid.hh
  dgfwriter.hh
  dgfyasp.hh
  entitykey.hh
  entitykey_inline.hh
  dgfoned.hh
  dgfgridfactory.hh
  macrogrid.hh
  gridptr.hh
  parser.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/io/file/dgfparser)
