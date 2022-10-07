# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

target_sources(dunegrid PRIVATE
    basic.cc  boundarydom.cc  boundaryseg.cc
    cube.cc  dim.cc  general.cc  gridparameter.cc
    interval.cc  periodicfacetrans.cc projection.cc simplex.cc  simplexgeneration.cc  vertex.cc)
#add_dune_mpi_flags(dgfparserblocks)

set(HEADERS
  basic.hh
  boundarydom.hh
  boundaryseg.hh
  cube.hh
  dim.hh
  general.hh
  gridparameter.hh
  interval.hh
  periodicfacetrans.hh
  polygon.hh
  polyhedron.hh
  projection.hh
  simplex.hh
  simplexgeneration.hh
  vertex.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/io/file/dgfparser/blocks)
