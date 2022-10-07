# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

set(HEADERS
  capabilities.hh
  commops.hh
  entity.hh
  enums.hh
  factory.hh
  function.hh
  geometry.hh
  gridview.hh
  hierarchical.hh
  idset.hh
  indexset.hh
  intersection.hh
  localview.hh
  mapper.hh
  numpy.hh
  numpycommdatahandle.hh
  pygridfunction.hh
  object.hh
  range.hh
  simplegridfunction.hh
  vtk.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/python/grid)
