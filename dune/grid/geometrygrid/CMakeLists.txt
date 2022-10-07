# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

set(HEADERS
  backuprestore.hh
  cachedcoordfunction.hh
  capabilities.hh
  cornerstorage.hh
  coordfunction.hh
  coordfunctioncaller.hh
  datahandle.hh
  declaration.hh
  entity.hh
  entityseed.hh
  geometry.hh
  grid.hh
  gridfamily.hh
  gridview.hh
  hostcorners.hh
  identity.hh
  idset.hh
  indexsets.hh
  intersection.hh
  intersectioniterator.hh
  iterator.hh
  persistentcontainer.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/geometrygrid)
