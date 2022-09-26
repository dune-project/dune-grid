# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

set(HEADERS
  backuprestore.hh
  coordinates.hh
  partitioning.hh
  structuredyaspgridfactory.hh
  torus.hh
  yaspgridentity.hh
  yaspgridentityseed.hh
  yaspgridgeometry.hh
  yaspgridhierarchiciterator.hh
  yaspgridindexsets.hh
  yaspgridintersection.hh
  yaspgridintersectioniterator.hh
  yaspgrididset.hh
  yaspgridleveliterator.hh
  yaspgridpersistentcontainer.hh
  ygrid.hh)

exclude_all_but_from_headercheck(backuprestore.hh torus.hh coordinates.hh ygrid.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/yaspgrid/)
