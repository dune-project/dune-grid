# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

set(HEADERS
  identitygridentity.hh
  identitygridentityseed.hh
  identitygridhierarchiciterator.hh
  identitygridgeometry.hh
  identitygridintersectioniterator.hh
  identitygridintersections.hh
  identitygridindexsets.hh
  identitygridleafiterator.hh
  identitygridleveliterator.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/identitygrid)
