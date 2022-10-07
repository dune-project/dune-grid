# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

set(HEADERS nulliteratorfactory.hh
  onedgridentity.hh
  onedgridentityseed.hh
  onedgridfactory.hh
  onedgridhieriterator.hh
  onedgridindexsets.hh
  onedgridleafiterator.hh
  onedgridleveliterator.hh
  onedgridlist.hh
  onedgridintersections.hh
  onedgridintersectioniterators.hh
  onedgridviews.hh)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/onedgrid/)

# Some of the OnedGrid headers do not pass make headercheck,
# and they don't need to necessarily. Remove this line when headers are fixed
exclude_from_headercheck( onedgridleafiterator.hh onedgridleveliterator.hh onedgridhieriterator.hh )

target_sources(dunegrid PRIVATE
  onedgrid.cc
  onedgridfactory.cc
  nulliteratorfactory.cc)
