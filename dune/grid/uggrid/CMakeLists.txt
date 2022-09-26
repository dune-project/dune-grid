# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

if(dune-uggrid_FOUND)
  target_sources(dunegrid PRIVATE
    uggridintersections.cc
    uggrid.cc
    uggridentity.cc
    boundaryextractor.cc
    boundaryextractor.hh
    uggridindexsets.cc
    uggridfactory.cc
    uggridgeometry.cc
    uggridhieriterator.cc)
endif()

exclude_all_but_from_headercheck(
  boundaryextractor.hh
  ug_undefs.hh)

set(HEADERS
  uggridfactory.hh
  uggridentityseed.hh
  uggridentity.hh
  uggridgeometry.hh
  uggridlocalgeometry.hh
  uggridhieriterator.hh
  uggridleveliterator.hh
  uggridviews.hh
  ugincludes.hh
  uggridintersections.hh
  uggridintersectioniterators.hh
  uggridindexsets.hh
  uggridleafiterator.hh
  uggridrenumberer.hh
  ug_undefs.hh
  uglbgatherscatter.hh
  ugmessagebuffer.hh
  ugwrapper.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/uggrid)

# tricks like undefAllMacros.pl don't have to be shipped, have they?
# Oliver Sander: It depends. undefAllMacros.pl creates the file ug_undefs.hh
# which contains an undef line for each macro encountered in the UG
# headers which DUNE includes.  This should be redone every time UG
# changes, and currently I do it and then check in the new ug_undefs.hh.
# An alternative would be to have the DUNE-build system call undefAllMacros
# for the UG it is configured with.  Then it needed to be shipped.")
