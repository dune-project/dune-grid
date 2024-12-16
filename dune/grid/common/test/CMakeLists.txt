# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

# Link all test targets in this directory against Dune::Grid
link_libraries(Dune::Grid)

dune_add_test(SOURCES scsgmappertest.cc)
dune_add_test(SOURCES mcmgmappertest.cc)
