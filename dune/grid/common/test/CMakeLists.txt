# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

dune_add_test(SOURCES scsgmappertest.cc LINK_LIBRARIES dunegrid)
dune_add_test(SOURCES mcmgmappertest.cc LINK_LIBRARIES dunegrid)
