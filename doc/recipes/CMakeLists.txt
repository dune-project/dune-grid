# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception


dune_add_test(SOURCES recipe-iterate-over-grid.cc LINK_LIBRARIES Dune::Grid)
dune_add_test(SOURCES recipe-integration.cc LINK_LIBRARIES Dune::Grid)
