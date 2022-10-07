# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

dune_add_test(SOURCES globalindexsettest.cc
              LINK_LIBRARIES dunegrid)

dune_add_test(SOURCES persistentcontainertest.cc
              LINK_LIBRARIES dunegrid)

dune_add_test(SOURCES structuredgridfactorytest.cc
              LINK_LIBRARIES dunegrid)

dune_add_test(SOURCES tensorgridfactorytest.cc
              LINK_LIBRARIES dunegrid)

dune_add_test(SOURCES vertexordertest.cc
              LINK_LIBRARIES dunegrid)
