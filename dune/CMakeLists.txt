# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_subdirectory(grid)
# if Python bindings are enabled, include necessary sub directories.
if( DUNE_ENABLE_PYTHONBINDINGS )
  add_subdirectory("python")
endif()
