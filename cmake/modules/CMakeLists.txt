# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

include(GNUInstallDirs)

install(FILES
  AddAlbertaFlags.cmake
  DuneGridMacros.cmake
  FindAlberta.cmake
  GridType.cmake
  UseUG.cmake
  DESTINATION ${DUNE_INSTALL_MODULEDIR})
