# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_subdirectory(dune)

dune_python_configure_bindings(
  PATH "."
  CMAKE_METADATA_FLAGS DUNE_OPTS_FILE HAVE_ALBERTA HAVE_DUNE_UGGRID
)
if(POLICY CMP0087)
  cmake_policy(SET CMP0087 NEW)
endif()
