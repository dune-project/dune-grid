# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

# make sure the example works and run in during tests
dune_python_add_test(NAME pyexample
                     SCRIPT example.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS quick)
dune_symlink_to_source_files(FILES example.py test2d_offset.dgf circle1storder.msh)

if(SKBUILD)
  install(FILES
    ../../../../doc/grids/dgf/test2d_offset.dgf ../../../../doc/grids/gmsh/circle1storder.msh
    DESTINATION python/dune/grid/tutorial)
endif()
