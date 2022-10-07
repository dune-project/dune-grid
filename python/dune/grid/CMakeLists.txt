# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_subdirectory(tutorial)
add_python_targets(grid
  __init__
  __main__
  _grids
  grid_generator
  core
  map
  datahandle
)
dune_add_pybind11_module(NAME _grid)
set_property(TARGET _grid PROPERTY LINK_LIBRARIES dunecommon dunegeometry dunegrid APPEND)
if(SKBUILD)
  install(TARGETS _grid LIBRARY DESTINATION python/dune/grid)
endif()
