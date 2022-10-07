# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

set(GRIDTYPE ONEDGRID)
set(GRIDDIM 1)

set(EXTRA_PROGRAMS dgf2dgf gmsh2dgf)

foreach(exe ${EXTRA_PROGRAMS})
  add_dgf_executable(${exe}  EXCLUDE_FROM_ALL ${exe}.cc)
  target_link_libraries(${exe} PUBLIC dunegrid)
endforeach(exe ${EXTRA_PROGRAMS})

if(dune-uggrid_FOUND)
  add_dgf_executable(gmsh2dgfugsimplex
    EXCLUDE_FROM_ALL
    gmsh2dgf.cc
    GRIDDIM 2 WORLDDIM 2
    GRIDTYPE UGGRID)
  target_link_libraries(gmsh2dgfugsimplex PUBLIC dunegrid)
  set(EXTRA_PROGRAMS ${EXTRA_PROGRAMS} gmsh2dgfugsimplex)
endif()
