set(GRIDTYPE ONEDGRID)
set(GRIDDIM 1)

set(EXTRA_PROGRAMS dgf2dgf gmsh2dgf)

foreach(exe ${EXTRA_PROGRAMS})
  add_dgf_executable(${exe}  EXCLUDE_FROM_ALL ${exe}.cc)
  target_link_libraries(${exe} PUBLIC dunegrid ${DUNE_LIBS})
endforeach(exe ${EXTRA_PROGRAMS})

if(dune-uggrid_FOUND)
  add_dgf_executable(gmsh2dgfugsimplex
    EXCLUDE_FROM_ALL
    gmsh2dgf.cc
    GRIDDIM 2 WORLDDIM 2
    GRIDTYPE UGGRID)
  add_dune_ug_flags(gmsh2dgfugsimplex)
  set(EXTRA_PROGRAMS ${EXTRA_PROGRAMS} gmsh2dgfugsimplex)
endif()
