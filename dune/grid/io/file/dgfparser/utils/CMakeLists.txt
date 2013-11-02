set(GRIDTYPE ONEDGRID)
set(GRIDDIM 1)

set(EXTRA_PROGRAMS dgf2dgf gmsh2dgf)

foreach(exe ${EXTRA_PROGRAMS})
  add_dgf_executable(${exe}  EXCLUDE_FROM_ALL ${exe}.cc)
  target_link_libraries(${exe} dunegrid ${DUNE_LIBS})
endforeach(exe ${EXTRA_PROGRAMS})

if(ALUGRID_FOUND)
  add_dgf_executable(gmsh2dgfalusimplex EXCLUDE_FROM_ALL gmsh2dgf.cc GRIDDIM 2 WORLDDIM 2
    GRIDTYPE ALUGRID_SIMPLEX)
  add_dune_alugrid_flags(gmsh2dgfalusimplex)
  set(EXTRA_PROGRAMS ${EXTRA_PROGRAMS} gmsh2dgfalusimplex)
endif(ALUGRID_FOUND)
