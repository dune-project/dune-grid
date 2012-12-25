set(GRIDTYPE ONEDGRID)
set(GRIDDIM 1)

set(EXTRA_PROGRAMS dgf2dgf gmsh2dgf)
add_dgf_executable(dgf2dgf dgf2dgf.cc)
add_dgf_executable(gmsh2dgf gmsh2dgf.cc)
foreach(exe ${EXTRA_PROGRAMS})
  target_link_libraries(${exe} dunegrid ${DUNE_LIBS})
endforeach(exe ${EXTRA_PROGRAMS})

if(ALUGRID_FOUND)
  add_dgf_executable(gmsh2dgfalusimplex gmsh2dgf.cc GRIDDIM 2 WORLDDIM 2
    GRIDTYPE ALUGRID_SIMPLEX)
  add_dune_alugrid_flags(gmsh2dgfalusimplex)
  set(EXTRA_PROGRAMS ${EXTRA_PROGRAMS} gmsh2dgfalusimplex)
endif(ALUGRID_FOUND)

# define HAVE_DUNE_GRID for the dgfparser
set_property(TARGET ${EXTRA_PROGRAMS} APPEND PROPERTY
  COMPILE_DEFINITIONS HAVE_DUNE_GRID=1)

message(AUTHOR_WARNING "TODO: Build ${EXTRA_PROGRAMS} only on demand")
