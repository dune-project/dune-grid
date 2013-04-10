set(GRIDTYPE ONEDGRID)
set(GRIDDIM 1)

set(EXTRA_PROGRAMS dgf2dgf gmsh2dgf)
add_dgf_executable(dgf2dgf dgf2dgf.cc)
add_dgf_executable(gmsh2dgf gmsh2dgf.cc)
foreach(exe ${EXTRA_PROGRAMS})
  target_link_libraries(${exe} dunegrid ${DUNE_LIBS})
endforeach(exe ${EXTRA_PROGRAMS})

message(AUTHOR_WARNING "TODO: Build ${EXTRA_PROGRAMS} only on demand")
