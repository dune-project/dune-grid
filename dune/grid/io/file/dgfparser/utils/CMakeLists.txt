message("Processig ${CMAKE_CURRENT_SOURCE_DIR}")
set(GRIDTYPE ONEDGRID)
set(GRIDDIM 1)

set(EXTRA_PROGRAMS dgf2dgf gmsh2dgf)
add_dgf_executable(dgf2dgf dgf2dgf.cc)
add_dgf_executable(gmsh2dgf gmsh2dgf.cc)
foreach(exe ${EXTRA_PROGRAMS})
    target_link_libraries(${exe} dunegrid ${DUNE_LIBS})
endforeach(exe ${EXTRA_PROGRAMS})
message("Adding linker flags")
 target_link_libraries(dgf2dgf dunegrid ${DUNE_LIBS})

if(ALUGRID_FOUND)
  add_dgf_executable(gmsh2dgfalusimplex gmsh2dgf.cc GRIDDIM 2 WORLDDIM 2
    GRIDTYPE ALUGRID_SIMPLEX)
  add_dune_alugrid_flags(gmsh2dgfalusimplex)
endif(ALUGRID_FOUND)


message(AUTHOR_WARNING "TODO: Build ${EXTRA_PROGRAMS} only on demand")
