set(GRIDTYPE ONEDGRID)
set(GRIDDIM 1)

set(EXTRA_PROGRAMS dgf2dgf gmsh2dgf)
add_executable(dgf2dgf dgf2dgf.cc)
add_executable(gmsh2dgf gmsh2dgf.cc)


if(ALUGRID_FOUND)
  list(APPEND EXTRA_PROGAMS gmsh2dgfalusimplex)
  add_executable(gmsh2dgfalusimplex gmsh2dgf.cc)
  set_property(TARGET gmsh2dgfalusimplex APPEND_STRING
    COMPILE_DEFINITIONS "-DALUGRID_SIMPLEX -DGRIDDIM=2")
  add_dune_alugrid_flags(gmsh2dgfalusimplex)
endif(ALUGRID_FOUND)

foreach(exe ${EXTRA_PROGRAMS})
    target_link_libraries(${exe} dunegrid ${DUNE_LIBS})
endforeach(exe ${EXTRA_PROGRAMS})

message(AUTHOR_WARNING "TODO: Build ${EXTRA_PROGRAMS} only on demand")
