
include(GridType)

set(DUNE_GRID_EXTRA_UTILS "" CACHE BOOL
  "Enable compilation and installation of extra utilities from the \"src\" subdirectory.")

find_package(ALUGrid)
find_package(Alberta)
include(UseUG)
find_package(Grape)
find_package(psurface)
find_package(AmiraMesh)

set(DEFAULT_DGF_GRIDDIM 1)
set(DEFAULT_DGF_WORLDDIM 1)
set(DEFAULT_DGF_GRIDTYPE ONEDGRID)
set(DGF_GRIDTYPES ONEDGRID ALUGRID_CONFORM ALUGRID_SIMPLEX ALBERTAGRID SGRID GEOGRID UGGRID)

dune_define_gridtype(GRID_CONFIG_H_BOTTOM GRIDTYPE ONEDGRID
  ASSERTION "(GRIDDIM == 1) && (WORLDDIM == 1)"
  DUNETYPE "Dune::OneDGrid"
  HEADERS "dune/grid/onedgrid.hh" "dune/grid/io/file/dgfparser/dgfoned.hh")

dune_define_gridtype(GRID_CONFIG_H_BOTTOM GRIDTYPE SGRID
  DUNETYPE "Dune::SGrid< dimgrid, dimworld >"
  HEADERS "dune/grid/sgrid.hh" "dune/grid/io/file/dgfparser/dgfs.hh")

dune_define_gridtype(GRID_CONFIG_H_BOTTOM GRIDTYPE YASPGRID
  ASSERTION "GRIDDIM == WORLDDIM"
  DUNETYPE "Dune::YaspGrid< dimgrid >"
  HEADERS "dune/grid/yaspgrid.hh" "dune/grid/io/file/dgfparser/dgfyasp.hh")

macro(add_dgf_flags target)
  cmake_parse_arguments(DGF "" "GRIDDIM;WORLDDIM;GRIDTYPE" "" ${ARGN})
  if(NOT DGF_GRIDDIM)
    set(DGF_GRIDDIM ${DEFAULT_DGF_GRIDDIM})
  endif(NOT DGF_GRIDDIM)
  if(NOT DGF_WORLDDIM)
    set(DGF_WORLDDIM ${DEFAULT_DGF_WORLDDIM})
  endif(NOT DGF_WORLDDIM)
  if(NOT DGF_GRIDTYPE)
    set(DGF_GRIDTYPE ${DEFAULT_DGF_GRIDTYPE})
  endif(NOT DGF_GRIDTYPE)

  set(replace_args "GRIDDIM.*" "GRIDDIM=${DGF_GRIDDIM}"
    "WORLDDIM.*" "WORLDDIM=${DGF_WORLDDIM}")
  foreach(grid ${DGF_GRIDTYPES})
    list(APPEND replace_args ${grid} ${DGF_GRIDTYPE})
  endforeach(grid ${DGF_GRIDTYPES})
  message("${target} ${replace_args}")
  replace_properties(TARGET ${target}
    PROPERTY COMPILE_DEFINITIONS
    ${replace_args})
endmacro(add_dgf_flags target)

macro(add_dgf_executable target)
  cmake_parse_arguments(DGF "" "GRIDDIM;WORLDDIM;GRIDTYPE" "" ${ARGN})
  add_executable(${target} ${DGF_UNPARSED_ARGUMENTS})
  add_dgf_flags(${target} ${ARGN})
endmacro(add_dgf_executable target)
