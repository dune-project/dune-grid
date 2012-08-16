
find_package(ALUGrid)
find_package(Alberta)
find_package(UG)
find_package(Grape)
find_package(psurface)

message(AUTHOR_WARNING "Using dumb grid macros with fixed includes. TODO: Implement all tests!")
#add_definitions("-DGRIDDIM=\$(GRIDDIM)" "-DWORLDDIM=\$(WORLDDIM)" "-D\$(GRIDTYPE)")
add_definitions("-DGRIDDIM=1" "-DWORLDDIM=1" "-DONEDGRID")

set(GRIDDIM 1)
set(WORLDDIM 1)
set(GRIDTYPE ONDEDGRID)

include(GridType)
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
