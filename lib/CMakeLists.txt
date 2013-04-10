message(AUTHOR_WARNING "TODO: Check whether dunealbertagrid library is correctly built from the albertagridlibs")

if(UG_FOUND)
  set(UGLIB _DUNE_TARGET_OBJECTS:uggrid_)
endif(UG_FOUND)

dune_add_library(dunegrid _DUNE_TARGET_OBJECTS:onedgrid_ ${UGLIB}
  _DUNE_TARGET_OBJECTS:dgfparser_  _DUNE_TARGET_OBJECTS:dgfparserblocks_ ADD_LIBS ${DUNE_LIBS})
add_dune_ug_flags(dunegrid)
add_dune_alugrid_flags(dunegrid)

foreach(_dim ${ALBERTA_WORLD_DIMS})
  dune_add_library(dunealbertagrid_${_dim}d _DUNE_TARGET_OBJECTS:albertagrid_${_dim}d_
    _DUNE_TARGET_OBJECTS:dgfparser_  _DUNE_TARGET_OBJECTS:dgfparserblocks_ ADD_LIBS ${DUNE_LIBS})
  add_dune_alberta_flags(dunealbertagrid_${_dim}d GRIDDIM ${_dim})
  list(APPEND DUNE_ALBERTA_LIBS dunealbertagrid_${dim}d)
endforeach(_dim "${ALBERTA_DIMS}")
