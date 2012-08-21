message(AUTHOR_WARNING "TODO: Check whether dunealbertagrid library is correctly built ftomm the albertagridlibs")

if(ALUGRID_FOUND)
  set(ALULIBS $<TARGET_OBJECTS:alu2d> $<TARGET_OBJECTS:alu3d>)
else(ALUGRID_FOUND)
  set(ALULIBS $<TARGET_OBJECTS:alu3d>)
endif(ALUGRID_FOUND)

if(UG_FOUND)
  set(UGLIB $<TARGET_OBJECTS:uggrid>)
endif(UG_FOUND)

dune_add_library(dunegrid $<TARGET_OBJECTS:onedgrid> ${UGLIB} ${ALULIBS}
  $<TARGET_OBJECTS:dgfparser>  $<TARGET_OBJECTS:dgfparserblocks> ADD_LIBS ${DUNE_LIBS})

foreach(_dim ${ALBERTA_WORLD_DIMS})
  dune_add_library(dunealbertagrid_${_dim}d $<TARGET_OBJECTS:albertagrid_${_dim}d>
    $<TARGET_OBJECTS:onedgrid> ${UGLIB} ${ALULIBS}
  $<TARGET_OBJECTS:dgfparser>  $<TARGET_OBJECTS:dgfparserblocks> ADD_LIBS ${DUNE_LIBS})
  list(APPEND DUNE_ALBERTA_LIBS dunealbertagrid_${dim}d)
endforeach(_dim "${ALBERTA_DIMS}")
