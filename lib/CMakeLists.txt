if(dune-uggrid_FOUND)
  set(UGLIB _DUNE_TARGET_OBJECTS:uggrid_)
endif()

if(Alberta_FOUND)
  foreach(_dim ${ALBERTA_WORLD_DIMS})
    list(APPEND ALBERTA_LIBS $<TARGET_OBJECTS:albertagrid${_dim}d>)
  endforeach(_dim)
endif()

# OBJECT is needed since CMake 3.0.0
# but it conflicts with shared libs
set(_OBJECT_FLAG "OBJECT")
if(BUILD_SHARED_LIBS)
  set(_OBJECT_FLAG "")
endif()

dune_add_library(dunegrid
  _DUNE_TARGET_OBJECTS:onedgrid_
  _DUNE_TARGET_OBJECTS:dgfparser_
  _DUNE_TARGET_OBJECTS:dgfparserblocks_
  ${UGLIB} ${ALBERTA_LIBS}
  ADD_LIBS ${DUNE_LIBS})
add_dune_ug_flags(dunegrid ${_OBJECT_FLAG} NO_LINK_DUNEGRID)

install(FILES dunegridam2cmake.lib
  DESTINATION ${CMAKE_INSTALL_BINDIR}/../lib)
