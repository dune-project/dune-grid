add_subdirectory(common)
add_subdirectory(alugrid)
add_subdirectory(onedgrid)
add_subdirectory(sgrid)
add_subdirectory(uggrid)
add_subdirectory(utility)
add_subdirectory(yaspgrid)
add_subdirectory(io)
add_subdirectory(albertagrid)
add_subdirectory(geometrygrid)
add_subdirectory(test)

if(UG_FOUND)
  set(UGLIB $<TARGET_OBJECTS:uggrid>)
endif(UG_FOUND)



#add_library(grid OBJECT $<TARGET_OBJECTS:onedgrid> ${UGLIB} ${ALULIBS}
# $<TARGET_OBJECTS:dgfparser>)
#  ADD_LIBS ${AMIRAMESHLIB} ${DUNE_LIBS})

set(HEADERS
  albertagrid.hh
  alugrid.hh
  sgrid.hh
  uggrid.hh
  yaspgrid.hh
  onedgrid.hh
  geometrygrid.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid)
