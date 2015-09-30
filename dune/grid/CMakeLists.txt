add_subdirectory(albertagrid)
add_subdirectory(alugrid)
add_subdirectory(common)
add_subdirectory(geometrygrid)
add_subdirectory(identitygrid)
add_subdirectory(io)
add_subdirectory(onedgrid)
add_subdirectory(test)
add_subdirectory(uggrid)
add_subdirectory(utility)
add_subdirectory(yaspgrid)

if(UG_FOUND)
  set(UGLIB $<TARGET_OBJECTS:uggrid>)
endif(UG_FOUND)

set(HEADERS
  albertagrid.hh
  alugrid.hh
  geometrygrid.hh
  identitygrid.hh
  onedgrid.hh
  sgrid.hh
  uggrid.hh
  yaspgrid.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid)
