add_subdirectory(albertagrid)
add_subdirectory(common)
add_subdirectory(geometrygrid)
add_subdirectory(identitygrid)
add_subdirectory(io)
add_subdirectory(onedgrid)
add_subdirectory(test)
add_subdirectory(uggrid)
add_subdirectory(utility)
add_subdirectory(yaspgrid)

if(dune-uggrid_FOUND)
  set(UGLIB $<TARGET_OBJECTS:uggrid>)
endif()

install(FILES
  albertagrid.hh
  geometrygrid.hh
  identitygrid.hh
  onedgrid.hh
  uggrid.hh
  yaspgrid.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid)
