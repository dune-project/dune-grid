set(GRIDS
  circle1storder.msh
  curved2d.msh
  pyramid4.msh
  telescope2ndorder.msh
  circle2ndorder.msh
  pyramid.geo
  telescope.geo
  circle.geo
  pyramid1storder.msh
  pyramid.msh
  telescope.msh
  curved2d.geo
  pyramid2ndorder.msh
  telescope1storder.msh
  sphere.msh
  oned-testgrid.msh)

install(FILES ${GRIDS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/doc/grids/gmsh)
