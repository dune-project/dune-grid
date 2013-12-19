set(GRIDS
  star.vrt
  tets.vrt
  withprism.vrt
  withpyramid.vrt
  star.cel
  tets.cel
  withprism.cel
  withpyramid.cel)

install(FILES ${GRIDS}
  DESTINATION ${CMAKE_INSTALL_DOCDIR}/grids/starcd)
