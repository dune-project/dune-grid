add_subdirectory(amiramesh)
add_subdirectory(test EXCLUDE_FROM_ALL)
add_subdirectory(dgfparser)
add_subdirectory(vtk)
add_subdirectory(gnuplot)

set(HEADERS
  amirameshreader.hh
  amirameshwriter.hh
  dgfparser.hh
  gmshreader.hh
  gmshwriter.hh
  gnuplot.hh
  starcdreader.hh
  vtk.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/io/file)
