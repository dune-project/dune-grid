add_subdirectory(test)
add_subdirectory(dgfparser)
add_subdirectory(vtk)
add_subdirectory(gnuplot)

set(HEADERS
  dgfparser.hh
  gmshreader.hh
  gmshwriter.hh
  gnuplot.hh
  printgrid.hh
  starcdreader.hh
  vtk.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/io/file)
