set(HEADERS
  b64enc.hh
  basicwriter.hh
  boundaryiterators.hh
  boundarywriter.hh
  common.hh
  corner.hh
  corneriterator.hh
  dataarraywriter.hh
  function.hh
  functionwriter.hh
  pointiterator.hh
  pvtuwriter.hh
  skeletonfunction.hh
  subsamplingvtkwriter.hh
  streams.hh
  volumeiterators.hh
  volumewriter.hh
  vtksequencewriter.hh
  vtksequencewriterbase.hh
  vtkwriter.hh
  vtuwriter.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/io/file/vtk)
