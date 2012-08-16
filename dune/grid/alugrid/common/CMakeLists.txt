set(HEADERS
  backuprestore.hh
  bndprojection.hh
  checkparallel.hh
  capabilities.hh
  declaration.hh
  defaultindexsets.hh
  geostorage.hh
  intersectioniteratorwrapper.hh
  interfaces.hh
  memory.hh
  objectfactory.hh
  persistentcontainer.hh
  transformation.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/alugrid/common)

