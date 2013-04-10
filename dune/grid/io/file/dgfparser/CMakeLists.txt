add_subdirectory(blocks)
add_subdirectory(test EXCLUDE_FROM_ALL)
add_subdirectory(utils)

dune_add_library(dgfparser OBJECT dgfparser.cc dgfug.cc)


set(HEADERS
  dgfalu.cc
  dgfexception.hh
  dgfalu.hh
  dgfug.hh
  dgfparser.hh
  dgfgeogrid.hh
  dgfwriter.hh
  dgfyasp.hh
  entitykey.hh
  dgfs.hh
  entitykey_inline.hh
  dgfoned.hh
  dgfgridfactory.hh
  macrogrid.hh
  gridptr.hh
  parser.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/io/file/dgfparser)
