add_subdirectory(blocks)
add_subdirectory(test)
add_subdirectory(utils)

dune_add_library(dgfparser OBJECT SOURCES dgfparser.cc dgfug.cc)

install(FILES
  dgfexception.hh
  dgfidentitygrid.hh
  dgfug.hh
  dgfparser.hh
  dgfgeogrid.hh
  dgfwriter.hh
  dgfyasp.hh
  entitykey.hh
  entitykey_inline.hh
  dgfoned.hh
  dgfgridfactory.hh
  macrogrid.hh
  gridptr.hh
  parser.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/io/file/dgfparser)
