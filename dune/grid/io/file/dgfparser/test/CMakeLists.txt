add_definitions(-DDUNE_GRID_EXAMPLE_GRIDS_PATH=${CMAKE_SOURCE_DIR}/doc/grids/\")

set(TESTS testsgrid testyasp testoned)

if(ALUGRID_FOUND)
  list(APPEND TESTS testalu)
endif(ALUGRID_FOUND)

if(ALBERTA_FOUND)
  list(APPEND TESTS testalberta)
endif(ALBERTA_FOUND)

if(UG_FOUND)
  list(APPEND TESTS testug)
endif(UG_FOUND)

foreach(_test ${TESTS})
  add_executable(${_test} main.cc)
  add_test(${_test} ${_test})
endforeach(_test ${TESTS})


if(ALUGRID_FOUND)
  add_dune_alugrid_flags(testalu)
  set_property(TARGET testalu APPEND_STRING
    COMPILE_DEFINITIONS "-DALUGRID_CUBE -DGRIDDIM=3")
endif(ALUGRID_FOUND)

if(ALBERTA_FOUND)
  add_dune_alberta_flags(testalberta GRID_DIM 2)
  set_property(TARGET testalberta APPEND_STRING
    COMPILE_DEFINITIONS "-DALBERTAGRID -DGRIDDIM=2 -DWORLDIM=3")
endif(ALBERTA_FOUND)

if(UG_FOUND)
  add_dune_ug_flags(testug)
  set_property(TARGET testug APPEND_STRING
    COMPILE_DEFINITIONS "-DUGGRID -DGRIDDIM=3")
endif(UG_FOUND)
