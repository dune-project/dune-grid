if(UG_FOUND)
  set(TESTS
    mcmgmappertest)
endif(UG_FOUND)

set(TESTS ${TESTS} scsgmappertest refinementtest)

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS})

foreach(_t ${TESTS})
  add_executable(${_t} ${_t}.cc)
endforeach(_t ${TESTS})

if(UG_FOUND)
  add_dune_ug_flags(mcmgmappertest)
endif(UG_FOUND)
