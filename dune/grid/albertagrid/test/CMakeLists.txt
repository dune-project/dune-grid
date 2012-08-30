set(GRIDDIM 2)

if(ALBERTA_FOUND)
  set(ALBERTA_TESTS
    test-alberta3d-refine)
  add_executable(test-alberta3d-refine test-alberta3d-refine.cc)
  add_dune_alberta_flags(test-alberta3d-refine GRIDDIM 3)
  add_test(test-alberta3d-refine test-alberta3d-refine)
endif(ALBERTA_FOUND)

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} "${ALBERTA_TESTS}")
