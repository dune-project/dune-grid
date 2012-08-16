set(GRIDDIM 2)

if(ALBERTA_FOUND)
  set(ALBERTA_TESTS
    test-alberta3d-refine)
  add_executable(test-alberta3d-refine test-alberta3d-refine.cc)
  add_dune_alberta_flags(test-alberta3d-refine)
  add_test(test-alberta3d-refine)
endif(ALBERTA_FOUND)
