dune_add_test(SOURCES scsgmappertest.cc)

if(UG_FOUND)
  dune_add_test(SOURCES mcmgmappertest.cc)
endif(UG_FOUND)
