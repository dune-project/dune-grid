dune_add_test(SOURCES scsgmappertest.cc)

dune_add_test(SOURCES mcmgmappertest.cc
              CMAKE_GUARD UG_FOUND)
