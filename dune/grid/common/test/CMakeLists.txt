dune_add_test(SOURCES scsgmappertest.cc)

dune_add_test(SOURCES mcmgmappertest.cc
              CMAKE_GUARD dune-uggrid_FOUND)
