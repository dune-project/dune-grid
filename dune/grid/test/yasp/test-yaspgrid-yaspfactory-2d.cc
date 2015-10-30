// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/tensorgridfactory.hh>

#include "test-yaspgrid.hh"

int main (int argc , char **argv) {
  try {
    // Initialize MPI, if present
    Dune::MPIHelper::instance(argc, argv);

    check_yasp(YaspFactory<2,Dune::EquidistantCoordinates<double,2> >::buildGrid(true, 0));
    check_yasp(YaspFactory<2,Dune::EquidistantOffsetCoordinates<double,2> >::buildGrid(true, 0));
    check_yasp(YaspFactory<2,Dune::TensorProductCoordinates<double,2> >::buildGrid(true, 0));

    // In 2D, also test refinement
    for (int refineOpt = 0; refineOpt <= 1; ++refineOpt) {
      check_yasp(YaspFactory<2,Dune::EquidistantCoordinates<double,2> >::buildGrid(refineOpt == 1, 1));
      check_yasp(YaspFactory<2,Dune::EquidistantOffsetCoordinates<double,2> >::buildGrid(refineOpt == 1, 1));
      check_yasp(YaspFactory<2,Dune::TensorProductCoordinates<double,2> >::buildGrid(refineOpt == 1, 1));
    }

    // And periodicity
//    check_yasp(YaspFactory<2,Dune::EquidistantCoordinates<double,2> >::buildGrid(true, 0, true));
//    check_yasp(YaspFactory<2,Dune::EquidistantOffsetCoordinates<double,2> >::buildGrid(true, 0, true));
//    check_yasp(YaspFactory<2,Dune::TensorProductCoordinates<double,2> >::buildGrid(true, 0, true));


  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
