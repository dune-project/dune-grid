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

    // In 2D, also test refinement
    for (int refineOpt = 0; refineOpt <= 1; ++refineOpt) {
      for (int refines = 0; refines <= 1; ++refines) {
        check_yasp(YaspFactory<2,Dune::EquidistantCoordinates<double,2> >::buildGrid(refineOpt == 1, refines));
        check_yasp(YaspFactory<2,Dune::EquidistantOffsetCoordinates<double,2> >::buildGrid(refineOpt == 1, refines));
        check_yasp(YaspFactory<2,Dune::TensorProductCoordinates<double,2> >::buildGrid(refineOpt == 1, refines));
      }
    }

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
