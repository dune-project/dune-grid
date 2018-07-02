// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/tensorgridfactory.hh>

#include "test-yaspgrid.hh"

int main (int argc , char **argv) {
  try {
    // Initialize MPI, if present
    const auto & mpiHelper = Dune::MPIHelper::instance(argc, argv);

    std::string testID =
      "yaspfactory-3d-np" + std::to_string(mpiHelper.size());

    check_yasp(testID + "equidistant",
               YaspFactory<3,Dune::EquidistantCoordinates<double,3> >::buildGrid());
    check_yasp(testID + "equidistantoffset",
               YaspFactory<3,Dune::EquidistantOffsetCoordinates<double,3> >::buildGrid());
    check_yasp(testID + "tensor",
               YaspFactory<3,Dune::TensorProductCoordinates<double,3> >::buildGrid());

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
