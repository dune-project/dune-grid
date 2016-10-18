// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/yaspgrid/backuprestore.hh>

#include "test-yaspgrid.hh"

int main (int argc , char **argv) {
  try {
    // Initialize MPI, if present
    Dune::MPIHelper::instance(argc, argv);

    // check the backup restore facility
    check_backuprestore(YaspFactory<2,Dune::EquidistantCoordinates<double,2> >::buildGrid());

    // Test again with refinement
    check_backuprestore(YaspFactory<2,Dune::EquidistantCoordinates<double,2> >::buildGrid(true, 1));

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
