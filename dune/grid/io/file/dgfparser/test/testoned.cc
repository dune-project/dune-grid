// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>
#include <iostream>
#include <string>

#include "../dgfoned.hh"
// use grid check
#include <dune/grid/test/gridcheck.cc>

using namespace Dune;

int main(int argc, char ** argv, char ** envp)
try
{
  std::cout << std::endl << "start OneDGrid test" << std::endl;

  // this method calls MPI_Init, if MPI is enabled
  MPIHelper::instance(argc,argv);

  {
    typedef OneDGrid GridType;
    std::string filename(SRCDIR "/examplegrid5.dgf");
    GridPtr<GridType> gridptr(filename);

    // run grid check to check grid
    gridcheck(*gridptr);
  }
  {
    typedef OneDGrid GridType;
    std::string filename(SRCDIR "/examplegrid6.dgf");
    GridPtr<GridType> gridptr(filename);

    // run grid check to check grid
    gridcheck(*gridptr);
  }
  return 0;
}
catch (Dune::Exception &e) {
  std::cerr << e << std::endl;
  return 1;
}
catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 1;
}
