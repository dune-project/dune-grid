// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>
#include <iostream>
#include <string>

#include "../dgfalberta.hh"
// use grid check
#include <dune/grid/test/gridcheck.cc>

using namespace Dune;

int main(int argc, char ** argv, char ** envp)
{
  // this method calls MPI_Init, if MPI is enabled
  MPIHelper::instance(argc,argv);

  {
    typedef AlbertaGrid<GRIDDIM,GRIDDIM> GridType;
    std::string filename;
    if(GRIDDIM == 2) filename += "examplegrid5.dgf";
    if(GRIDDIM == 3) filename += "examplegrid7.dgf";
    GridPtr<GridType> gridptr(filename);

    // run grid check to check grid
    gridcheck(*gridptr);
  }
  return 0;
}
