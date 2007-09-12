// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>
#include <iostream>
#include <string>

#include "../dgfalu.hh"
// use grid check
#include <dune/grid/test/gridcheck.cc>

using namespace Dune;
#if HAVE_GRAPE
#include <dune/grid/io/visual/grapegriddisplay.hh>
template <class GridType>
void test(GridType& grid) {
  //GrapeGridDisplay<GridType> disp(grid);
  //disp.display();
}
#else
template <class GridType>
void test(GridType& grid) {}
#endif

int main(int argc, char ** argv, char ** envp)
try {

  std::cout << std::endl << "start ALUGrid test" << std::endl;
  // this method calls MPI_Init, if MPI is enabled
  MPIHelper::instance(argc,argv);

  {
    typedef ALUSimplexGrid<2,2> GridType;
    std::string filename(SRCDIR "/examplegrid5.dgf");
    GridPtr<GridType> gridptr(filename);

    // run grid check to check grid
    gridcheck(*gridptr);

    test(*gridptr);
  }

  {
    typedef ALUConformGrid<2,2> GridType;
    std::string filename(SRCDIR "/examplegrid5.dgf");
    GridPtr<GridType> gridptr(filename);

    // run grid check to check grid
    gridcheck(*gridptr);

    test(*gridptr);
  }

  {
    typedef ALUCubeGrid<3,3> GridType;
    std::string filename(SRCDIR "/examplegrid6.dgf");
    GridPtr<GridType> gridptr(filename);

    // run grid check to check grid
    gridcheck(*gridptr);

    test(*gridptr);
  }

  {
    typedef ALUSimplexGrid<3,3> GridType;
    std::string filename(SRCDIR "/examplegrid6.dgf");
    GridPtr<GridType> gridptr(filename);

    // run grid check to check grid
    gridcheck(*gridptr);

    test(*gridptr);
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
