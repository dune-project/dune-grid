// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>
#include <dune/grid/test/gridcheck.cc>
#include <dune/grid/common/gridpart.hh>
#include "../gridtype.hh"
using namespace Dune;
#if HAVE_GRAPE
#include <dune/grid/io/visual/grapegriddisplay.hh>
void test(GridType& grid) {
  gridcheck(grid);
  // HierarchicGridPart<GridType> part(grid);
  // GrapeGridDisplay<GridType> disp(part);
  // disp.display();
}
#else
void test(GridType& grid) {
  gridcheck(grid);
}
#endif

int main(int argc, char ** argv, char ** envp)
try {
  // this method calls MPI_Init, if MPI is enabled
  MPIHelper & mpiHelper = MPIHelper::instance(argc,argv);
  int myrank = mpiHelper.rank();

  if (argc<2) {
    std::cerr << "supply grid file as parameter!" << std::endl;
    return 1;
  }

  std::cout << "tester: start grid reading" << std::endl;

  // create Grid from DGF parser
  GridPtr<GridType> grid(argv[1], mpiHelper.getCommunicator() );

  // display
  // if (myrank <= 0) test(*grid);
  // refine
  std::cout << "tester: refine grid" << std::endl;
  grid->globalRefine(DGFGridInfo<GridType>::refineStepsForHalf());
  // display
  if (myrank <= 0) test(*grid);

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
