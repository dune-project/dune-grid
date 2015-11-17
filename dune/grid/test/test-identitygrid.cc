// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define DISABLE_DEPRECATED_METHOD_CHECK 1

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/identitygrid.hh>

#include "gridcheck.hh"
#include "checkintersectionit.hh"

using namespace Dune;

// test IdentityGrid for given dimension
template <int dim>
void testDim()
{
  typedef YaspGrid<dim> GridType;
  std::array<int,dim> n;
  std::fill(n.begin(), n.end(), 1 << (5 - dim));
  Dune::FieldVector<double,dim> extension(1.0);

  GridType grid(extension,n);

  grid.globalRefine(1);

  IdentityGrid<GridType> identityGrid(grid);

  gridcheck(identityGrid);
  checkIntersectionIterator(identityGrid);
}

int main (int argc, char *argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);
  testDim<1>();
  testDim<2>();
  testDim<3>();

  return 0;
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Exception e)
{
  std::cout << e << std::endl;
  return 1;
}
