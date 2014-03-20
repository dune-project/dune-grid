// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define DISABLE_DEPRECATED_METHOD_CHECK 1

#include <dune/grid/sgrid.hh>
#include <dune/grid/identitygrid.hh>

#include "gridcheck.cc"
#include "checkintersectionit.cc"

using namespace Dune;

// test IdentityGrid for given dimension
template <int dim>
void testDim()
{
  typedef SGrid<dim,dim> GridType;
  int n[dim];
  double h[dim];

  for (int i=0; i<dim; ++i)
  {
    n[i] = 1;
    h[i] = 1.0;
  }

  GridType grid(n,h);

  grid.globalRefine(1);

  IdentityGrid<GridType> identityGrid(grid);

  gridcheck(identityGrid);
  checkIntersectionIterator(identityGrid);
}

int main (int argc, char *argv[]) try
{
  testDim<1>();
  testDim<2>();
  testDim<3>();
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Exception e)
{
  std::cout << e << std::endl;
  return 1;
}
