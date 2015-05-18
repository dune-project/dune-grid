// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>

#include <iostream>

#include <dune/grid/sgrid.hh>

#include "gridcheck.hh"
#include "checkgeometryinfather.hh"
#include "checkintersectionit.hh"
#include "checkpartition.hh"

template<int d, int w>
void runtest()
{
  int n[] = { 5, 5, 5, 5 };
  double h[] = { 1.0, 2.0, 3.0, 4.0 };

  std::cout << std::endl << "SGrid<" << d << "," << w << ">" << std::endl;
  Dune::SGrid<d,w> g(n, h);
  gridcheck(g);

  g.globalRefine(1);
  checkGeometryInFather(g);
  checkIntersectionIterator(g);
  checkPartitionType( g.leafGridView() );
  // check geometry lifetime
  checkGeometryLifetime( g.leafGridView() );

  std::cout << std::endl;
}

int main () {
  try {
    runtest<1,1>();
    runtest<2,2>();
    runtest<3,3>();
    //    runtest<4,4>();
    runtest<1,3>();
    runtest<2,3>();

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
