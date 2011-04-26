// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id: test-ug.cc 3316 2006-12-01 13:45:43Z carsten $

#include <config.h>

#if UG_LGMDOMAIN

#include <iostream>

/*

   Instantiate UG-Grid and feed it to the generic gridcheck()

 */

#include <dune/grid/uggrid.hh>

#include "gridcheck.cc"
#include "checkgeometryinfather.cc"
#include "checkintersectionit.cc"

// Test parallel interface if a parallel UG is used
#ifdef ModelP
#include <mpi.h>
#include "checkcommunicate.cc"
#endif

int main (int argc , char **argv) try
{
#ifdef ModelP
  // initialize MPI
  MPI_Init(&argc,&argv);
#endif

  // ////////////////////////////////////////////////////////////////////////
  //  Do the standard grid test for a 2d UGGrid
  // ////////////////////////////////////////////////////////////////////////
  // extra-environment to check destruction
  {

    std::cout << "Testing UGGrid<3> with LGM domain" << std::endl;

    Dune::UGGrid<3> grid3d;

    grid3d.createLGMGrid("couplex2.lgm");

    // check macro grid
    gridcheck(grid3d);

    // create hybrid grid
    grid3d.mark(1, grid3d.leafbegin<0>());
    grid3d.adapt();

    gridcheck(grid3d);

    grid3d.globalRefine(1);
    gridcheck(grid3d);

    // check the method geometryInFather()
    checkGeometryInFather(grid3d);

    // check the intersection iterator
    checkIntersectionIterator(grid3d);

#ifdef ModelP
    // check communication interface
    checkCommunication(grid3d,-1,Dune::dvverb);
    for(int l=0; l<=grid3d.maxLevel(); ++l)
      checkCommunication(grid3d,l,Dune::dvverb);
#endif
  }

#ifdef ModelP
  // Terminate MPI
  MPI_Finalize();
#endif

  return 0;
}
catch (Dune::Exception& e) {
  std::cerr << e << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}

#else  // UG_LGMDOMAIN
int main () {
  return 77;
}
#endif // UG_LGMDOMAIN
