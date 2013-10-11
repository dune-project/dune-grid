// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include <config.h>

#include <iostream>

#include <dune/grid/yaspgrid.hh>

#include "gridcheck.cc"
#include "checkcommunicate.cc"
#include "checkgeometryinfather.cc"
#include "checkintersectionit.cc"
#include "checkadaptation.cc"
#include "checkpartition.cc"

int rank;

template <int dim>
void check_yasp(bool p0=false) {
  typedef Dune::FieldVector<double,dim> fTupel;

  std::cout << std::endl << "YaspGrid<" << dim << ">";
  if (p0) std::cout << " periodic\n";
  std::cout << std::endl << std::endl;

  fTupel Len; Len = 1.0;
  Dune::array<int,dim> s;
  std::fill(s.begin(), s.end(), 2);
  s[0] = 6;
  std::bitset<dim> p(0);
  p[0] = p0;
  int overlap = 1;

#if HAVE_MPI
  Dune::YaspGrid<dim> grid(MPI_COMM_WORLD,Len,s,p,overlap);
#else
  Dune::YaspGrid<dim> grid(Len,s,p,overlap);
#endif

  gridcheck(grid);

  grid.globalRefine(2);

  gridcheck(grid);

  // check communication interface
  checkCommunication(grid,-1,Dune::dvverb);
  for(int l=0; l<=grid.maxLevel(); ++l)
    checkCommunication(grid,l,Dune::dvverb);

  // check geometry lifetime
  checkGeometryLifetime( grid.leafView() );
  // check the method geometryInFather()
  checkGeometryInFather(grid);
  // check the intersection iterator and the geometries it returns
  checkIntersectionIterator(grid);
  // check grid adaptation interface
  checkAdaptRefinement(grid);
  checkPartitionType( grid.leafView() );

  // test operator<<
  std::cout << grid << std::endl;
}

int main (int argc , char **argv) {
  try {
#if HAVE_MPI
    // initialize MPI
    MPI_Init(&argc,&argv);

    // get own rank
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif

    check_yasp<1>();
    //check_yasp<1>(true);
    check_yasp<2>();
    //check_yasp<2>(true);
    check_yasp<3>();
    //check_yasp<3>(true);
    //check_yasp<4>();

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

#if HAVE_MPI
  // Terminate MPI
  MPI_Finalize();
#endif

  return 0;
}
