// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <sstream>
#include <string>

#include <dune/common/mpihelper.hh>

#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include "gridcheck.cc"

#include "checkgeometryinfather.cc"
#include "checkintersectionit.cc"
#include "checkcommunicate.cc"
#include "checktwists.cc"

using namespace Dune;


template< class GridType >
void checkALUSerial ( GridType &grid, int maxLevel = 2 )
{
  // be careful, each global refine create 4 x maxlevel elements
  if( grid.comm().rank() == 0 )
    std::cout << ">>> Checking macro grid..." << std::endl;
  gridcheck( grid );

  for( int i = 0; i < maxLevel; ++i )
  {
    grid.globalRefine( DGFGridInfo< GridType >::refineStepsForHalf() );
    if( grid.comm().rank() == 0 )
      std::cout << ">>> Checking grid refined " << (i+1) << " times..." << std::endl;
    gridcheck( grid );
  }

  // check the method geometryInFather()
  checkGeometryInFather(grid);

#if 1
  // check the intersection iterator and the geometries it returns
  checkIntersectionIterator(grid);
#endif
  checkTwists( grid.leafView(), NoMapTwist() );
}


int main ( int argc, char **argv )
try
{
  // this method calls MPI_Init, if MPI is enabled
  MPIHelper &mpihelper = MPIHelper::instance( argc, argv );
  int myrank = mpihelper.rank();

  // extra-environment to check destruction
  {
    typedef ALUSimplexGrid< 2, 2 > GridType;

    {
      const std::string filename( "simplex-testgrid-2-2.dgf" );

      GridPtr< GridType > gridPtr( filename );
      GridType &grid = *gridPtr;
      grid.loadBalance();

      if (myrank == 0)
        std::cout << "Check serial grid" << std::endl;
      checkALUSerial( grid, 2 );
    }
  };

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
