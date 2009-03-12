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

#include "basicunitcube.hh"

using namespace Dune;


template< class GridType >
void checkALUSerial ( GridType &grid, int maxLevel = 2 )
{
  // be careful, each global refine create 8 x maxlevel elements
  if( grid.comm().rank() == 0 )
    std::cout << ">>> Checking macro grid..." << std::endl;
  gridcheck( grid );

  for( int i = 0; i < maxLevel; ++i )
  {
    grid.globalRefine( DGFGridInfo<GridType> :: refineStepsForHalf() );
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
}


int main ( int argc, char **argv )
try
{
  // this method calls MPI_Init, if MPI is enabled
  MPIHelper &mpihelper = MPIHelper::instance( argc, argv );
  int myrank = mpihelper.rank();
  int mysize = mpihelper.size();

  // extra-environment to check destruction
  {
    typedef ALUCubeGrid< 3, 3 > GridType;

    {
      if( myrank == 0 )
        std::cout << ">>> Checking empty grid..." << std::endl;
      GridType grid;
      checkALUSerial( grid );
    }

    {
      if( myrank == 0 )
        std::cerr << ">>> Checking twisted unit cube..." << std::endl;
      GridFactory< GridType > factory;
      BasicUnitCube< GridType::dimension >::insertVertices( factory );
      BasicUnitCube< GridType::dimension >::insertCubes( factory );
      GridType *grid = factory.createGrid();
      checkALUSerial( *grid );
      delete grid;
    }

    {
      std::string filename;
      if( argc > 1 )
        filename = argv[ 1 ];
      else if( mysize <= 2 )
        filename = "cube-testgrid-3-3.dgf";
      else
        filename = "cube-testgrid-3-3-large.dgf";

      GridPtr< GridType > gridPtr( filename );
      GridType &grid = *gridPtr;
      grid.loadBalance();

      if( myrank == 0 )
        std::cout << ">>> Checking serial grid..." << std::endl;
      checkALUSerial( grid, 1 );
    }
  }

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
