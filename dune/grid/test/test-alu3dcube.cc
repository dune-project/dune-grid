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

template <class GridType>
void makeNonConfGrid(GridType &grid,int level,int adapt) {
  int myrank = grid.comm().rank();
  grid.loadBalance();
  grid.globalRefine(level);
  grid.loadBalance();
  for (int i=0; i<adapt; i++)
  {
    if (myrank==0)
    {
      typedef typename GridType :: template Codim<0> ::
      template Partition<Interior_Partition> :: LeafIterator LeafIterator;

      LeafIterator endit = grid.template leafend<0,Interior_Partition>   ();
      int nr = 0;
      int size = grid.size(0);
      for(LeafIterator it    = grid.template leafbegin<0,Interior_Partition> ();
          it != endit ; ++it,nr++ )
      {
        grid.mark( 1, *it );
        if (nr>size*0.2) break;
      }
    }
    grid.adapt();
    grid.postAdapt();
    grid.loadBalance();
  }
}

template< class GridType >
void checkALUSerial ( GridType &grid, int maxLevel = 1 )
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

  // check also non-conform grids
  makeNonConfGrid(grid,0,1);
  gridcheck(grid);

  // check the method geometryInFather()
  checkGeometryInFather(grid);

  // check the intersection iterator and the geometries it returns
  checkIntersectionIterator(grid);
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

    if( argc < 2 )
    {
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
