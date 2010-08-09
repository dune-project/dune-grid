// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#define NEW_SUBENTITY_NUMBERING 1

#include <iostream>
#include <sstream>

#ifndef GRIDDIM
#define GRIDDIM ALBERTA_DIM
#endif

#include <dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>

#include "basicunitcube.hh"

#include "../common/boundaryprojection.hh"

#include "gridcheck.cc"
#include "checkgeometryinfather.cc"
#include "checkgeometry.cc"
#include "checkintersectionit.cc"
#include "checkcommunicate.cc"
#include "checkiterators.cc"
#include "checktwists.cc"
#include "check-albertareader.cc"
#include "checkadaptation.cc"
#include "checkpartition.cc"


template< int dim, int dimworld >
struct EnableLevelIntersectionIteratorCheck< Dune::AlbertaGrid< dim, dimworld > >
{
  static const bool v = false;
};



template <class GridType >
void markOne ( GridType & grid , int num , int ref )
{
  typedef typename GridType::template Codim<0> :: LeafIterator LeafIterator;

  int count = 0;

  const LeafIterator end = grid.template leafend< 0 >();
  for( LeafIterator it = grid.template leafbegin< 0 >(); it != end ; ++it )
  {
    if( num == count )
      grid.mark( ref, *it );
    ++count;
  }

  grid.preAdapt();
  grid.adapt();
  grid.postAdapt();
}


template< class Grid >
void checkProjectedUnitCube ()
{
  typedef Dune::CircleBoundaryProjection< Grid::dimensionworld > Projection;
  std::cout << ">>> Checking projected unit cube..." << std::endl;
  Dune::GridFactory< Grid > gridFactory;
  BasicUnitCube< Grid::dimension >::insertVertices( gridFactory, -1.0, 1.0 );
  BasicUnitCube< Grid::dimension >::insertSimplices( gridFactory );
  gridFactory.insertBoundaryProjection( new Projection );
  gridFactory.markLongestEdge();
  Grid *grid = gridFactory.createGrid();
  for( int i = 0; i < 2; ++i )
  {
    grid->globalRefine( Grid::dimension );
    gridcheck( *grid );
    checkGeometry( grid->leafView() );
  }
  delete grid;
}


int main ( int argc, char **argv )
try {
  const int dim = GRIDDIM;

  typedef Dune::AlbertaGrid< dim > GridType;

  std::cout << "Testing " << GridType::typeName() << "..." << std::endl;

  std::string filename;
  if( argc <= 1 )
  {
    checkAlbertaReader< GridType >();

#if ALBERTA_DIM == GRIDDIM
    if( dim < 3 )
      checkProjectedUnitCube< GridType >();
#endif

    /* use grid-file appropriate for dimensions */
    std::ostringstream sfilename;
    const int dimWorld = GridType::dimensionworld;
    sfilename << DUNE_GRID_EXAMPLE_GRIDS_PATH "dgf/simplex-testgrid-" << GridType::dimension << "-"
              << std::min( dimWorld, 3) << ".dgf";
    filename = sfilename.str();
  }
  else
    filename = argv[ 1 ];

  std::cout << std::endl << GridType::typeName() << " with grid file: " << filename << std::endl << std::endl;
  {
    factorEpsilon = 5e2;

    Dune::GridPtr<GridType> gridPtr( filename );
    GridType & grid = *gridPtr;

    // extra-environment to check destruction

    std::cout << ">>> Checking macro grid..." << std::endl;

    gridcheck(grid); // check macro grid

    // check grid adaptation interface
    checkAdaptation( grid );

    checkPartitionType( grid.leafView() );

    checkGeometry( grid.leafView() );
    checkIterators( grid.leafView() );
    checkIntersectionIterator(grid,true);
    checkTwists( grid.leafView(), NoMapTwist() );
    for(int i=0; i<1; i++)
    {
      std::cout << ">>> Refining grid and checking again..." << std::endl;
      grid.globalRefine( 1 );
      gridcheck(grid);
      checkGeometry( grid.leafView() );
      checkIterators( grid.leafView() );
      checkIntersectionIterator(grid,true);
      checkTwists( grid.leafView(), NoMapTwist() );
    }

    // check dgf grid width half refinement
    const int stepsForHalf = Dune::DGFGridInfo< GridType >::refineStepsForHalf();
    std::cout << ">>> Refining grid (" << stepsForHalf
              << " times) and checking again..." << std::endl;
    grid.globalRefine( stepsForHalf );
    gridcheck(grid);
    checkGeometry( grid.leafView() );
    checkIterators( grid.leafView() );
    checkIntersectionIterator(grid,true);
    checkTwists( grid.leafView(), NoMapTwist() );

    for(int i=0; i<2; i++)
    {
      std::cout << ">>> Refining one element and checking again..." << std::endl;
      markOne(grid,0,dim);
      gridcheck(grid);
      checkGeometry( grid.leafView() );
      checkIterators( grid.leafView() );
    }

    checkGeometryInFather(grid);
    checkIntersectionIterator(grid,true);
    checkTwists( grid.leafView(), NoMapTwist() );

    checkCommunication(grid, -1, Dune::dvverb);
  };

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
catch( ... )
{
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
