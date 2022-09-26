// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <sstream>

#ifndef GRIDDIM
#define GRIDDIM ALBERTA_DIM
#endif

#include <dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>

#include <doc/grids/gridfactory/testgrids.hh>

#include "../common/boundaryprojection.hh"

#include "gridcheck.hh"
#include "checkgeometryinfather.hh"
#include "checkgeometry.hh"
#include "checkintersectionit.hh"
#include "checkcommunicate.hh"
#include "checkiterators.hh"
#include "checktwists.hh"
#include "check-albertareader.hh"
#include "checkadaptation.hh"
#include "checkpartition.hh"
#include "checkgridfactory.hh"

#include <doc/grids/gridfactory/testgrids.hh>


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

template< class Grid, int dim >
void addToGridFactory ( Dune::GridFactory< Grid > &factory, Dune::Dim< dim > );

template< class Grid >
void addToGridFactory ( Dune::GridFactory< Grid > &factory, Dune::Dim< 1 > )
{
  Dune::TestGrids::unitLine.addToGridFactory( factory,
      [] ( const typename Dune::TestGrid< 1 >::Vertex &x ) {
        typename Dune::TestGrid< 1 >::Vertex y = x;
        y *= 2.0; y -= typename Dune::TestGrid< 1 >::Vertex( 1.0 );
        return y;
        }
      );
}

template< class Grid >
void addToGridFactory ( Dune::GridFactory< Grid > &factory, Dune::Dim< 2 > )
{
  Dune::TestGrids::kuhn2d.addToGridFactory( factory,
      [] ( const typename Dune::TestGrid< 2 >::Vertex &x ) {
        typename Dune::TestGrid< 2 >::Vertex y = x;
        y *= 2.0; y -= typename Dune::TestGrid< 2 >::Vertex( 1.0 );
        return y;
        }
      );
}

template< class Grid >
void addToGridFactory ( Dune::GridFactory< Grid > &factory, Dune::Dim< 3 > )
{
  Dune::TestGrids::kuhn3d.addToGridFactory( factory,
      [] ( const typename Dune::TestGrid< 3 >::Vertex &x ) {
        typename Dune::TestGrid< 3 >::Vertex y = x;
        y *= 2.0; y -= typename Dune::TestGrid< 3 >::Vertex( 1.0 );
        return y;
        }
      );
}

template< class Grid >
void checkProjectedUnitCube ()
{
  typedef Dune::CircleBoundaryProjection< Grid::dimensionworld > Projection;
  std::cout << ">>> Checking projected unit cube..." << std::endl;
  Dune::GridFactory< Grid > gridFactory;
  addToGridFactory( gridFactory, Dune::Dim< Grid::dimensionworld > () );
  gridFactory.insertBoundaryProjection( new Projection );
  gridFactory.markLongestEdge();
  auto grid = gridFactory.createGrid();
  for( int i = 0; i < 2; ++i )
  {
    grid->globalRefine( Grid::dimension );
    gridcheck( *grid );
  }
}


int main ( int argc, char **argv )
try {
  const int dim = GRIDDIM;

  typedef Dune::AlbertaGrid< dim > GridType;

  std::cout << "Testing " << GridType::typeName() << "..." << std::endl;

#if ALBERTA_DIM == 2 && GRIDDIM == 2
  std::cout << "Check GridFactory ..." <<std::endl;
  Dune::checkGridFactory< GridType >( Dune::TestGrids::kuhn2d );
#endif // #if ALBERTA_DIM == 2 && GRIDDIM == 2

#if ALBERTA_DIM == 3 && GRIDDIM == 3
  std::cout << "Check GridFactory ..." <<std::endl;
  Dune::checkGridFactory< GridType >( Dune::TestGrids::kuhn3d );
#endif // #if ALBERTA_DIM == 3 && GRIDDIM == 3

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

#if ALBERTA_DIM == 3 && GRIDDIM == 3
    // FS#1234: The recursive bisection exhausts the execution stack using the other grid.
    filename = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH "dgf/grid3Y.dgf");
#endif
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

    checkPartitionType( grid.leafGridView() );

    checkIterators( grid.leafGridView() );
    checkIntersectionIterator(grid,true);
    checkTwists( grid.leafGridView(), NoMapTwist() );
    for(int i=0; i<1; i++)
    {
      std::cout << ">>> Refining grid and checking again..." << std::endl;
      grid.globalRefine( 1 );
      gridcheck(grid);
      checkIterators( grid.leafGridView() );
      checkIntersectionIterator(grid,true);
      checkTwists( grid.leafGridView(), NoMapTwist() );
    }

    // check dgf grid width half refinement
    const int stepsForHalf = Dune::DGFGridInfo< GridType >::refineStepsForHalf();
    std::cout << ">>> Refining grid (" << stepsForHalf
              << " times) and checking again..." << std::endl;
    grid.globalRefine( stepsForHalf );
    gridcheck(grid);
    checkIterators( grid.leafGridView() );
    checkIntersectionIterator(grid,true);
    checkTwists( grid.leafGridView(), NoMapTwist() );

    for(int i=0; i<2; i++)
    {
      std::cout << ">>> Refining one element and checking again..." << std::endl;
      markOne(grid,0,dim);
      gridcheck(grid);
      checkIterators( grid.leafGridView() );
    }

    checkGeometryInFather(grid);
    checkIntersectionIterator(grid,true);
    checkTwists( grid.leafGridView(), NoMapTwist() );

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
