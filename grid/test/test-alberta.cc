// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

/*

   Instantiate Alberta-Grid and feed it to the generic gridcheck()

   Note: Albert needs the defines DIM and DIM_OF_WORLD on the
   commandline anyway thus we can use them to select the correct class

 */

#include <iostream>
#include <sstream>

#ifndef GRIDDIM
#define GRIDDIM ALBERTA_DIM
#endif

#include <dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>
#include <dune/grid/albertagrid/albertareader.hh>

#include "gridcheck.cc"
#include "checkgeometryinfather.cc"
#include "checkintersectionit.cc"
#include "checkcommunicate.cc"
#include "checkiterators.cc"
#include "checkadaptation.cc"


template <class GridType >
void markOne ( GridType & grid , int num , int ref )
{
  typedef typename GridType::template Codim<0> :: LeafIterator LeafIterator;

  int count = 0;

  LeafIterator endit = grid.template leafend  <0> ();
  for(LeafIterator it = grid.template leafbegin<0> (); it != endit ; ++it )
  {
    if(num == count) grid.mark( ref, *it );
    count++;
  }

  grid.preAdapt();
  grid.adapt();
  grid.postAdapt();
}


template< class Grid >
void testAlbertaReader ()
{
  std::cout << ">>> Checking AlbertaReader..." << std::endl;

  std::ostringstream filename;
  filename << "grid-" << Grid::dimension << "-" << Grid::dimensionworld << ".amc";

  AlbertaReader< Grid > reader;
  GridFactory< Grid > factory;
  reader.readGrid( filename.str(), factory );

  // create grid and just check the macro grid
  Grid *grid = factory.createGrid();
  gridcheck( *grid );
  GridFactory< Grid >::destroyGrid( grid );
}


int main ()
try {
  const int dim = GRIDDIM;

  typedef Dune::AlbertaGrid< dim > GridType;

  std::cout << "Testing " << GridType::typeName() << "..." << std::endl;

#if DUNE_ALBERTA_VERSION >= 0x200
  testAlbertaReader< GridType >();
#endif

  /* use grid-file appropriate for dimensions */
  std::ostringstream filename;
  filename << "simplex-testgrid-" << GridType::dimension << "-" << GridType::dimensionworld << ".dgf";

  std::cout << std::endl << GridType::typeName() << " with grid file: " << filename.str() << std::endl << std::endl;

  // extra-environment to check destruction
  {
    factorEpsilon = 5e2;

    Dune::GridPtr<GridType> gridPtr(filename.str());
    GridType & grid = *gridPtr;

    // check adaptation interface
    checkAdaptation( grid );

    std::cout << ">>> Checking macro grid..." << std::endl;
    gridcheck(grid); // check macro grid
    checkIterators( grid.leafView() );
    checkIntersectionIterator(grid,true);
    for(int i=0; i<1; i++)
    {
      std::cout << ">>> Refining grid and checking again..." << std::endl;
      grid.globalRefine( 1 );
      gridcheck(grid);
      checkIterators( grid.leafView() );
      checkIntersectionIterator(grid,true);
    }

    // check dgf grid width half refinement
    const int stepsForHalf = DGFGridInfo< GridType >::refineStepsForHalf();
    std::cout << ">>> Refining grid (" << stepsForHalf
              << " times) and checking again..." << std::endl;
    grid.globalRefine( stepsForHalf );
    gridcheck(grid);
    checkIterators( grid.leafView() );
    checkIntersectionIterator(grid,true);

    for(int i=0; i<2; i++)
    {
      std::cout << ">>> Refining one element and checking again..." << std::endl;
      markOne(grid,0,dim);
      gridcheck(grid);
      checkIterators( grid.leafView() );
    }

    checkGeometryInFather(grid);
    checkIntersectionIterator(grid,true);

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
