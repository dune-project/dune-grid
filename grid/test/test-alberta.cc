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

#include "gridcheck.cc"
#include "checkgeometryinfather.cc"
#include "checkintersectionit.cc"
#include "checkcommunicate.cc"
#include "checkiterators.cc"
#include "checktwists.cc"
#include "check-albertareader.cc"
#include "checkadaptation.cc"


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


int main ( int argc, char **argv )
try {
  const int dim = GRIDDIM;

  typedef Dune::AlbertaGrid< dim > GridType;

  std::cout << "Testing " << GridType::typeName() << "..." << std::endl;

  checkAlbertaReader< GridType >();

  std::string filename;
  if( argc <= 1 )
  {
    /* use grid-file appropriate for dimensions */
    std::ostringstream sfilename;
    sfilename << "simplex-testgrid-" << GridType::dimension << "-" << GridType::dimensionworld << ".dgf";
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

    checkIterators( grid.leafView() );
    checkIntersectionIterator(grid,true);
    checkTwists( grid.leafView(), NoMapTwist() );
    for(int i=0; i<1; i++)
    {
      std::cout << ">>> Refining grid and checking again..." << std::endl;
      grid.globalRefine( 1 );
      gridcheck(grid);
      checkIterators( grid.leafView() );
      checkIntersectionIterator(grid,true);
      checkTwists( grid.leafView(), NoMapTwist() );
    }

    // check dgf grid width half refinement
    const int stepsForHalf = DGFGridInfo< GridType >::refineStepsForHalf();
    std::cout << ">>> Refining grid (" << stepsForHalf
              << " times) and checking again..." << std::endl;
    grid.globalRefine( stepsForHalf );
    gridcheck(grid);
    checkIterators( grid.leafView() );
    checkIntersectionIterator(grid,true);
    checkTwists( grid.leafView(), NoMapTwist() );

    for(int i=0; i<2; i++)
    {
      std::cout << ">>> Refining one element and checking again..." << std::endl;
      markOne(grid,0,dim);
      gridcheck(grid);
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
