// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#define ENABLE_NORMAL_CHECK 1

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
#include <dune/grid/io/file/dgfparser/dgfalberta.hh>

#include "gridcheck.cc"
#include "checkgeometryinfather.cc"
#include "checkintersectionit.cc"
#include "checkcommunicate.cc"



template <class GridType >
void markOne ( GridType & grid , int num , int ref )
{
  typedef typename GridType::template Codim<0> :: LeafIterator LeafIterator;

  int count = 0;

  LeafIterator endit = grid.template leafend  <0> ();
  for(LeafIterator it = grid.template leafbegin<0> (); it != endit ; ++it )
  {
    if(num == count) grid.mark( ref, it );
    count++;
  }

  grid.preAdapt();
  grid.adapt();
  grid.postAdapt();
}

int main () {
  try {
    const int dim      = GRIDDIM;
    const int dimworld = GRIDDIM;

    /* use grid-file appropriate for dimensions */
    std::ostringstream filename;

    filename << "simplex-testgrid-" << dim
             << "-" << dimworld << ".dgf";

    std::cout << std::endl << "AlbertaGrid<" << dim
              << "," << dimworld
              << "> with grid file: " << filename.str()
              << std::endl << std::endl;
    {
      factorEpsilon = 5e2;

      typedef Dune::AlbertaGrid<dim,dimworld> GridType;

      Dune::GridPtr<GridType> gridPtr(filename.str());
      GridType & grid = *gridPtr;

      // extra-environment to check destruction

      gridcheck(grid); // check macro grid
      for(int i=0; i<1; i++)
      {
        // check single refinement
        grid.globalRefine( 1 );
        gridcheck(grid);
        checkIntersectionIterator(grid,true);
      }

      // check dgf grid width half refinement
      grid.globalRefine( DGFGridInfo<GridType> :: refineStepsForHalf() );
      gridcheck(grid);
      checkIntersectionIterator(grid,true);

      for(int i=0; i<2; i++)
      {
        markOne(grid,0,dim);
        gridcheck(grid);
      }

      checkGeometryInFather(grid);
      checkIntersectionIterator(grid,true);

      checkCommunication(grid, -1, Dune::dvverb);
    };
  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
