// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_CHECK_ALBERTAREADER_HH
#define DUNE_GRID_TEST_CHECK_ALBERTAREADER_HH

#include <iostream>

// AlbertaReader needs ALBERTA 3.0 or newer
// (otherwise just pass this check)
#if HAVE_ALBERTA

#include <dune/grid/albertagrid/albertareader.hh>

#include "gridcheck.hh"

template< class Grid >
void checkAlbertaReader ()
{
  std::cout << ">>> Checking AlbertaReader..." << std::endl;

  std::ostringstream filename;
  filename << DUNE_GRID_EXAMPLE_GRIDS_PATH << "amc/grid-" << Grid::dimension << "-" << Grid::dimensionworld << ".amc";

  Dune::AlbertaReader< Grid > reader;
  Dune::GridFactory< Grid > factory;
  reader.readGrid( filename.str(), factory );

  // create grid and just check the macro grid
  Grid *grid = factory.createGrid();
  grid->globalRefine( 2 );
  gridcheck( *grid );
  //GridFactory< Grid >::destroyGrid( grid );
  delete grid;
}

#else

template< class Grid >
void checkAlbertaReader ()
{
  std::cerr << "Warning: Skipping AlbertaReader check, "
            << "because ALBERTA is not available." << std::endl;
}

#endif // #if HAVE_ALBERTA

#endif // DUNE_GRID_TEST_CHECK_ALBERTAREADER_HH
