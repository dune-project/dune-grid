// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// AlbertaReader needs ALBERTA 2.0 or newer
// (otherwise just pass this check)
#if HAVE_ALBERTA

#include <dune/grid/albertagrid/albertareader.hh>

#include "gridcheck.cc"

template< class Grid >
void checkAlbertaReader ()
{
  std::cout << ">>> Checking AlbertaReader..." << std::endl;

  std::ostringstream filename;
  filename << SRCDIR << "grid-" << Grid::dimension << "-" << Grid::dimensionworld << ".amc";

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
            << "because ALBERTA 2.0 or newer is not available." << std::endl;
}

#endif // #if HAVE_ALBERTA
