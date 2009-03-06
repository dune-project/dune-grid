// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// AlbertaReader needs ALBERTA 2.0 or newer
// (otherwise just pass this check)
#if HAVE_ALBERTA
#if DUNE_ALBERTA_VERSION >= 0x200

#include <dune/grid/albertagrid/albertareader.hh>

#include "gridcheck.cc"

template< class Grid >
void checkAlbertaReader ()
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
  //GridFactory< Grid >::destroyGrid( grid );
  delete grid;
}

#else

template< class Grid >
void checkAlbertaReader ()
{
  std::cerr << "Warning: Skipping AlbertaReader check, "
            << "because AlbertaReader is not available for ALBERTA 1.2." << std::endl;
}

#endif // #if DUNE_ALBERTA_VERSION >= 0x200

#else

template< class Grid >
void checkAlbertaReader ()
{
  std::cerr << "Warning: Skipping AlbertaReader check, "
            << "because ALBERTA 2.0 or newer is not available." << std::endl;
}

#endif // #if HAVE_ALBERTA
