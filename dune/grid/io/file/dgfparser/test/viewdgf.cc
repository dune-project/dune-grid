// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <string>

#include <dune/grid/io/visual/grapegriddisplay.hh>

template< class Grid >
void display ( const Grid &grid )
{
  Dune::GrapeGridDisplay< Grid > grape( grid, grid.comm().rank() );
  grape.display();
}

int main ( int argc, char **argv )
try
{
  Dune::MPIHelper::instance( argc, argv );

  if( argc < 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <dgffile>" << std::endl;
    return 1;
  }

  typedef Dune::GridSelector::GridType GridType;
  Dune::GridPtr< GridType > gridptr( argv[ 1 ] );
  display( *gridptr );
  return 0;
}
catch( const Dune::Exception &exception )
{
  std::cerr << exception << std::endl;
  return 1;
}
catch( ... )
{
  std :: cerr << "Generic exception!" << std::endl;
  return 1;
}
