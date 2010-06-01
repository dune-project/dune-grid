// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <string>

#include <dune/grid/io/visual/grapegriddisplay.hh>

using namespace Dune;

int main ( int argc, char **argv )
try
{
  const Dune :: MPIHelper &mpi = MPIHelper :: instance( argc, argv );

  if( argc < 2 )
  {
    std :: cerr << "Usage: " << argv[ 0 ] << " <dgffile>" << std :: endl;
    return 1;
  }

  GridPtr< GridType > gridptr( argv[ 1 ] );
  GrapeGridDisplay< GridType > grape( *gridptr, mpi.rank() );
  grape.display();
  return 0;
}
catch( const Dune :: Exception &exception )
{
  std :: cerr << exception << std::endl;
  return 1;
}
catch( ... )
{
  std :: cerr << "Generic exception!" << std::endl;
  return 1;
}
