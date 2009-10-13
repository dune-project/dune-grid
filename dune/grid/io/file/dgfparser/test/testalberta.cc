// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#define NEW_SUBENTITY_NUMBERING 1

#include <iostream>
#include <string>

#if HAVE_GRAPE
#include <dune/grid/io/visual/grapegriddisplay.hh>
#endif

#include <dune/grid/albertagrid/dgfparser.hh>
#include <dune/grid/test/gridcheck.cc>

using namespace Dune;

template< int dim, int dimworld >
struct EnableLevelIntersectionIteratorCheck< AlbertaGrid< dim, dimworld > >
{
  static const bool v = false;
};

int main(int argc, char ** argv, char ** envp)
try
{
  MPIHelper::instance( argc, argv );

  typedef AlbertaGrid< ALBERTA_DIM, ALBERTA_DIM > GridType;

  std::string filename;
  if( argc > 1 )
    filename = argv[ 1 ];
  else if( ALBERTA_DIM == 2 )
    filename = SRCDIR "/examplegrid5.dgf";
  else if( ALBERTA_DIM == 3 )
    filename = SRCDIR "/examplegrid6.dgf";

  GridPtr<GridType> gridptr( filename );
  if( argc > 2 )
    gridptr->globalRefine( atoi( argv[ 2 ] ) );

  gridcheck( *gridptr );

#if HAVE_GRAPE
  GrapeGridDisplay< GridType > grape( *gridptr );
  grape.display();
#endif

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << std::endl << "Error: " << e << std::endl;
  return 1;
}
catch (...)
{
  std::cerr << std::endl << "Error: Unknown exception caught" << std::endl;
  return 1;
}
