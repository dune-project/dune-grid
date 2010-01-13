// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <memory>

#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#endif

#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/io/file/gmshreader.hh>

#include <dune/grid/io/visual/grapegriddisplay.hh>

#ifndef ALBERTA_DIM
#define ALBERTA_DIM 2
#endif

#ifndef GRIDDIM
#define GRIDDIM ALBERTA_DIM
#endif

template <class GridType>
void checkGmshReader(const char* filename, const int refinements)
{
  std::auto_ptr< GridType > grid( Dune::GmshReader< GridType >::read( filename ) );
  if( refinements > 0 )
    grid->globalRefine( refinements );

  Dune::GrapeGridDisplay< GridType > grape( *grid );
  grape.display();
}

int main ( int argc, char **argv )
try
{
  if( argc < 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <gmshfile> [refinements]" << std::endl;
    return 1;
  }

  int refinements = 0;
  if( argc >= 3 )
    refinements = atoi( argv[2] );

#if HAVE_ALBERTA
  std::cout << "Checking AlbertaGrid< " << GRIDDIM << " >..." << std::endl;
  //checkGmshReader< Dune::AlbertaGrid< GRIDDIM > > ( argv[1], refinements );
#endif

#if HAVE_ALUGRID
  std::cout << "Checking ALUGrid \n";
  checkGmshReader< Dune::ALUSimplexGrid< GRIDDIM, GRIDDIM > > ( argv[1], refinements );
#endif

#if HAVE_UG && (GRIDDIM == 3)
  std::cout << "Checking UG \n";
  checkGmshReader< Dune::UGGrid< 3 > > ( argv[1], refinements );
#endif

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << "Error: " << e << std::endl;
  return 1;
}
catch( ... )
{
  std::cerr << "Error: Generic exception." << std::endl;
  return 2;
}
