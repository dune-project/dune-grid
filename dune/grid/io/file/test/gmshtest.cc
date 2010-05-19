// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"

#include <dune/common/mpihelper.hh>

// dune grid includes
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#endif
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif

// alberta related stuff
#ifndef ALBERTA_DIM
#define ALBERTA_DIM 2
#endif

#ifndef GRIDDIM
#define GRIDDIM ALBERTA_DIM
#endif

// grape include
#if HAVE_GRAPE
#include <dune/grid/io/visual/grapegriddisplay.hh>
#endif

// include VTK writer
#include <dune/grid/io/file/vtk/vtkwriter.hh>

// include gmshreader
#include <dune/grid/io/file/gmshreader.hh>

using namespace Dune;

template <typename GridType>
void testReadingGrid( const std::string& filename, int refinements )
{
  // Read the grid with insertion of boundary segments (they may be the default ones)
  std::auto_ptr<GridType> grid( GmshReader<GridType>::read( filename, true, true ) );

  // load balancing and refinement
  grid->loadBalance();
  if ( refinements > 0 )
    grid->globalRefine( refinements );

  // grape output
#if HAVE_GRAPE && USEGRAPE
  Dune::GrapeGridDisplay<GridType> grape( *grid );
  grape.display();
#endif // #if HAVE_GRAPE

  // vtk output
  std::ostringstream vtkName;
  vtkName << filename << "-" << refinements;
  Dune::VTKWriter<typename GridType::LeafGridView> vtkWriter( grid->leafView() );
  vtkWriter.write( vtkName.str() );
}

int main( int argc, char** argv )
try
{
  Dune::MPIHelper::instance( argc, argv );
  int refinements = 0;

  if ( argc > 1 )
    refinements = atoi( argv[1] );

  const std::string path( "../../../../../doc/grids/gmsh/" );
  std::string curved2d( path ); curved2d += "curved2d.msh";
  std::string pyramid(  path ); pyramid  += "pyramid.msh";

  // test reading of unstructured grids
#if HAVE_UG
  std::cout << "reading UGGrid<2>" << std::endl;
  testReadingGrid<UGGrid<2> >( curved2d, refinements );

  std::cout << "reading UGGrid<3>" << std::endl;
  testReadingGrid<UGGrid<3> >( pyramid, refinements );
#endif

#if HAVE_ALBERTA
#if ALBERTA_DIM==2
  std::cout << "reading AlbertaGrid<2>" << std::endl;
  testReadingGrid<AlbertaGrid<2> >( curved2d, refinements );
#endif
#if ALBERTA_DIM==3
  std::cout << "reading AlbertaGrid<3>" << std::endl;
  testReadingGrid<AlbertaGrid<3> >( pyramid, refinements );
#endif
#endif

#if HAVE_ALUGRID
  std::cout << "reading ALUSimplexGrid<2,2>" << std::endl;
  testReadingGrid<ALUSimplexGrid<2,2> >( curved2d, refinements );

  std::cout << "reading ALUSimplexGrid<3,3>" << std::endl;
  testReadingGrid<ALUSimplexGrid<3,3> >( pyramid, refinements );
#endif

#if !defined(HAVE_UG) && !defined(HAVE_ALUGRID) && !(HAVE_ALBERTA)
  // signal 'skipped' to the test system
  return 77;
#endif

  return 0;

}
catch ( Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
catch ( ... )
{
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
