// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"
#define DISABLE_DEPRECATED_METHOD_CHECK 1

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
#include <dune/grid/onedgrid.hh>

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

#include <dune/grid/test/gridcheck.cc>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/gmshreader.hh>

using namespace Dune;

#if HAVE_ALBERTA
template< int dim, int dimworld >
struct EnableLevelIntersectionIteratorCheck< Dune::AlbertaGrid< dim, dimworld > >
{
  static const bool v = false;
};
#endif

template <typename GridType>
void testReadingGrid( const std::string& filename, int refinements )
{
  // Read the grid with insertion of boundary segments (they may be the default ones)
  std::auto_ptr<GridType> grid( GmshReader<GridType>::read( filename, true, true ) );

  // load balancing and refinement
  grid->loadBalance();
  if ( refinements > 0 )
    grid->globalRefine( refinements );

  // Do some tests to make sure the grid has been properly read
  gridcheck(*grid);

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

  const std::string path = std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "gmsh/";
  std::string curved2d( path ); curved2d += "curved2d.msh";
  std::string circ2nd(  path ); circ2nd  += "circle2ndorder.msh";
  std::string unitsquare_quads_2x2(path);  unitsquare_quads_2x2 += "unitsquare_quads_2x2.msh";
  std::string sphere(   path ); sphere   += "sphere.msh";
  std::string pyramid(  path ); pyramid  += "pyramid.msh";
  std::string pyr2nd(   path ); pyr2nd   += "pyramid2ndorder.msh";

  // test reading of unstructured grids
#if HAVE_UG
  std::cout << "reading UGGrid<2>" << std::endl;
  testReadingGrid<UGGrid<2> >( curved2d, refinements );

  std::cout << "reading UGGrid<2> with second order boundary approximation" << std::endl;
  testReadingGrid<UGGrid<2> >( circ2nd, refinements );

  std::cout << "reading UGGrid<2>" << std::endl;
  testReadingGrid<UGGrid<2> >( unitsquare_quads_2x2, refinements );

  std::cout << "reading hybrid UGGrid<2>" << std::endl;
  testReadingGrid<UGGrid<2> >( path + "hybrid-testgrid-2d.msh", refinements );

  std::cout << "reading UGGrid<3>" << std::endl;
  testReadingGrid<UGGrid<3> >( pyramid, refinements );

  std::cout << "reading UGGrid<3> with second order boundary approximation" << std::endl;
  testReadingGrid<UGGrid<3> >( pyr2nd,  refinements );

  std::cout << "reading hybrid UGGrid<3>" << std::endl;
  testReadingGrid<UGGrid<3> >( path + "hybrid-testgrid-3d.msh", refinements );
#endif

#if HAVE_ALBERTA
#if ALBERTA_DIM==2
  std::cout << "reading AlbertaGrid<2>" << std::endl;
  testReadingGrid<AlbertaGrid<2> >( curved2d, refinements );
#endif
#if ALBERTA_DIM==3
  std::cout << "reading AlbertaGrid<2>" << std::endl;
  testReadingGrid<AlbertaGrid<2> >( sphere, refinements );
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

  std::cout << "reading OneDGrid" << std::endl;
  testReadingGrid<OneDGrid>( path + "oned-testgrid.msh", refinements );

  return 0;

}
catch ( Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
catch (std::exception &e) {
  std::cerr << e.what << std::endl;
  return 1;
}
catch ( ... )
{
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
