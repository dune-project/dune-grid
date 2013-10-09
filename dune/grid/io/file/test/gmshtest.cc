// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"
#define DISABLE_DEPRECATED_METHOD_CHECK 1

#include <dune/common/parallel/mpihelper.hh>

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
#include <dune/grid/io/file/gmshwriter.hh>

using namespace Dune;

#if HAVE_ALBERTA
template< int dim, int dimworld >
struct EnableLevelIntersectionIteratorCheck< Dune::AlbertaGrid< dim, dimworld > >
{
  static const bool v = false;
};
#endif

template <typename GridType>
void testReadingAndWritingGrid( const std::string& filename, const std::string& outFilename, int refinements )
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

  // Test writing
  Dune::GmshWriter<typename GridType::LeafGridView> writer( grid->leafView() );
  writer.write( outFilename );

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
  std::string sphere(   path ); sphere    += "sphere.msh";
  std::string pyramid(  path ); pyramid   += "pyramid.msh";
  std::string pyr2nd(   path ); pyr2nd    += "pyramid2ndorder.msh";
  std::string hybrid_2d( path); hybrid_2d += "hybrid-testgrid-2d.msh";
  std::string hybrid_3d( path); hybrid_3d += "hybrid-testgrid-3d.msh";
  std::string oned(      path); oned += "oned-testgrid.msh";

  // test reading and writing of unstructured grids
#if HAVE_UG
  std::cout << "reading and writing UGGrid<2>" << std::endl;
  testReadingAndWritingGrid<UGGrid<2> >( curved2d, curved2d+".UGGrid_2_-gmshtest-write.msh", refinements );

  std::cout << "reading and writing UGGrid<2> with second order boundary approximation" << std::endl;
  testReadingAndWritingGrid<UGGrid<2> >( circ2nd, circ2nd+".UGGrid_2_-gmshtest-write.msh", refinements );

  std::cout << "reading and writing UGGrid<2>" << std::endl;
  testReadingAndWritingGrid<UGGrid<2> >( unitsquare_quads_2x2, unitsquare_quads_2x2+".UGGrid_2_-gmshtest-write.msh", refinements );

  std::cout << "reading and writing hybrid UGGrid<2>" << std::endl;
  testReadingAndWritingGrid<UGGrid<2> >( hybrid_2d, hybrid_2d+".UGGrid_2_-gmshtest-write.msh", refinements );

  std::cout << "reading and writing UGGrid<3>" << std::endl;
  testReadingAndWritingGrid<UGGrid<3> >( pyramid, pyramid+".UGGrid_3_-gmshtest-write.msh", refinements );

  std::cout << "reading and writing UGGrid<3> with second order boundary approximation" << std::endl;
  testReadingAndWritingGrid<UGGrid<3> >( pyr2nd,  pyr2nd+".UGGrid_3_-gmshtest-write.msh", refinements );

  std::cout << "reading and writing hybrid UGGrid<3>" << std::endl;
  testReadingAndWritingGrid<UGGrid<3> >( hybrid_3d, hybrid_3d+".UGGrid_3_-gmshtest-write.msh", refinements );
#endif

#if HAVE_ALBERTA
#if ALBERTA_DIM==2
  std::cout << "reading and writing AlbertaGrid<2>" << std::endl;
  testReadingAndWritingGrid<AlbertaGrid<2> >( curved2d, curved2d+".AlbertaGrid_2_-gmshtest-write.msh", refinements );
#endif
#if ALBERTA_DIM==3
  std::cout << "reading and writing AlbertaGrid<2>" << std::endl;
  testReadingAndWritingGrid<AlbertaGrid<2> >( sphere, sphere+".AlbertaGrid_2_-gmshtest-write.msh", refinements );
  std::cout << "reading and writing AlbertaGrid<3>" << std::endl;
  testReadingAndWritingGrid<AlbertaGrid<3> >( pyramid, pyramid+".AlbertaGrid_3_-gmshtest-write.msh", refinements );
#endif
#endif

#if HAVE_ALUGRID
  std::cout << "reading and writing ALUSimplexGrid<2,2>" << std::endl;
  testReadingAndWritingGrid<ALUSimplexGrid<2,2> >( curved2d, curved2d+".ALUSimplexGrid_2_2_-gmshtest-write.msh", refinements );

  std::cout << "reading and writing ALUSimplexGrid<3,3>" << std::endl;
  testReadingAndWritingGrid<ALUSimplexGrid<3,3> >( pyramid, pyramid+".ALUSimplexGrid_3_3_-gmshtest-write.msh", refinements );
#endif

  std::cout << "reading and writing OneDGrid" << std::endl;
  testReadingAndWritingGrid<OneDGrid>( oned, oned+".OneDGrid-gmshtest-write.msh", refinements );


  return 0;

}
catch ( Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
catch (std::exception &e) {
  std::cerr << e.what() << std::endl;
  return 1;
}
catch ( ... )
{
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
