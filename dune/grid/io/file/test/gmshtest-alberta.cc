// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>

#include <memory>

#include <dune/common/mpihelper.hh>

#undef HAVE_ALBERTA

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

#if HAVE_GRAPE
#include <dune/grid/io/visual/grapegriddisplay.hh>
#endif // #if HAVE_GRAPE

#ifndef ALBERTA_DIM
#define ALBERTA_DIM 2
#endif
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#ifndef GRIDDIM
#define GRIDDIM ALBERTA_DIM
#endif

template <class GridType>
void checkGmshReader(const char* filename, const int refinements)
{
  std::auto_ptr< GridType > grid( Dune::GmshReader< GridType >::read( filename ) );
  grid->loadBalance();
  if( refinements > 0 )
    grid->globalRefine( refinements );

#if HAVE_GRAPE
  Dune::GrapeGridDisplay< GridType > grape( *grid );
  grape.display();
#endif // #if HAVE_GRAPE

  std::ostringstream vtkName;
  vtkName << filename << "-" << refinements;
  Dune::VTKWriter< typename GridType::LeafGridView > vtkWriter( grid->leafView() );
  vtkWriter.write( vtkName.str() );
}

void runChecks(const std::string& filename, unsigned refinements) {
  std::cout << "running gmshtest with gmsh files \"" << filename << "\" and "
            << refinements << " refinements." << std::endl;

#if HAVE_ALBERTA
  std::cout << "Checking AlbertaGrid< " << GRIDDIM << " >..." << std::endl;
  checkGmshReader<Dune::AlbertaGrid<GRIDDIM> >(filename, refinements);
#endif

#if HAVE_ALUGRID
  std::cout << "Checking ALUGrid \n";
  checkGmshReader<Dune::ALUSimplexGrid<GRIDDIM, GRIDDIM> >(filename,
                                                           refinements);
#endif

#if HAVE_UG && (GRIDDIM == 3)
  std::cout << "Checking UG \n";
  checkGmshReader<Dune::UGGrid<3> >(filename, refinements);
#endif
}

int main ( int argc, char **argv )
try
{
  Dune::MPIHelper::instance( argc, argv );
  int refinements = 0;

  if( argc < 2 )
  {
    const std::string path("../../../../../doc/grids/gmsh/");

    std::string curved2d( path );
    curved2d += "curved2d.msh";
    runChecks(curved2d, refinements);

    std::string pyramid( path );
    pyramid += "pyramid.msh";
    runChecks(pyramid, refinements);

    return 0;
  }

  if( argc >= 3 )
    refinements = atoi( argv[2] );

  runChecks(argv[1], refinements);

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
