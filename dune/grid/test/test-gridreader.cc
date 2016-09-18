#include <config.h>

#include <iostream>
#include <sstream>

#if !HAVE_ALBERTA && !HAVE_UG
error "Alberta or UGGrid must be available for this test"
#endif

#ifndef GRIDDIM
#define GRIDDIM ALBERTA_DIM
#endif

#ifndef WORLDDIM
#define WORLDDIM ALBERTA_DIM
#endif

#include <dune/common/parallel/mpihelper.hh>

#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/albertareader.hh>
#endif

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/io/file/amirameshreader.hh>
#include <dune/grid/io/file/starcdreader.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/dgfparser.hh>

template <typename Grid>
void testReader()
{
  const int dim = Grid::dimension;
  const int dow = Grid::dimensionworld;

#if HAVE_ALBERTA
  {
  std::cout << ">>> Checking AlbertaReader..." << std::endl;
  typedef Dune::AlbertaReader< Grid > Reader;

  std::ostringstream filename;
  filename << DUNE_GRID_EXAMPLE_GRIDS_PATH << "amc/grid-" << dim << "-" << dow << ".amc";
  auto grid = Reader::read(filename.str());
  }
#endif


  {
  std::cout << ">>> Checking DgfReader..." << std::endl;
  typedef Dune::DgfReader< Grid > Reader;

  std::ostringstream filename;
  filename << DUNE_GRID_EXAMPLE_GRIDS_PATH << "dgf/simplex-testgrid-3d.dgf";
  auto grid = Reader::read(filename.str());
  }


  {
  std::cout << ">>> Checking GmshReader..." << std::endl;
  typedef Dune::GmshReader< Grid > Reader;

  std::ostringstream filename;
  filename << DUNE_GRID_EXAMPLE_GRIDS_PATH << "gmsh/simplex-testgrid-3d.msh";
  auto grid = Reader::read(filename.str());
  }


  {
  std::cout << ">>> Checking StarCDReader..." << std::endl;
  typedef Dune::StarCDReader< Grid > Reader;

  std::ostringstream filename;
  filename << DUNE_GRID_EXAMPLE_GRIDS_PATH << "starcd/tets"; // NOTE: give filename without extension
  auto grid = Reader::read(filename.str());
  }


#if HAVE_AMIRAMESH
  {
  std::cout << ">>> Checking AmiraMeshReader..." << std::endl;
  typedef Dune::AmiraMeshReader< Grid > Reader;

  std::ostringstream filename;
  filename << DUNE_GRID_EXAMPLE_GRIDS_PATH << "amiramesh/simplex-testgrid-3d.am";
  auto grid = Reader::read(filename.str());
  }
#endif
}


template <typename Grid>
void testReaderFactory()
{
  const int dim = Grid::dimension;
  const int dow = Grid::dimensionworld;

#if HAVE_ALBERTA
  {
  std::cout << ">>> Checking AlbertaReader..." << std::endl;
  Dune::GridFactory<Grid> factory;
  typedef Dune::AlbertaReader< Grid > Reader;

  std::ostringstream filename;
  filename << DUNE_GRID_EXAMPLE_GRIDS_PATH << "amc/grid-" << dim << "-" << dow << ".amc";
  Reader::read(factory, filename.str());
  auto grid = std::unique_ptr<Grid>{ factory.createGrid() };
  }
#endif


  {
  std::cout << ">>> Checking DgfReader..." << std::endl;
  Dune::GridFactory<Grid> factory;
  typedef Dune::DgfReader< Grid > Reader;

  std::ostringstream filename;
  filename << DUNE_GRID_EXAMPLE_GRIDS_PATH << "dgf/simplex-testgrid-3d.dgf";
//   Reader::read(factory, filename.str()); // NOTE: not supported
//   auto grid = std::unique_ptr<Grid>{ factory.createGrid() };
  }


  {
  std::cout << ">>> Checking GmshReader..." << std::endl;
  Dune::GridFactory<Grid> factory;
  typedef Dune::GmshReader< Grid > Reader;

  std::ostringstream filename;
  filename << DUNE_GRID_EXAMPLE_GRIDS_PATH << "gmsh/simplex-testgrid-3d.msh";
  Reader::read(factory, filename.str());
  auto grid = std::unique_ptr<Grid>{ factory.createGrid() };
  }


  {
  std::cout << ">>> Checking StarCDReader..." << std::endl;
  Dune::GridFactory<Grid> factory;
  typedef Dune::StarCDReader< Grid > Reader;

  std::ostringstream filename;
  filename << DUNE_GRID_EXAMPLE_GRIDS_PATH << "starcd/tets"; // NOTE: give filename without extension
  Reader::read(factory, filename.str());
  auto grid = std::unique_ptr<Grid>{ factory.createGrid() };
  }


#if HAVE_AMIRAMESH
  {
  std::cout << ">>> Checking AmiraMeshReader..." << std::endl;
  Dune::GridFactory<Grid> factory;
  typedef Dune::AmiraMeshReader< Grid > Reader;

  std::ostringstream filename;
  filename << DUNE_GRID_EXAMPLE_GRIDS_PATH << "amiramesh/simplex-testgrid-3d.am";
  Reader::read(factory, filename.str());
  auto grid = std::unique_ptr<Grid>{ factory.createGrid() };
  }
#endif

}


int main ( int argc, char **argv ) try
{
  Dune::MPIHelper::instance( argc, argv );

  const int dim = GRIDDIM;
  const int dow = WORLDDIM;

  static_assert( dim == dow && dow == 3, "dim=dow=3 is assumed for this test!" );

#if HAVE_ALBERTA
  std::cout << "AlbertaGrid:\n----------------------------------" << std::endl;
  testReader< Dune::AlbertaGrid<dim,dow> >();
  testReaderFactory< Dune::AlbertaGrid<dim,dow> >();
#endif

#if HAVE_UG
  std::cout << "UGGrid:\n----------------------------------" << std::endl;
  testReader< Dune::UGGrid<dim> >();
  testReaderFactory< Dune::UGGrid<dim> >();
#endif

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
catch( ... )
{
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
