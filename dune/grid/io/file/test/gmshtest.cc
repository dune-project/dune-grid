// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"

#include <dune/common/mpihelper.hh>

#include <dune/grid/sgrid.hh>
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#endif
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif

#include <dune/grid/io/file/gmshreader.hh>

using namespace Dune;

template <class GridType>
void testReadingGrid(const std::string& filename) {

  // Read the grid
  std::auto_ptr<GridType> grid(GmshReader<GridType>::read(filename));

}

int main(int argc, char** argv) try {
  // initialize MPI if neccessary
  Dune::MPIHelper::instance(argc, argv);

  const std::string path("../../../../../doc/grids/gmsh/");
  std::string curved2d( path );
  curved2d += "curved2d.msh";
  std::string pyramid( path );
  pyramid += "pyramid.msh";

  // Test whether unstructured grids can be read and written
#if HAVE_UG
  std::cout << "reading UGGrid<2>" << std::endl;
  testReadingGrid<UGGrid<2> >( curved2d );

  std::cout << "reading UGGrid<3>" << std::endl;
  testReadingGrid<UGGrid<3> >( pyramid );
#endif

#if HAVE_ALBERTA
#if ALBERTA_DIM==2
  std::cout << "reading AlbertaGrid<2>" << std::endl;
  testReadingGrid<AlbertaGrid<2> >( curved2d );
#endif

#if ALBERTA_DIM==3
  std::cout << "reading AlbertaGrid<3>" << std::endl;
  testReadingGrid<AlbertaGrid<3> >( pyramid );
#endif
#endif

#if HAVE_ALUGRID
  std::cout << "reading ALUSimplexGrid<2,2>" << std::endl;
  testReadingGrid<ALUSimplexGrid<2,2> >( curved2d );

  std::cout << "reading ALUSimplexGrid<3,3>" << std::endl;
  testReadingGrid<ALUSimplexGrid<3,3> >( pyramid );
#endif

  return 0;

}
catch (Dune::Exception &e) {
  std::cerr << e << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
