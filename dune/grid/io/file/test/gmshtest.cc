// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <dune/grid/sgrid.hh>
#ifdef ENABLE_UG
#include <dune/grid/uggrid.hh>
#endif
#ifdef ENABLE_ALBERTA
#include <dune/grid/albertagrid.hh>
#endif
#ifdef ENABLE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif

#include <dune/grid/io/file/gmshreader.hh>

using namespace Dune;

template <class GridType>
void testReadingGrid(const std::string& filename) {

  // Read the grid
  std::auto_ptr<GridType> grid(GmshReader<GridType>::read(filename));

}

int main() try {

  const std::string path("../../../../../doc/grids/gmsh/");
  std::string curved2d( path );
  curved2d += "curved2d.msh";
  std::string pyramid( path );
  pyramid += "pyramid.msh";

  // Test whether unstructured grids can be read and written
#ifdef ENABLE_UG
  std::cout << "reading UGGrid<2>" << std::endl;
  testReadingGrid<UGGrid<2> >( curved2d );

  std::cout << "reading UGGrid<3>" << std::endl;
  testReadingGrid<UGGrid<3> >( pyramid );
#endif

#ifdef ENABLE_ALBERTA
  std::cout << "reading AlbertaGrid<2>" << std::endl;
  testReadingGrid<AlbertaGrid<2> >( curved2d );

  std::cout << "reading AlbertaGrid<3>" << std::endl;
  testReadingGrid<AlbertaGrid<3> >( pyramid );
#endif

#ifdef ENABLE_ALUGRID
  //     std::cout << "reading ALUSimplexGrid<2,2>" << std::endl;
  //     testReadingGrid<ALUSimplexGrid<2,2> >( cruved2d );

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
