// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <dune/grid/sgrid.hh>
#ifdef HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#ifdef HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#endif
#ifdef HAVE_ALUGRID
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

  // Test whether unstructured grids can be read and written
#ifdef HAVE_UG
  std::cout << "reading UGGrid<2>" << std::endl;
  testReadingGrid<UGGrid<2> >("../../../../doc/grids/gmsh/curved2d.msh");

  std::cout << "reading UGGrid<3>" << std::endl;
  testReadingGrid<UGGrid<3> >("../../../../doc/grids/gmsh/pyramid.msh");
#endif

#ifdef HAVE_ALBERTA
  std::cout << "reading AlbertaGrid<2>" << std::endl;
  testReadingGrid<AlbertaGrid<2> >("../../../../doc/grids/gmsh/curved2d.msh");

  std::cout << "reading AlbertaGrid<3>" << std::endl;
  testReadingGrid<AlbertaGrid<3> >("../../../../doc/grids/gmsh/pyramid.msh");
#endif

#ifdef HAVE_ALUGRID
  //     std::cout << "reading ALUSimplexGrid<2,2>" << std::endl;
  //     testReadingGrid<ALUSimplexGrid<2,2> >("../../../../doc/grids/gmsh/curved2d.msh");

  std::cout << "reading ALUSimplexGrid<3,3>" << std::endl;
  testReadingGrid<ALUSimplexGrid<3,3> >("../../../../doc/grids/gmsh/pyramid.msh");
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
