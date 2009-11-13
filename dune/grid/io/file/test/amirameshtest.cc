// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

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

#include <dune/grid/io/file/amirameshreader.hh>
#include <dune/grid/io/file/amirameshwriter.hh>

using namespace Dune;

template <class GridType>
void testReadingUnstructuredGrid(const std::string& filename) {

  // Read the grid
  std::auto_ptr<GridType> grid(AmiraMeshReader<GridType>::read(filename));

  // Write the grid into a tmp file
  LeafAmiraMeshWriter<GridType> amiramesh(*grid);
  amiramesh.write("tmp.grid", true);

  /** \todo Ideally we should check here whether the two files differ.
      However this is not easy as DUNE may permute the vertices of each element */

}

void testWritingUniformData() {

  int n[3] = { 10, 10, 10 };
  double h[3] = {3.0, 2.0, 1.0 };

  // /////////////////////////////////////
  //   Test writing of 2d uniform grid
  // /////////////////////////////////////
  SGrid<2,2> sgrid2d(n,h);

  // create data buffer
  std::vector<double> vertexdata2d(sgrid2d.size(2));

  SGrid<2,2>::Codim<2>::LeafIterator vIt2d    = sgrid2d.leafbegin<2>();
  SGrid<2,2>::Codim<2>::LeafIterator vEndIt2d = sgrid2d.leafend<2>();

  for (; vIt2d!=vEndIt2d; ++vIt2d)
    vertexdata2d[sgrid2d.leafIndexSet().index(*vIt2d)] = vIt2d->geometry().corner(0).two_norm();

  array<unsigned int, 2> n2;
  n2[0] = n[0]+1;
  n2[1] = n[1]+1;

  // write data
  AmiraMeshWriter<SGrid<2,2>::LeafGridView> amiramesh2d;
  amiramesh2d.addUniformData(sgrid2d.leafView(), n2, vertexdata2d);
  amiramesh2d.write("sgrid2d.am");

  // /////////////////////////////////////
  //   Test writing of 3d uniform grid
  // /////////////////////////////////////
  SGrid<3,3> sgrid3d(n,h);

  // create data buffer
  std::vector<double> vertexdata3d(sgrid3d.size(3));

  SGrid<3,3>::Codim<3>::LeafIterator vIt    = sgrid3d.leafbegin<3>();
  SGrid<3,3>::Codim<3>::LeafIterator vEndIt = sgrid3d.leafend<3>();

  for (; vIt!=vEndIt; ++vIt)
    vertexdata3d[sgrid3d.leafIndexSet().index(*vIt)] = vIt->geometry().corner(0).two_norm();

  array<unsigned int, 3> n3;
  n3[0] = n[0]+1;
  n3[1] = n[1]+1;
  n3[2] = n[2]+1;

  // write data
  AmiraMeshWriter<SGrid<3,3>::LeafGridView> amiramesh3d;
  amiramesh3d.addUniformData(sgrid3d.leafView(), n3, vertexdata3d);
  amiramesh3d.write("sgrid3d.am");

}

int main() try {

  // Test whether unstructured grids can be read and written
#if HAVE_UG
  std::cout << "reading UGGrid<2>" << std::endl;
  testReadingUnstructuredGrid<UGGrid<2> >("../../../../doc/grids/amiramesh/hybrid-testgrid-2d.am");

  std::cout << "reading UGGrid<3>" << std::endl;
  testReadingUnstructuredGrid<UGGrid<3> >("../../../../doc/grids/amiramesh/hybrid-testgrid-3d.am");
#endif

#if HAVE_ALBERTA
  std::cout << "reading AlbertaGrid<2>" << std::endl;
  testReadingUnstructuredGrid<AlbertaGrid<2> >("../../../../doc/grids/amiramesh/simplex-testgrid-2d.am");

  std::cout << "reading AlbertaGrid<3>" << std::endl;
  testReadingUnstructuredGrid<AlbertaGrid<3> >("../../../../doc/grids/amiramesh/simplex-testgrid-3d.am");
#endif

#if HAVE_ALUGRID
  std::cout << "reading ALUSimplexGrid<2,2>" << std::endl;
  testReadingUnstructuredGrid<ALUSimplexGrid<2,2> >("../../../../doc/grids/amiramesh/simplex-testgrid-2d.am");

  std::cout << "reading ALUSimplexGrid<3,3>" << std::endl;
  testReadingUnstructuredGrid<ALUSimplexGrid<3,3> >("../../../../doc/grids/amiramesh/simplex-testgrid-3d.am");

  std::cout << "reading ALUCubeGrid<3,3>" << std::endl;
  testReadingUnstructuredGrid<ALUCubeGrid<3,3> >("../../../../doc/grids/amiramesh/cube-testgrid-3d.am");
#endif

  // Test whether writing uniform data works
  testWritingUniformData();

  return 0;

}
catch (Dune::Exception &e) {
  std::cerr << e << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
