// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <dune/grid/yaspgrid.hh>
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#endif

#include <dune/grid/io/file/amirameshreader.hh>
#include <dune/grid/io/file/amirameshwriter.hh>

using namespace Dune;

template <class GridType>
void testReadingUnstructuredGrid(const std::string& filename) {

  // Read the grid
  std::unique_ptr<GridType> grid = AmiraMeshReader<GridType>::read(filename);

  // Write the grid into a tmp file
  LeafAmiraMeshWriter<GridType> amiramesh(*grid);
  amiramesh.write("tmp.grid", true);

  /** \todo Ideally we should check here whether the two files differ.
      However this is not easy as DUNE may permute the vertices of each element */

  // Remove the following tests once that AmiraMeshReader::read returns a std::unique_ptr
  GridType* gridPtr = AmiraMeshReader<GridType>::read(fileName);
  delete gridPtr;
  std::shared_ptr<GridType> gridShared = AmiraMeshReader<GridType>::read(fileName);
}

void testWritingUniformData() {

  std::array<int,2>     n = { {10, 10} };
  FieldVector<double,2> h = {3.0, 2.0};

  // /////////////////////////////////////
  //   Test writing of 2d uniform grid
  // /////////////////////////////////////
  YaspGrid<2> grid2d(h,n);

  // create data buffer
  std::vector<FieldVector<double,1> > vertexdata2d(grid2d.size(2));

  YaspGrid<2>::Codim<2>::LeafIterator vIt2d    = grid2d.leafbegin<2>();
  YaspGrid<2>::Codim<2>::LeafIterator vEndIt2d = grid2d.leafend<2>();

  for (; vIt2d!=vEndIt2d; ++vIt2d)
    vertexdata2d[grid2d.leafIndexSet().index(*vIt2d)] = vIt2d->geometry().corner(0).two_norm();

  std::array<unsigned int, 2> n2;
  n2[0] = n[0]+1;
  n2[1] = n[1]+1;

  // write data
  AmiraMeshWriter<YaspGrid<2>::LeafGridView> amiramesh2d;
  amiramesh2d.addUniformData(grid2d.leafGridView(), n2, vertexdata2d);
  amiramesh2d.write("structuredgrid2d.am");

  // /////////////////////////////////////
  //   Test writing of 3d uniform grid
  // /////////////////////////////////////

  YaspGrid<3> grid3d({3.0, 2.0, 1.0}, { {10, 10, 10} });

  // create data buffer
  std::vector<FieldVector<double,1> > vertexdata3d(grid3d.size(3));

  YaspGrid<3>::Codim<3>::LeafIterator vIt    = grid3d.leafbegin<3>();
  YaspGrid<3>::Codim<3>::LeafIterator vEndIt = grid3d.leafend<3>();

  for (; vIt!=vEndIt; ++vIt)
    vertexdata3d[grid3d.leafIndexSet().index(*vIt)] = vIt->geometry().corner(0).two_norm();

  std::array<unsigned int, 3> n3;
  n3[0] = n[0]+1;
  n3[1] = n[1]+1;
  n3[2] = n[2]+1;

  // write data
  AmiraMeshWriter<YaspGrid<3>::LeafGridView> amiramesh3d;
  amiramesh3d.addUniformData(grid3d.leafGridView(), n3, vertexdata3d);
  amiramesh3d.write("yaspgrid3d.am");

}

int main() try {

  // Test whether unstructured grids can be read and written
#if HAVE_UG
  std::cout << "reading UGGrid<2>" << std::endl;
  testReadingUnstructuredGrid<UGGrid<2> >(std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "amiramesh/hybrid-testgrid-2d.am");

  std::cout << "reading UGGrid<3>" << std::endl;
  testReadingUnstructuredGrid<UGGrid<3> >(std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "amiramesh/hybrid-testgrid-3d.am");
#endif

#if HAVE_ALBERTA
  std::cout << "reading AlbertaGrid<2>" << std::endl;
  testReadingUnstructuredGrid<AlbertaGrid<2> >(std::string(DUNE_GRID_EXAMPLE_GRIDS_PATH) + "amiramesh/simplex-testgrid-2d.am");
#endif

  // Test whether writing uniform data works
  testWritingUniformData();

  return 0;

}
catch (Dune::Exception &e) {
  std::cerr << e << std::endl;
  return 1;
} catch (std::exception &e) {
  std::cerr << e.what() << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
