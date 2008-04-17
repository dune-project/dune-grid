// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/amirameshwriter.hh>

using namespace Dune;

int main() try {

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
    vertexdata2d[sgrid2d.leafIndexSet().index(*vIt2d)] = vIt2d->geometry()[0].two_norm();

  array<unsigned int, 2> n2;
  n2[0] = n[0]+1;
  n2[1] = n[1]+1;
  // write data
  AmiraMeshWriter<SGrid<2,2>,SGrid<2,2>::Traits::LeafIndexSet>::writeUniformData(sgrid2d, n2, vertexdata2d, "sgrid2d.am");

  // /////////////////////////////////////
  //   Test writing of 3d uniform grid
  // /////////////////////////////////////
  SGrid<3,3> sgrid3d(n,h);

  // create data buffer
  std::vector<double> vertexdata3d(sgrid3d.size(3));

  SGrid<3,3>::Codim<3>::LeafIterator vIt    = sgrid3d.leafbegin<3>();
  SGrid<3,3>::Codim<3>::LeafIterator vEndIt = sgrid3d.leafend<3>();

  for (; vIt!=vEndIt; ++vIt)
    vertexdata3d[sgrid3d.leafIndexSet().index(*vIt)] = vIt->geometry()[0].two_norm();

  array<unsigned int, 3> n3;
  n3[0] = n[0]+1;
  n3[1] = n[1]+1;
  n3[2] = n[2]+1;

  // write data
  AmiraMeshWriter<SGrid<3,3>,SGrid<3,3>::Traits::LeafIndexSet>::writeUniformData(sgrid3d, n3, vertexdata3d, "sgrid3d.am");

  return 0;

}
catch (Dune::Exception &e) {
  std::cerr << e << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
