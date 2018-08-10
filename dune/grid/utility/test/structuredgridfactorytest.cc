// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \file
    \brief A unit test for the StructuredGridFactory
 */

#include <config.h>

#include <cassert>
#include <iostream>
#include <memory>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/test/gridcheck.hh>

using namespace Dune;


int main (int argc , char **argv)
try {

  // this method calls MPI_Init, if MPI is enabled
  MPIHelper::instance(argc,argv);

  // /////////////////////////////////////////////////////////////////////////////
  //   Test 1d grids
  // /////////////////////////////////////////////////////////////////////////////

  // Test creation of 1d cube grids
  std::array<unsigned int,1> elements1d;
  elements1d.fill(4);

  std::shared_ptr<OneDGrid> onedCubeGrid = StructuredGridFactory<OneDGrid>::createCubeGrid(FieldVector<double,1>(0),
                                                                                      FieldVector<double,1>(1),
                                                                                      elements1d);

  assert(onedCubeGrid->size(1) == static_cast< int >(elements1d[0]+1));
  assert(onedCubeGrid->size(0) == static_cast< int >(elements1d[0]));

  gridcheck(*onedCubeGrid);


  // Test creation of 1d simplex grids
  std::shared_ptr<OneDGrid> onedSimplexGrid = StructuredGridFactory<OneDGrid>::createSimplexGrid(FieldVector<double,1>(0),
                                                                                            FieldVector<double,1>(1),
                                                                                            elements1d);

  assert(onedCubeGrid->size(1) == static_cast< int >(elements1d[0]+1));
  assert(onedCubeGrid->size(0) == static_cast< int >(elements1d[0]));

  gridcheck(*onedSimplexGrid);

  // Test creation of 1d YaspGrid
  std::shared_ptr<YaspGrid<1> > yaspGrid1d =
    StructuredGridFactory<YaspGrid<1> >::createCubeGrid(FieldVector<double,1>(0),
                                                        FieldVector<double,1>(1),
                                                        elements1d);

  assert(yaspGrid1d->size(1) == static_cast< int >(elements1d[0]+1));
  assert(yaspGrid1d->size(0) == static_cast< int >(elements1d[0]));

  gridcheck(*yaspGrid1d);

  {
    std::shared_ptr<YaspGrid<1, EquidistantOffsetCoordinates<double,1> > > yaspGridOff1d
      = StructuredGridFactory<YaspGrid<1, EquidistantOffsetCoordinates<double,1> > >::createCubeGrid(FieldVector<double,1>(-1.),
                                                                FieldVector<double,1>(1),
                                                                elements1d);

    assert(yaspGridOff1d->size(1) == static_cast< int >(elements1d[0]+1));
    assert(yaspGridOff1d->size(0) == static_cast< int >(elements1d[0]));

    gridcheck(*yaspGridOff1d);
  }

  // /////////////////////////////////////////////////////////////////////////////
  //   Test 2d grids
  // /////////////////////////////////////////////////////////////////////////////

  std::array<unsigned int,2> elements2d;
  elements2d.fill(4);
  unsigned int numVertices2d = (elements2d[0]+1) * (elements2d[1]+1);
  unsigned int numCubes2d    = elements2d[0] * elements2d[1];

  // Test creation of 2d YaspGrid
  std::shared_ptr<YaspGrid<2> > yaspGrid2d
    = StructuredGridFactory<YaspGrid<2> >::createCubeGrid(FieldVector<double,2>(0),
                                                          FieldVector<double,2>(1),
                                                          elements2d);

  assert(yaspGrid2d->size(2) == static_cast< int >(numVertices2d));
  assert(yaspGrid2d->size(0) == static_cast< int >(numCubes2d));

  gridcheck(*yaspGrid2d);

  {
    std::shared_ptr<YaspGrid<2, EquidistantOffsetCoordinates<double,2> > > yaspGridOff2d
      = StructuredGridFactory<YaspGrid<2, EquidistantOffsetCoordinates<double,2> > >::createCubeGrid(FieldVector<double,2>(0),
                                                              FieldVector<double,2>(1),
                                                              elements2d);

      assert(yaspGridOff2d->size(2) == static_cast< int >(numVertices2d));
      assert(yaspGridOff2d->size(0) == static_cast< int >(numCubes2d));

      gridcheck(*yaspGridOff2d);
  }

  // Test creation of 2d cube grid using UG
#if HAVE_UG
  typedef UGGrid<2> QuadrilateralGridType;

  std::shared_ptr<QuadrilateralGridType> quadrilateralGrid = StructuredGridFactory<QuadrilateralGridType>::createCubeGrid(FieldVector<double,2>(0),
                                                                                                                     FieldVector<double,2>(1),
                                                                                                                     elements2d);

  assert(quadrilateralGrid->size(2) == numVertices2d);
  assert(quadrilateralGrid->size(0) == numCubes2d);

  gridcheck(*quadrilateralGrid);
#ifdef ModelP  // parallel UGGrid can only have one grid at a time
  quadrilateralGrid.reset();
#endif
#else
  std::cout << "WARNING: 2d cube grids not tested because no suitable grid implementation is available!" << std::endl;
#endif


  // Test creation of 2d triangle grid using UG
#if HAVE_UG
  typedef UGGrid<2> TriangleGridType;

  std::shared_ptr<TriangleGridType> triangleGrid = StructuredGridFactory<TriangleGridType>::createSimplexGrid(FieldVector<double,2>(0),
                                                                                                         FieldVector<double,2>(1),
                                                                                                         elements2d);

  assert(triangleGrid->size(2) == numVertices2d);
  assert(triangleGrid->size(0) == 2*numCubes2d);    // each cube gets split into 2 triangles

  gridcheck(*triangleGrid);
#ifdef ModelP  // parallel UGGrid can only have one grid at a time
  triangleGrid.reset();
#endif
#else
  std::cout << "WARNING: 2d simplicial grids not tested because no suitable grid implementation is available!" << std::endl;
#endif

  // /////////////////////////////////////////////////////////////////////////////
  //   Test 3d grids
  // /////////////////////////////////////////////////////////////////////////////

  // Test creation of 3d cube grids
#if HAVE_UG
  typedef UGGrid<3> HexahedralGridType;

  std::array<unsigned int,3> elements3d;
  elements3d.fill(4);
  unsigned int numVertices3d = (elements3d[0]+1) * (elements3d[1]+1) * (elements3d[2]+1);
  unsigned int numCubes3d    = elements3d[0] * elements3d[1] * elements3d[2];

  std::shared_ptr<HexahedralGridType> hexahedralGrid = StructuredGridFactory<HexahedralGridType>::createCubeGrid(FieldVector<double,3>(0),
                                                                                                            FieldVector<double,3>(1),
                                                                                                            elements3d);

  assert(hexahedralGrid->size(3) == numVertices3d);
  assert(hexahedralGrid->size(0) == numCubes3d);

  gridcheck(*hexahedralGrid);
#ifdef ModelP  // parallel UGGrid can only have one grid at a time
  hexahedralGrid.reset();
#endif
#else
  std::cout << "WARNING: 3d cube grids not tested because no suitable grid implementation is available!" << std::endl;
#endif

  // Test creation of 3d simplex grids
#if HAVE_UG
  typedef UGGrid<3> TetrahedralGridType;

  std::shared_ptr<TetrahedralGridType> tetrahedralGrid = StructuredGridFactory<TetrahedralGridType>::createSimplexGrid(FieldVector<double,3>(0),
                                                                                                                  FieldVector<double,3>(1),
                                                                                                                  elements3d);

  assert(tetrahedralGrid->size(3) == numVertices3d);
  assert(tetrahedralGrid->size(0) == 6*numCubes3d);    // each cube gets split into 6 tetrahedra

  gridcheck(*tetrahedralGrid);
#ifdef ModelP  // parallel UGGrid can only have one grid at a time
  tetrahedralGrid.reset();
#endif
#else
  std::cout << "WARNING: 3d simplicial grids not tested because no suitable grid implementation is available!" << std::endl;
#endif

  return 0;

}
catch (Exception &e) {
  std::cerr << e << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
