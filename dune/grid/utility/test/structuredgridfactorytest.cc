// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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
#include <dune/common/test/testsuite.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/test/gridcheck.hh>

using namespace Dune;

Dune::TestSuite checkNEntities(const char* gridName, const char* entityName, int n, int nExpected)
{
  Dune::TestSuite t;
  t.check(nExpected == n)
    << gridName << ": has " << n << " " << entityName << ", but expected " << nExpected;
  return t;
}

int main (int argc , char **argv)
try {
  Dune::TestSuite t;

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

  t.subTest(checkNEntities("onedCubeGrid", "vertices", onedCubeGrid->size(1), static_cast<int>(elements1d[0]+1)));
  t.subTest(checkNEntities("onedCubeGrid", "elements", onedCubeGrid->size(0), static_cast<int>(elements1d[0])));

  gridcheck(*onedCubeGrid);


  // Test creation of 1d simplex grids
  std::shared_ptr<OneDGrid> onedSimplexGrid = StructuredGridFactory<OneDGrid>::createSimplexGrid(FieldVector<double,1>(0),
                                                                                            FieldVector<double,1>(1),
                                                                                            elements1d);

  t.subTest(checkNEntities("onedSimplexGrid", "vertices", onedSimplexGrid->size(1), static_cast<int>(elements1d[0]+1)));
  t.subTest(checkNEntities("onedSimplexGrid", "elements", onedSimplexGrid->size(0), static_cast<int>(elements1d[0])));

  gridcheck(*onedSimplexGrid);

  // Test creation of 1d YaspGrid
  std::shared_ptr<YaspGrid<1> > yaspGrid1d =
    StructuredGridFactory<YaspGrid<1> >::createCubeGrid(FieldVector<double,1>(0),
                                                        FieldVector<double,1>(1),
                                                        elements1d);

  t.subTest(checkNEntities("yaspGrid1d", "vertices", yaspGrid1d->size(1),static_cast<int>(elements1d[0]+1)));
  t.subTest(checkNEntities("yaspGrid1d", "elements", yaspGrid1d->size(1),static_cast<int>(elements1d[0]+1)));

  gridcheck(*yaspGrid1d);

  {
    std::shared_ptr<YaspGrid<1, EquidistantOffsetCoordinates<double,1> > > yaspGridOff1d
      = StructuredGridFactory<YaspGrid<1, EquidistantOffsetCoordinates<double,1> > >::createCubeGrid(FieldVector<double,1>(-1.),
                                                                FieldVector<double,1>(1),
                                                                elements1d);

    t.subTest(checkNEntities("yaspGridOff1d", "vertices", yaspGridOff1d->size(1),static_cast<int>(elements1d[0]+1)));
    t.subTest(checkNEntities("yaspGridOff1d", "elements", yaspGridOff1d->size(1),static_cast<int>(elements1d[0]+1)));

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

  t.subTest(checkNEntities("yaspGrid2d", "vertices", yaspGrid2d->size(2), static_cast<int>(numVertices2d)));
  t.subTest(checkNEntities("yaspGrid2d", "elements", yaspGrid2d->size(0), static_cast<int>(numCubes2d)));

  gridcheck(*yaspGrid2d);

  {
    std::shared_ptr<YaspGrid<2, EquidistantOffsetCoordinates<double,2> > > yaspGridOff2d
      = StructuredGridFactory<YaspGrid<2, EquidistantOffsetCoordinates<double,2> > >::createCubeGrid(FieldVector<double,2>(0),
                                                              FieldVector<double,2>(1),
                                                              elements2d);

    t.subTest(checkNEntities("yaspGridOff2d", "vertices", yaspGridOff2d->size(2), static_cast<int>(numVertices2d)));
    t.subTest(checkNEntities("yaspGridOff2d", "elements", yaspGridOff2d->size(0), static_cast<int>(numCubes2d)));

    gridcheck(*yaspGridOff2d);
  }

  // Test creation of 2d cube grid using UG
#if HAVE_DUNE_UGGRID
  typedef UGGrid<2> QuadrilateralGridType;

  std::shared_ptr<QuadrilateralGridType> quadrilateralGrid = StructuredGridFactory<QuadrilateralGridType>::createCubeGrid(FieldVector<double,2>(0),
                                                                                                                     FieldVector<double,2>(1),
                                                                                                                     elements2d);

  t.subTest(checkNEntities("quadrilateralGrid", "vertices", quadrilateralGrid->size(2), static_cast<int>(numVertices2d)));
  t.subTest(checkNEntities("quadrilateralGrid", "elements", quadrilateralGrid->size(0), static_cast<int>(numCubes2d)));

  gridcheck(*quadrilateralGrid);
#else
  std::cout << "WARNING: 2d cube grids not tested because no suitable grid implementation is available!" << std::endl;
#endif


  // Test creation of 2d triangle grid using UG
#if HAVE_DUNE_UGGRID
  typedef UGGrid<2> TriangleGridType;

  std::shared_ptr<TriangleGridType> triangleGrid = StructuredGridFactory<TriangleGridType>::createSimplexGrid(FieldVector<double,2>(0),
                                                                                                         FieldVector<double,2>(1),
                                                                                                         elements2d);

  t.subTest(checkNEntities("triangleGrid", "vertices", triangleGrid->size(2), static_cast<int>(numVertices2d)));
  // each cube gets split into 2 triangles
  t.subTest(checkNEntities("triangleGrid", "elements", triangleGrid->size(0), static_cast<int>(2*numCubes2d)));

  gridcheck(*triangleGrid);
#else
  std::cout << "WARNING: 2d simplicial grids not tested because no suitable grid implementation is available!" << std::endl;
#endif

  // /////////////////////////////////////////////////////////////////////////////
  //   Test 3d grids
  // /////////////////////////////////////////////////////////////////////////////

  // Test creation of 3d cube grids
#if HAVE_DUNE_UGGRID
  typedef UGGrid<3> HexahedralGridType;

  std::array<unsigned int,3> elements3d;
  elements3d.fill(4);
  unsigned int numVertices3d = (elements3d[0]+1) * (elements3d[1]+1) * (elements3d[2]+1);
  unsigned int numCubes3d    = elements3d[0] * elements3d[1] * elements3d[2];

  std::shared_ptr<HexahedralGridType> hexahedralGrid = StructuredGridFactory<HexahedralGridType>::createCubeGrid(FieldVector<double,3>(0),
                                                                                                            FieldVector<double,3>(1),
                                                                                                            elements3d);

  t.subTest(checkNEntities("hexahedralGrid", "vertices", hexahedralGrid->size(3), numVertices3d));
  t.subTest(checkNEntities("hexahedralGrid", "elements", hexahedralGrid->size(0), numCubes3d));

  gridcheck(*hexahedralGrid);
#else
  std::cout << "WARNING: 3d cube grids not tested because no suitable grid implementation is available!" << std::endl;
#endif

  // Test creation of 3d simplex grids
#if HAVE_DUNE_UGGRID
  typedef UGGrid<3> TetrahedralGridType;

  std::shared_ptr<TetrahedralGridType> tetrahedralGrid = StructuredGridFactory<TetrahedralGridType>::createSimplexGrid(FieldVector<double,3>(0),
                                                                                                                  FieldVector<double,3>(1),
                                                                                                                  elements3d);

  t.subTest(checkNEntities("tetrahedralGrid", "vertices", tetrahedralGrid->size(3), numVertices3d));
  // each cube gets split into 6 tetrahedra
  t.subTest(checkNEntities("tetrahedralGrid", "elements", tetrahedralGrid->size(0), 6*numCubes3d));

  gridcheck(*tetrahedralGrid);
#else
  std::cout << "WARNING: 3d simplicial grids not tested because no suitable grid implementation is available!" << std::endl;
#endif

  return t.exit();
}
catch (Exception &e) {
  std::cerr << e << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
