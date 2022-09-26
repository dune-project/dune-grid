// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>

#include <vector>
#include <memory>

#include <dune/grid/onedgrid.hh>

#include "gridcheck.hh"
#include "checkgeometryinfather.hh"
#include "checkintersectionit.hh"
#include "checkadaptation.hh"

using namespace Dune;

std::unique_ptr<OneDGrid> testFactory()
{
  GridFactory<OneDGrid> factory;

  // Insert vertices
  std::vector<FieldVector<double,1> > vertexPositions = {0.6, 1.0, 0.2, 0.0, 0.4, 0.3, 0.7};

  for (auto&& pos : vertexPositions)
    factory.insertVertex(pos);

  // Insert elements
  factory.insertElement(GeometryTypes::line, {6,1});
  factory.insertElement(GeometryTypes::line, {4,0});
  factory.insertElement(GeometryTypes::line, {0,6});
  factory.insertElement(GeometryTypes::line, {5,4});
  factory.insertElement(GeometryTypes::line, {3,2});
  factory.insertElement(GeometryTypes::line, {2,5});

  // Insert boundary segments
  factory.insertBoundarySegment({1});
  factory.insertBoundarySegment({3});

  // Create the grid
  std::unique_ptr<OneDGrid> grid(factory.createGrid());

  // //////////////////////////////////////////////////////////////
  //   Test whether the vertex numbering is in insertion order
  // //////////////////////////////////////////////////////////////

  const OneDGrid::LevelGridView::IndexSet& levelIndexSet = grid->levelGridView(0).indexSet();
  const OneDGrid::LeafGridView::IndexSet&  leafIndexSet  = grid->leafGridView().indexSet();

  for (const auto& vertex : vertices(grid->levelGridView(0)))
  {
    auto idx = levelIndexSet.index(vertex);
    FieldVector<double,1> p = vertex.geometry().corner(0);
    if ( (vertexPositions[idx] - p).two_norm() > 1e-6 )
      DUNE_THROW(GridError, "Vertex with level index " << idx << " should have position " << vertexPositions[idx]
                                                       << " but has position " << p << ".");

    // leaf index should be the same
    if (idx != leafIndexSet.index(vertex))
      DUNE_THROW(GridError, "Newly created OneDGrids should have matching level- and leaf vertex indices.");
  }

  // //////////////////////////////////////////////////////////////
  //   Test whether the element numbering is in insertion order
  // //////////////////////////////////////////////////////////////

  std::vector<FieldVector<double,1> > elementCenters(6);    // a priori knowledge: this is where the element centers should be
  elementCenters[0] = 0.85;
  elementCenters[1] = 0.5;
  elementCenters[2] = 0.65;
  elementCenters[3] = 0.35;
  elementCenters[4] = 0.1;
  elementCenters[5] = 0.25;

  for (const auto& element : elements(grid->levelGridView(0)))
  {
    auto idx = levelIndexSet.index(element);
    FieldVector<double,1> p = element.geometry().center();
    if ( (elementCenters[idx] - p).two_norm() > 1e-6 )
      DUNE_THROW(GridError, "Element with index " << idx << " should have center " << elementCenters[idx]
                                                  << " but has center " << p << ".");

    // leaf index should be the same
    if (idx != leafIndexSet.index(element))
      DUNE_THROW(GridError, "Newly created OneDGrids should have matching level- and leaf element indices.");
  }

  // /////////////////////////////////////////////////////////////////////////
  //   Test whether the boundary segment numbering is in insertion order
  // /////////////////////////////////////////////////////////////////////////

  typedef  OneDGrid::LeafGridView GridView;
  const GridView gridView = grid->leafGridView();

  for (const auto& element : elements(grid->levelGridView(0)))
  {
    for (const auto& intersection : intersections(gridView, element))
    {

      if (intersection.boundary()) {

        if (intersection.boundarySegmentIndex()==0 && std::abs(intersection.geometry().corner(0)[0] - 1) > 1e-6)
          DUNE_THROW(GridError, "BoundarySegment with index 0 should have position 1,"
                     << " but has " << intersection.geometry().corner(0) << ".");

        if (intersection.boundarySegmentIndex()==1 && std::abs(intersection.geometry().corner(0)[0]) > 1e-6)
          DUNE_THROW(GridError, "BoundarySegment with index 1 should have position 0,"
                     << " but has " << intersection.geometry().corner(0) << ".");

      }

    }

  }

  // //////////////////////////////////////////////////////////////
  //   return the grid for further tests
  // //////////////////////////////////////////////////////////////
  return grid;
}

void testOneDGrid(OneDGrid& grid)
{
  // check macro grid
  gridcheck(grid);

  // create hybrid grid
  grid.mark(1, * grid.leafGridView().begin<0>());
  grid.preAdapt();
  grid.adapt();
  grid.postAdapt();
  checkIntersectionIterator(grid);

  // check the grid again
  gridcheck(grid);

  grid.globalRefine(1);
  gridcheck(grid);

  // check geometry lifetime
  checkGeometryLifetime( grid.leafGridView() );

  // check the method geometryInFather()
  checkGeometryInFather(grid);

  // check the intersection iterator
  checkIntersectionIterator(grid);

  checkAdaptation( grid );
}

int main () try
{
  // Create a OneDGrid using the grid factory and test it
  std::unique_ptr<Dune::OneDGrid> factoryGrid(testFactory());

  testOneDGrid(*factoryGrid);

  // Create a OneDGrid with an array of vertex coordinates and test it
  std::vector<double> coords = {-1,
                                -0.4,
                                 0.1,
                                 0.35,
                                 0.38,
                                 1};

  Dune::OneDGrid coordsGrid(coords);

  testOneDGrid(coordsGrid);

  // Create a uniform OneDGrid and test it
  Dune::OneDGrid uniformGrid(7,       // Number of elements
                             -0.5,    // Left boundary
                             2.3      // Right boundary
                             );

  testOneDGrid(uniformGrid);

  // Test a uniform grid with RefinementType set to COPY
  Dune::OneDGrid uniformGrid2(7,       // Number of elements
                              -0.5,    // Left boundary
                              2.3      // Right boundary
                              );

  uniformGrid2.setRefinementType(OneDGrid::COPY);

  testOneDGrid(uniformGrid2);


  // everything okay
  return 0;

}
catch (Dune::Exception &e) {
  std::cerr << e << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
