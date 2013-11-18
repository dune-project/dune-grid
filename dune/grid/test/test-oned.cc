// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// $Id$

#include <config.h>

#include <vector>
#include <memory>

#include <dune/grid/onedgrid.hh>

#include "gridcheck.cc"
#include "checkgeometryinfather.cc"
#include "checkintersectionit.cc"
#include "checkadaptation.cc"

using namespace Dune;

OneDGrid* testFactory()
{
  GridFactory<OneDGrid> factory;

  // Insert vertices
  std::vector<FieldVector<double,1> > vertexPositions(7);
  vertexPositions[0][0] = 0.6;    factory.insertVertex(vertexPositions[0]);
  vertexPositions[1][0] = 1.0;    factory.insertVertex(vertexPositions[1]);
  vertexPositions[2][0] = 0.2;    factory.insertVertex(vertexPositions[2]);
  vertexPositions[3][0] = 0.0;    factory.insertVertex(vertexPositions[3]);
  vertexPositions[4][0] = 0.4;    factory.insertVertex(vertexPositions[4]);
  vertexPositions[5][0] = 0.3;    factory.insertVertex(vertexPositions[5]);
  vertexPositions[6][0] = 0.7;    factory.insertVertex(vertexPositions[6]);

  // Insert elements
  GeometryType segment(GeometryType::simplex,1);
  std::vector<unsigned int> v(2);
  v[0] = 6;  v[1] = 1;   factory.insertElement(segment, v);
  v[0] = 4;  v[1] = 0;   factory.insertElement(segment, v);
  v[0] = 0;  v[1] = 6;   factory.insertElement(segment, v);
  v[0] = 5;  v[1] = 4;   factory.insertElement(segment, v);
  v[0] = 3;  v[1] = 2;   factory.insertElement(segment, v);
  v[0] = 2;  v[1] = 5;   factory.insertElement(segment, v);

  // Insert boundary segments
  std::vector<unsigned int> s(1);
  s[0] = 1;   factory.insertBoundarySegment(s);
  s[0] = 3;   factory.insertBoundarySegment(s);

  // Create the grid
  OneDGrid* grid = factory.createGrid();

  // //////////////////////////////////////////////////////////////
  //   Test whether the vertex numbering is in insertion order
  // //////////////////////////////////////////////////////////////
  OneDGrid::Codim<1>::LevelIterator vIt    = grid->lbegin<1>(0);
  OneDGrid::Codim<1>::LevelIterator vEndIt = grid->lend<1>(0);

  const OneDGrid::LevelGridView::IndexSet& levelIndexSet = grid->levelView(0).indexSet();
  const OneDGrid::LeafGridView::IndexSet&  leafIndexSet  = grid->leafGridView().indexSet();

  for (; vIt!=vEndIt; ++vIt) {
    unsigned int idx = levelIndexSet.index(*vIt);
    FieldVector<double,1> p = vIt->geometry().corner(0);
    if ( (vertexPositions[idx] - p).two_norm() > 1e-6 )
      DUNE_THROW(GridError, "Vertex with level index " << idx << " should have position " << vertexPositions[idx]
                                                       << " but has position " << p << ".");

    // leaf index should be the same
    if (idx != leafIndexSet.index(*vIt))
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

  OneDGrid::Codim<0>::LevelIterator eIt    = grid->lbegin<0>(0);
  OneDGrid::Codim<0>::LevelIterator eEndIt = grid->lend<0>(0);

  for (; eIt!=eEndIt; ++eIt) {
    unsigned int idx = levelIndexSet.index(*eIt);
    FieldVector<double,1> p = eIt->geometry().center();
    if ( (elementCenters[idx] - p).two_norm() > 1e-6 )
      DUNE_THROW(GridError, "Element with index " << idx << " should have center " << elementCenters[idx]
                                                  << " but has center " << p << ".");

    // leaf index should be the same
    if (idx != leafIndexSet.index(*eIt))
      DUNE_THROW(GridError, "Newly created OneDGrids should have matching level- and leaf element indices.");
  }

  // /////////////////////////////////////////////////////////////////////////
  //   Test whether the boundary segment numbering is in insertion order
  // /////////////////////////////////////////////////////////////////////////

  typedef  OneDGrid::LeafGridView GridView;
  typedef  GridView::IntersectionIterator IntersectionIterator;
  const GridView gridView = grid->leafGridView();

  for (eIt = grid->lbegin<0>(0); eIt!=eEndIt; ++eIt) {
    IntersectionIterator iIt = gridView.ibegin(*eIt);
    const IntersectionIterator iEndIt = gridView.iend(*eIt);

    for (; iIt!=iEndIt; ++iIt) {

      if (iIt->boundary()) {

        if (iIt->boundarySegmentIndex()==0 && std::abs(iIt->geometry().corner(0)[0] - 1) > 1e-6)
          DUNE_THROW(GridError, "BoundarySegment with index 0 should have position 1,"
                     << " but has " << iIt->geometry().corner(0) << ".");

        if (iIt->boundarySegmentIndex()==1 && std::abs(iIt->geometry().corner(0)[0]) > 1e-6)
          DUNE_THROW(GridError, "BoundarySegment with index 1 should have position 0,"
                     << " but has " << iIt->geometry().corner(0) << ".");

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
  grid.mark(1, * grid.leafbegin<0>());
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
  std::auto_ptr<Dune::OneDGrid> factoryGrid(testFactory());

  testOneDGrid(*factoryGrid.get());

  // Create a OneDGrid with an array of vertex coordinates and test it
  std::vector<double> coords(6);
  coords[0] = -1;
  coords[1] = -0.4;
  coords[2] = 0.1;
  coords[3] = 0.35;
  coords[4] = 0.38;
  coords[5] = 1;

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
