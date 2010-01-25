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
  FieldVector<double,1> pos;
  pos[0] = 0.6;    factory.insertVertex(pos);
  pos[0] = 1.0;    factory.insertVertex(pos);
  pos[0] = 0.2;    factory.insertVertex(pos);
  pos[0] = 0.0;    factory.insertVertex(pos);
  pos[0] = 0.4;    factory.insertVertex(pos);
  pos[0] = 0.3;    factory.insertVertex(pos);
  pos[0] = 0.7;    factory.insertVertex(pos);

  // Insert elements
  GeometryType segment(GeometryType::simplex,1);
  std::vector<unsigned int> v(2);
  v[0] = 6;  v[1] = 1;   factory.insertElement(segment, v);
  v[0] = 0;  v[1] = 6;   factory.insertElement(segment, v);
  v[0] = 4;  v[1] = 0;   factory.insertElement(segment, v);
  v[0] = 5;  v[1] = 4;   factory.insertElement(segment, v);
  v[0] = 2;  v[1] = 5;   factory.insertElement(segment, v);
  v[0] = 3;  v[1] = 2;   factory.insertElement(segment, v);

  // Create the grid
  return factory.createGrid();
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
