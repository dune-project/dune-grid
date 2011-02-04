// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/** \file
    \brief A unit test for the BasicGeometry class
 */


#include <config.h>

#include <ostream>
#include <iostream>

#include <dune/common/exceptions.hh>

#include <dune/grid/genericgeometry/geometry.hh>



using namespace Dune;

void fail(int &result) {
  result = 1;
}
void pass(int &result) {
  if(result == 77) result = 0;
}

/** \brief Test the interface of the given BasicGeometry object */
template <class TestGeometry>
void testBasicGeometry(const TestGeometry& geometry, int &result)
{
  const int dim = TestGeometry::mydimension;

  geometry.normal(0, FieldVector<double,dim>(0));

  pass(result);
}

int main (int argc , char **argv) try
{
  // 77 means "SKIP"
  int result = 77;

  // Test a line segment
  // ...

  // Test a triangle
  // ...

  {   // Test a quadrilateral
    std::cout << "Testing quadrilaterals..." << std::endl;

    std::vector<FieldVector<double, 2> > corners(4);

    corners[0][0] = 0.5;   corners[0][1] = 0.0;
    corners[1][0] = 1.0;   corners[1][1] = 0.0;
    corners[2][0] = 0.5;   corners[2][1] = 0.25;
    corners[3][0] = 0.75;  corners[3][1] = 0.25;

    typedef GenericGeometry::BasicGeometry<2, GenericGeometry::DefaultGeometryTraits<double,2,2> > ElementGeometry;
    ElementGeometry insideGeometry( GenericGeometry::topologyId( GeometryType(GeometryType::cube,2 )), corners );

    testBasicGeometry(insideGeometry, result);
  }   // quadrilateral

  // Test a tetrahedron
  // ...

  // Test a pyramid
  // ...

  // Test a prism
  // ...

  // Test a hexahedron
  // ...

  return result;
}
catch (Dune::Exception& e) {
  std::cerr << e << std::endl;
  throw;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  throw;
}
