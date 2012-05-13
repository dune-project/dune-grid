// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \file
 * \brief Unit tests for the virtual refinement code
 */

#include "config.h"

#include <iostream>
#include <ostream>

#include <dune/geometry/genericgeometry/topologytypes.hh>
#include <dune/geometry/type.hh>

#include <dune/grid/common/virtualrefinement.hh>
#include <dune/grid/test/checkgeometry.cc>

using namespace Dune;

void pass(int &result) {
  // 77 means 'SKIP'
  if(result == 77) result = 0;
}
void fail(int &result) {
  result = 1;
}
void collect(int &result, bool passed)
{
  if(passed) pass(result);else fail(result);
}

/** \brief Test virtual refinement for an element with a run-time type
 */
template <class ct, int dim>
void testVirtualRefinement(int &result, const Dune::GeometryType& elementType,
                           const Dune::GeometryType& coerceTo, int refinement)
{
  std::cout << "Checking virtual refinement " << elementType << " -> "
            << coerceTo << " level " << refinement << std::endl;

  typedef Dune::VirtualRefinement<dim, ct> Refinement;
  typedef typename Refinement::ElementIterator eIterator;

  // Make a virtual refinement of the reference element
  Refinement & elementRefinement =
    Dune::buildRefinement<dim, ct>(elementType,coerceTo);

  eIterator eSubEnd = elementRefinement.eEnd(refinement);
  eIterator eSubIt  = elementRefinement.eBegin(refinement);

  for (; eSubIt != eSubEnd; ++eSubIt)
    std::cout << eSubIt.coords() << std::endl;

  pass(result);
}


template <unsigned topologyId, class ct, unsigned coerceToId, int dim>
void testStaticRefinementGeometry(int &result, int refinement)
{
  std::cout << "Checking static refinement geometry "
            << GeometryType(topologyId, dim) << " -> "
            << GeometryType(coerceToId, dim) << " level " << refinement
            << std::endl;

  typedef Dune::StaticRefinement<topologyId, ct, coerceToId, dim> Refinement;
  typedef typename Refinement::ElementIterator eIterator;

  eIterator eSubEnd = Refinement::eEnd(refinement);
  eIterator eSubIt  = Refinement::eBegin(refinement);

  for (; eSubIt != eSubEnd; ++eSubIt)
    // Call the standard test for geometries
    collect(result, checkGeometry(eSubIt.geometry()));

}



int main(int argc, char** argv) try
{
  using GenericGeometry::Point;

  typedef GenericGeometry::Prism  <Point>    Line;

  typedef GenericGeometry::Prism  <Line>     Square;
  typedef GenericGeometry::Pyramid<Line>     Triangle;

  typedef GenericGeometry::Prism  <Square>   Cube;
  typedef GenericGeometry::Pyramid<Square>   Pyramid;
  typedef GenericGeometry::Prism  <Triangle> Prism;
  typedef GenericGeometry::Pyramid<Triangle> Tet;

  // 77 means 'SKIP'
  int result = 77;

  GeometryType gt1, gt2;

  // test segment
  gt1.makeLine();
  gt2.makeLine();
  for (int refinement=0; refinement<3; refinement++) {
    testVirtualRefinement<double,1>(result, gt1, gt2, refinement);
    testStaticRefinementGeometry<Line::id,double,Line::id,1>
      (result, refinement);
  }

  // test triangle
  gt1.makeTriangle();
  gt2.makeTriangle();
  for (int refinement=0; refinement<3; refinement++) {
    testVirtualRefinement<double,2>(result, gt1, gt2, refinement);
    testStaticRefinementGeometry<Triangle::id,double,Triangle::id,2>
      (result, refinement);
  }

  // test quadrilateral
  gt1.makeQuadrilateral();
  gt2.makeQuadrilateral();
  for (int refinement=0; refinement<3; refinement++) {
    testVirtualRefinement<double,2>(result, gt1, gt2, refinement);
    // testStaticRefinementGeometry<Square::id,double,Square::id,2>
    //     (result, refinement);
  }
  gt2.makeTriangle();
  for (int refinement=0; refinement<3; refinement++) {
    testVirtualRefinement<double,2>(result, gt1, gt2, refinement);
    testStaticRefinementGeometry<Square::id,double,Triangle::id,2>
      (result, refinement);
  }

  // test tetrahedron
  gt1.makeTetrahedron();
  gt2.makeTetrahedron();
  for (int refinement=0; refinement<3; refinement++) {
    testVirtualRefinement<double,3>(result, gt1, gt2, refinement);
    testStaticRefinementGeometry<Tet::id,double,Tet::id,3>
      (result, refinement);
  }

  // test pyramid
  gt1.makePyramid();
  gt2.makeTetrahedron();
  for (int refinement=0; refinement<3; refinement++) {
    testVirtualRefinement<double,3>(result, gt1, gt2, refinement);
    // testStaticRefinementGeometry<Pyramid::id,double,Tet::id,3>
    //     (result, refinement);
  }

  // test prism
  gt1.makePrism();
  gt2.makeTetrahedron();
  for (int refinement=0; refinement<3; refinement++) {
    testVirtualRefinement<double,3>(result, gt1, gt2, refinement);
    // testStaticRefinementGeometry<Prism::id,double,Tet::id,3>
    //     (result, refinement);
  }

  // test hexahedron
  gt1.makeHexahedron();
  gt2.makeHexahedron();
  for (int refinement=0; refinement<3; refinement++) {
    testVirtualRefinement<double,3>(result, gt1, gt2, refinement);
    // testStaticRefinementGeometry<Cube::id,double,Cube::id,3>
    //    (result, refinement);
  }
  gt1.makeHexahedron();
  gt2.makeTetrahedron();
  for (int refinement=0; refinement<3; refinement++) {
    testVirtualRefinement<double,3>(result, gt1, gt2, refinement);
    testStaticRefinementGeometry<Cube::id,double,Tet::id,3>
      (result, refinement);
  }

  return result;

}
catch (Exception &e) {
  std::cerr << e << std::endl;
  throw;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  throw;
}
