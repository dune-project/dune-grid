// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/** \file
 * \brief Unit tests for the virtual refinement code
 */

#include "config.h"

#include <dune/geometry/type.hh>
#include <dune/grid/common/virtualrefinement.hh>
#include <dune/grid/test/checkgeometry.cc>

using namespace Dune;

/** \brief Test virtual refinement for an element with a run-time type
 */
template <class ct, int dim>
void testVirtualRefinement(const Dune::GeometryType& elementType, int refinement)
{
  typedef Dune::VirtualRefinement<dim, ct> Refinement;
  typedef typename Refinement::ElementIterator eIterator;

  // Make a virtual refinement of the reference element
  Refinement & elementRefinement = Dune::buildRefinement<dim, ct>(elementType,elementType);

  eIterator eSubEnd = elementRefinement.eEnd(refinement);
  eIterator eSubIt  = elementRefinement.eBegin(refinement);

  for (; eSubIt != eSubEnd; ++eSubIt)
    std::cout << eSubIt.coords() << std::endl;

}


template <class ct, int dim>
void testSimplexRefinement(int refinement)
{
  // Create a simplex-to-simplex refinement
  typedef Dune::StaticRefinement<Dune::GenericGeometry::SimplexTopology<dim>::type::id,
      ct,
      Dune::GenericGeometry::SimplexTopology<dim>::type::id,
      dim> Refinement;
  typedef typename Refinement::ElementIterator eIterator;

  eIterator eSubEnd = Refinement::eEnd(refinement);
  eIterator eSubIt  = Refinement::eBegin(refinement);

  for (; eSubIt != eSubEnd; ++eSubIt)
    // Call the standard test for geometries
    checkGeometry(eSubIt.geometry());

}



int main(int argc, char** argv) try
{
  GeometryType simplex;

  // test segment
  simplex.makeSimplex(1);
  for (int refinement=0; refinement<3; refinement++) {
    testVirtualRefinement<double,1>(simplex,refinement);
    testSimplexRefinement<double,1>(refinement);
  }

  // test triangle
  simplex.makeSimplex(2);
  for (int refinement=0; refinement<3; refinement++) {
    testVirtualRefinement<double,2>(simplex,refinement);
    testSimplexRefinement<double,2>(refinement);
  }

  // test tetrahedron
  simplex.makeSimplex(3);
  for (int refinement=0; refinement<3; refinement++) {
    testVirtualRefinement<double,3>(simplex,refinement);
    testSimplexRefinement<double,3>(refinement);
  }

  return 0;

}
catch (Exception &e) {
  std::cerr << e << std::endl;
  return 1;
} catch (...) {
  std::cerr << "Generic exception!" << std::endl;
  return 2;
}
