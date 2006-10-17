// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CHECK_INTERSECTIONITERATOR_CC
#define DUNE_CHECK_INTERSECTIONITERATOR_CC

#include <set>

#include <dune/grid/common/quadraturerules.hh>

/** \file
    \brief Tests for the IntersectionIterator
 */

/** \brief Helper routine: Test a general geometry
 */

template <class GeometryImp>
void checkGeometry(const GeometryImp& geometry)
{
  using namespace Dune;

  // Get the dimensions and types
  const int dim      = GeometryImp::mydimension;
  const int dimworld = GeometryImp::coorddimension;

  typedef typename GeometryImp::ctype ctype;

  // Get the corresponding reference element
  const ReferenceElement<double,dim>& refElement
    = ReferenceElements<double, dim>::general(geometry.type());

  // Check whether the number of corners is correct
  if (geometry.corners() != refElement.size(dim))
    DUNE_THROW(GridError, "Geometry has wrong number of corners!");

  // check consistency between operator[] and global()
  for (int i=0; i<refElement.size(dim); i++) {
    FieldVector<double,dim> localPos = refElement.position(i,dim);
    if ( (geometry[i] - geometry.global(localPos)).infinity_norm() > 1e-6)
      DUNE_THROW(GridError, "Methods operator[] and global() are inconsistent!");
  }

  // Compute the element center just to have an argument for the following methods
  FieldVector<ctype, dimworld> center(0);
  for (int i=0; i<geometry.corners(); i++)
    center += geometry[i];

  center /= geometry.corners();

#ifdef DUNE_UGGRID_HH
#warning Test for quadrilateral intersections disabled!
  if (!(geometry.type().isQuadrilateral() && dimworld==3)) {
#endif
  // The geometry center in local coordinates
  FieldVector<ctype, dim> localCenter = geometry.local(center);

  // Check whether center is within the intersection
  // This implicitly assumes convex intersections
  if (!geometry.checkInside(localCenter))
    DUNE_THROW(GridError, "Center of geometry is not within geometry!");

  // Back to global coordinates to check for correctness
  FieldVector<ctype, dimworld> worldCenter = geometry.global(localCenter);
  if ((center-worldCenter).infinity_norm() > 1e-6)
    DUNE_THROW(GridError, "local() and global() are not inverse to each other!");

  // The integration element at the element center
  ctype intElement = geometry.integrationElement(localCenter);
  if (intElement <=0)
    DUNE_THROW(GridError, "nonpositive integration element found!");

#if 0
  // This method exists in the interface, but it is not expected to work
  // unless dim==dimworld
  const FieldMatrix<ctype, Geometry::mydimension, Geometry::mydimension> jacobi
    = intersectionGlobal.jacobianInverseTransposed(localCenter);
#endif
#ifdef DUNE_UGGRID_HH
}
#endif


}

/** \brief Test the IntersectionIterator
 */
template <class GridType, class IndexSet, class ElementIterator, class IntersectionIterator>
void checkIntersectionIter(const GridType & grid, const IndexSet& indexSet,
                           const ElementIterator & eIt,
                           IntersectionIterator & iIt, const IntersectionIterator &iEndIt,
                           bool isConforming)
{
  using namespace Dune;

  typedef typename GridType::ctype ctype;

  for (;iIt!=iEndIt; ++iIt)
  {
    typedef typename IntersectionIterator::Entity EntityType;
    typedef typename EntityType::EntityPointer EntityPointer;

    assert(eIt == iIt.inside());

    // /////////////////////////////////////////////////////////////
    //   Check the consistency of numberInSelf, numberInNeighbor
    //   and the indices of the subface between.
    // /////////////////////////////////////////////////////////////
    if ( isConforming && iIt.neighbor() )
    {
      assert(iIt.level() >= 0);
      EntityPointer outside = iIt.outside();
      int numberInSelf     = iIt.numberInSelf();
      int numberInNeighbor = iIt.numberInNeighbor();

      assert(indexSet.template subIndex<1>(*eIt, numberInSelf)
             == indexSet.template subIndex<1>(*outside, numberInNeighbor));

#ifndef DUNE_UGGRID_HH
      assert(grid.localIdSet().template subId<1>(*eIt, numberInSelf)
             == grid.localIdSet().template subId<1>(*outside, numberInNeighbor));

      assert(grid.globalIdSet().template subId<1>(*eIt, numberInSelf)
             == grid.globalIdSet().template subId<1>(*outside, numberInNeighbor));
#else
#warning Test disabled, as UG does not have Face IDs
#endif
    }

    // //////////////////////////////////////////////////////////
    //   Check the geometry returned by intersectionGlobal()
    // //////////////////////////////////////////////////////////
    typedef typename IntersectionIterator::Geometry Geometry;
    const Geometry& intersectionGlobal = iIt.intersectionGlobal();

    checkGeometry(intersectionGlobal);

    // //////////////////////////////////////////////////////////
    //   Check the geometry returned by intersectionSelfLocal()
    // //////////////////////////////////////////////////////////

    const typename IntersectionIterator::LocalGeometry& intersectionSelfLocal = iIt.intersectionSelfLocal();
    checkGeometry(intersectionSelfLocal);

    //  Check the consistency of intersectionSelfLocal() and intersectionGlobal

    if (intersectionSelfLocal.corners() != intersectionGlobal.corners())
      DUNE_THROW(GridError, "Geometry of intersection is inconsistent from left hand side and global view!");

    // (Ab)use a quadrature rule as a set of test points
    const int interDim = IntersectionIterator::LocalGeometry::mydimension;
    const QuadratureRule<double, interDim>& quad
      = QuadratureRules<double, interDim>::rule(intersectionSelfLocal.type(), 2);

    for (size_t i=0; i<quad.size(); i++) {

      FieldVector<double,GridType::dimension> globalPos = intersectionGlobal.global(quad[i].position());
      FieldVector<double,GridType::dimension> localPos  = eIt->geometry().global(intersectionSelfLocal.global(quad[i].position()));

      if ( (globalPos - localPos).infinity_norm() > 1e-6)
        DUNE_THROW(GridError, "global( intersectionSelfLocal(global() ) is not the same as intersectionGlobal.global() at " << quad[i].position() << "!");

    }

    // ////////////////////////////////////////////////////////////////
    //   Check the geometry returned by intersectionNeighborLocal()
    // ////////////////////////////////////////////////////////////////

    if (iIt.neighbor() )
    {

      const typename IntersectionIterator::LocalGeometry& intersectionNeighborLocal = iIt.intersectionNeighborLocal();

      checkGeometry(intersectionNeighborLocal);

      if (intersectionSelfLocal.corners() != intersectionNeighborLocal.corners())
        DUNE_THROW(GridError, "Geometry of intersection is inconsistent from left and right hand side!");

      // (Ab)use a quadrature rule as a set of test points
      const int interDim = IntersectionIterator::LocalGeometry::mydimension;
      const QuadratureRule<double, interDim>& quad
        = QuadratureRules<double, interDim>::rule(intersectionNeighborLocal.type(), 2);

      for (size_t i=0; i<quad.size(); i++) {

        FieldVector<double,GridType::dimension> globalPos = intersectionGlobal.global(quad[i].position());
        FieldVector<double,GridType::dimension> localPos  = iIt.outside()->geometry().global(intersectionNeighborLocal.global(quad[i].position()));

        if ( (globalPos - localPos).infinity_norm() > 1e-6)
          DUNE_THROW(GridError, "global( intersectionNeighborLocal(global() ) is not the same as intersectionGlobal.global() at " << quad[i].position() << "!");

      }

    }

  }

}

/** \brief Test both IntersectionIterators
 */
template <class GridType>
void checkIntersectionIterator(const GridType& grid, bool skipLevelIntersectionTest = false) {

  using namespace Dune;

  typedef typename GridType::ctype ctype;

  // Loop over all levels
  if(skipLevelIntersectionTest)
  {
    std::cerr<<"WARNING: skip test of LevelIntersectionIterator! \n";
  }
  else
  {
    for (int i=0; i<=grid.maxLevel(); i++)
    {

      typedef typename GridType::template Codim<0>::LevelIterator ElementIterator;
      ElementIterator eIt    = grid.template lbegin<0>(i);
      ElementIterator eEndIt = grid.template lend<0>(i);

      for (; eIt!=eEndIt; ++eIt)
      {
        // check Level IntersectionIterator
        typedef typename GridType::template Codim<0>::Entity EntityType;
        typedef typename EntityType::LevelIntersectionIterator IntersectionIterator;

        IntersectionIterator iIt    = eIt->ilevelbegin();
        IntersectionIterator iEndIt = eIt->ilevelend();
        bool isConforming = Dune::Capabilities::isLevelwiseConforming < GridType > :: v ;

        // /////////////////////////////////////////////////////////
        //   Check the types defined by the iterator
        // /////////////////////////////////////////////////////////
        IsTrue< SameType<
                typename IntersectionIterator::ctype,
                typename GridType::ctype>::value == true >::yes();

        IsTrue<static_cast<int>(IntersectionIterator::dimension)
            == static_cast<int>(GridType::dimension)>::yes();

        IsTrue<static_cast<int>(IntersectionIterator::dimensionworld)
            == static_cast<int>(GridType::dimensionworld)>::yes();

        checkIntersectionIter(grid,grid.levelIndexSet(i),eIt,iIt,iEndIt,isConforming);
      }
    }
  }

  // test leaf intersection iterator
  {
    typedef typename GridType::template Codim<0>::LeafIterator ElementIterator;

    ElementIterator eEndIt = grid.template leafend<0>();
    for (ElementIterator eIt = grid.template leafbegin<0>(); eIt!=eEndIt; ++eIt)
    {
      typedef typename GridType::template Codim<0>::Entity EntityType;
      typedef typename EntityType::LeafIntersectionIterator IntersectionIterator;

      IntersectionIterator iIt    = eIt->ileafbegin();
      IntersectionIterator iEndIt = eIt->ileafend();
      bool isConforming = Dune :: Capabilities ::isLeafwiseConforming < GridType > :: v ;

      // /////////////////////////////////////////////////////////
      //   Check the types defined by the iterator
      // /////////////////////////////////////////////////////////
      IsTrue< SameType<
              typename IntersectionIterator::ctype,
              typename GridType::ctype>::value == true >::yes();

      IsTrue<static_cast<int>(IntersectionIterator::dimension)
          == static_cast<int>(GridType::dimension)>::yes();

      IsTrue<static_cast<int>(IntersectionIterator::dimensionworld)
          == static_cast<int>(GridType::dimensionworld)>::yes();

      checkIntersectionIter(grid,grid.leafIndexSet(),eIt,iIt,iEndIt,isConforming);
    }
  }
}

#endif
