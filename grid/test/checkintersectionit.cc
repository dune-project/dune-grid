// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CHECK_INTERSECTIONITERATOR_CC
#define DUNE_CHECK_INTERSECTIONITERATOR_CC

#include <cmath>

#include <dune/grid/common/quadraturerules.hh>

/** \file
    \brief Tests for the IntersectionIterator
 */

struct CheckIntersectionIteratorErrorState
{
  unsigned int sumNormalsNonZero;

  CheckIntersectionIteratorErrorState ()
    : sumNormalsNonZero( 0 )
  {}
};



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

  // Use a quadrature rule to create a few test points for the following checks
  const QuadratureRule<double, dim>& quad
    = QuadratureRules<double, dim>::rule(geometry.type(), 2);

  for (size_t i=0; i<quad.size(); i++) {

    const FieldVector<double,dim>& testPoint = quad[i].position();

    // Check whether point is within the intersection
    if (!geometry.checkInside(testPoint))
    {
      std :: cerr << "Test point (" << testPoint << ") not within geometry."
                  << std :: endl;
      //DUNE_THROW(GridError, "Test point is not within geometry!");
    }

    // Transform to global coordinates
    FieldVector<ctype, dimworld> global = geometry.global(testPoint);

    // The back to local coordinates
    FieldVector<ctype, dim> local = geometry.local(global);

    // check for correctness
    if ((testPoint-local).infinity_norm() > 1e-6)
      DUNE_THROW(GridError, "local() and global() are not inverse to each other!");

    // The integration element at the element center
    ctype intElement = geometry.integrationElement(testPoint);
    if (intElement <=0)
      DUNE_THROW(GridError, "nonpositive integration element found!");

#if 0
    // This method exists in the interface, but it is not expected to work
    // unless dim==dimworld
    const FieldMatrix<ctype, dim, dim> jacobi
      = intersectionGlobal.jacobianInverseTransposed(testPoint);
#endif
  }
}


// Check whether the normal is orthogonal to the intersection, i.e.,
// whether (J^-1 * n) = 0. Here J is the jacobian of the intersection
// geometry (intersectionGlobal) and n is a normal.
template< class ctype, int dimworld, int facedim >
inline bool checkNormal ( const Dune :: FieldVector< ctype, dimworld > &normal,
                          const Dune :: FieldMatrix< ctype, dimworld, facedim > &jit )
{
  Dune :: FieldVector< ctype, facedim > x( ctype( 0 ) );
  jit.umtv( normal, x );
  return (x.infinity_norm() <= 1e-8);
}


/** \brief Test the IntersectionIterator
 */
template <class GridViewType, class ErrorState >
void checkIntersectionIterator(const GridViewType& view,
                               const typename GridViewType::template Codim<0>::Iterator& eIt,
                               ErrorState &errorState )
{
  using namespace Dune;

  typedef typename GridViewType::Grid GridType;
  typedef typename GridViewType::IntersectionIterator IntersectionIterator;

  const GridType& grid = view.grid();
  const bool checkOutside = (grid.name() != "AlbertaGrid");
  const typename GridViewType::IndexSet& indexSet = view.indexSet();

  typedef typename GridType :: ctype ctype;

  const int dim      = GridType :: dimension;
  const int dimworld = GridType :: dimensionworld;

  FieldVector< ctype,dimworld > sumNormal( ctype( 0 ) );

  // /////////////////////////////////////////////////////////
  //   Check the types defined by the iterator
  // /////////////////////////////////////////////////////////
  dune_static_assert((is_same<
                        typename IntersectionIterator::ctype,
                        typename GridType::ctype>::value),
                     "IntersectionIterator has wrong ctype");

  dune_static_assert((is_same<
                        typename IntersectionIterator::Intersection,
                        typename GridViewType::Intersection>::value),
                     "IntersectionIterator has wrong Intersection type");

  dune_static_assert((static_cast<int>(IntersectionIterator::dimension)
                      == static_cast<int>(GridType::dimension)),"IntersectionIterator has wrong dimension");

  dune_static_assert((static_cast<int>(IntersectionIterator::dimensionworld)
                      == static_cast<int>(GridType::dimensionworld)),"IntersectionIterator has wrong dimensionworld");

  IntersectionIterator iIt    = view.ibegin(*eIt);
  IntersectionIterator iEndIt = view.iend(*eIt);

  bool hasBoundaryIntersection = false;

  for (;iIt!=iEndIt; ++iIt)
  {
    // //////////////////////////////////////////////////////////////////////
    //   Compute the integral of the outer normal over the whole element.
    //   This has to be zero.
    // //////////////////////////////////////////////////////////////////////
    const int interDim = IntersectionIterator::LocalGeometry::mydimension;
    const QuadratureRule<double, interDim>& quad
      = QuadratureRules<double, interDim>::rule(iIt->intersectionSelfLocal().type(), interDim);

    typedef typename IntersectionIterator::Entity EntityType;
    typedef typename EntityType::EntityPointer EntityPointer;

    assert(eIt == iIt->inside());

    // check that boundary id has positive value
    if( iIt->boundary() )
    {
      // entity has boundary intersections
      hasBoundaryIntersection = true;

      if( iIt->boundaryId() <= 0 )
      {
        DUNE_THROW(GridError, "boundary id has non-positive value (" << iIt->boundaryId() << ") !");
      }
    }

    // //////////////////////////////////////////////////////////////////////
    //   Check whether the 'has-intersection-with'-relation is symmetric
    // //////////////////////////////////////////////////////////////////////

    if (iIt->neighbor() && checkOutside )
    {
      EntityPointer outside = iIt->outside();
      bool insideFound = false;

      IntersectionIterator outsideIIt    = view.ibegin(*outside);
      IntersectionIterator outsideIEndIt = view.iend(*outside);

      for (; outsideIIt!=outsideIEndIt; ++outsideIIt) {

        if (outsideIIt->neighbor() && outsideIIt->outside() == iIt->inside()) {

          if (outsideIIt->numberInSelf() != iIt->numberInNeighbor())
            DUNE_THROW(GridError, "outside()->outside() == inside(), but with incorrect numbering!");
          else
            insideFound = true;

        }

      }

      if (!insideFound)
        DUNE_THROW(GridError, "Could not find inside() through intersection iterator of outside()!");

    }
    else if (!checkOutside)
    {
      static bool called = false;
      if(!called)
      {
        derr << "WARNING: skip reverse intersection iterator test for " << grid.name() << "!"<< std::endl;
        called = true;
      }
    }

    // /////////////////////////////////////////////////////////////
    //   Check the consistency of numberInSelf, numberInNeighbor
    //   and the indices of the subface between.
    // /////////////////////////////////////////////////////////////
    if ( GridViewType::conforming && iIt->neighbor() )
    {
      EntityPointer outside = iIt->outside();
      const int numberInSelf     = iIt->numberInSelf();
      const int numberInNeighbor = iIt->numberInNeighbor();

      const unsigned int iIdx = indexSet.template subIndex< 1 >( *eIt, numberInSelf );
      const unsigned int oIdx = indexSet.template subIndex< 1 >( *outside, numberInNeighbor );

      if( iIdx != oIdx )
      {
        std :: cerr << "Error: index of conforming intersection differs when "
                    << "obtained from inside and outside." << std :: endl;
        std :: cerr << "       inside index = " << iIdx
                    << ", outside index = " << oIdx << std :: endl;
        assert( false );
      }

      assert(grid.localIdSet().template subId<1>(*eIt, numberInSelf)
             == grid.localIdSet().template subId<1>(*outside, numberInNeighbor));

      assert(grid.globalIdSet().template subId<1>(*eIt, numberInSelf)
             == grid.globalIdSet().template subId<1>(*outside, numberInNeighbor));
    }

    // //////////////////////////////////////////////////////////
    //   Check the geometry returned by intersectionGlobal()
    // //////////////////////////////////////////////////////////
    typedef typename IntersectionIterator::Geometry Geometry;
    const Geometry& intersectionGlobal = iIt->intersectionGlobal();

    checkGeometry(intersectionGlobal);

    // //////////////////////////////////////////////////////////
    //   Check the geometry returned by intersectionSelfLocal()
    // //////////////////////////////////////////////////////////

    const typename IntersectionIterator::LocalGeometry& intersectionSelfLocal = iIt->intersectionSelfLocal();
    checkGeometry(intersectionSelfLocal);

    //  Check the consistency of intersectionSelfLocal() and intersectionGlobal

    if (intersectionSelfLocal.corners() != intersectionGlobal.corners())
      DUNE_THROW(GridError, "Geometry of intersection is inconsistent from left hand side and global view!");

    // Use a quadrature rule as a set of test points
    for (size_t i=0; i<quad.size(); i++)
    {
      const FieldVector< ctype, dim-1 > &pt = quad[ i ].position();
      const FieldMatrix< ctype, dimworld, dim-1 > &jit
        = intersectionGlobal.jacobianInverseTransposed( pt );

      // Check outer normal
      const FieldVector< ctype, dimworld > normal = iIt->outerNormal( pt );

      if( !checkNormal( normal, jit ) )
      {
        std :: cerr << "Error: outerNormal "
                    << "is not orthogonal to the intersection." << std :: endl;
        std :: cerr << "       outerNormal = " << normal << std :: endl;
        assert( false );
      }

      // Check integration outer normal
      const FieldVector< ctype, dimworld > intNormal = iIt->integrationOuterNormal( pt );
      sumNormal.axpy( quad[ i ].weight(), intNormal );

      const ctype det = intersectionGlobal.integrationElement( pt );
      if( std :: abs( det - intNormal.two_norm() ) > 1e-8 )
      {
        std :: cerr << "Error: integrationOuterNormal yields wrong length."
                    << std :: endl;
        std :: cerr << "       |integrationOuterNormal| = " << intNormal.two_norm()
                    << ", integrationElement = " << det << std :: endl;
        assert( false );
      }

      if( !checkNormal( intNormal, jit ) )
      {
        std :: cerr << "Error: integrationOuterNormal "
                    << "is not orthogonal to the intersection." << std :: endl;
        std :: cerr << "       integrationOuterNormal = " << intNormal
                    << std :: endl;
        assert( false );
      }

      if( std :: abs( normal * intNormal - normal.two_norm() * intNormal.two_norm() ) > 1e-8 )
      {
        std :: cerr << "Error: integrationOuterNormal "
                    << "does not point in the direction of outerNormal." << std :: endl;
        std :: cerr << "       integrationOuterNormal = " << intNormal
                    << ", normal = " << normal << std :: endl;
        assert( false );
      }

      // Check unit outer normal
      const FieldVector< ctype, dimworld > unitNormal = iIt->unitOuterNormal( pt );

      if( std :: abs( ctype( 1 ) - unitNormal.two_norm() ) > 1e-8 )
      {
        std :: cerr << "Error: unitOuterNormal yields wrong length." << std :: endl;
        std :: cerr << "       |unitOuterNormal| = " << unitNormal.two_norm()
                    << std :: endl;
        assert( false );
      }

      if( !checkNormal( unitNormal, jit ) )
      {
        std :: cerr << "Error: unitOuterNormal "
                    << "is not orthogonal to the intersection." << std :: endl;
        std :: cerr << "       unitOuterNormal = " << unitNormal
                    << std :: endl;
        assert( false );
      }

      if( std :: abs( normal * unitNormal - normal.two_norm() * unitNormal.two_norm() ) > 1e-8 )
      {
        std :: cerr << "Error: unitOuterNormal "
                    << "does not point in the direction of outerNormal." << std :: endl;
        std :: cerr << "       unitOuterNormal = " << unitNormal
                    << ", normal = " << normal << std :: endl;
        assert( false );
      }

      // check intersectionSelfLocal
      FieldVector<double,dimworld> globalPos = intersectionGlobal.global(quad[i].position());
      FieldVector<double,dimworld> localPos  = eIt->geometry().global(intersectionSelfLocal.global(quad[i].position()));

      if ( (globalPos - localPos).infinity_norm() > 1e-6)
        DUNE_THROW(GridError, "global( intersectionSelfLocal(global() ) is not the same as intersectionGlobal.global() at " << quad[i].position() << "!");

    }

    // ////////////////////////////////////////////////////////////////
    //   Check the geometry returned by intersectionNeighborLocal()
    // ////////////////////////////////////////////////////////////////

    if (iIt->neighbor() )
    {

      const typename IntersectionIterator::LocalGeometry& intersectionNeighborLocal = iIt->intersectionNeighborLocal();

      checkGeometry(intersectionNeighborLocal);

      if (intersectionSelfLocal.corners() != intersectionNeighborLocal.corners())
        DUNE_THROW(GridError, "Geometry of intersection is inconsistent from left and right hand side!");

      // (Ab)use a quadrature rule as a set of test points
      const int interDim = IntersectionIterator::LocalGeometry::mydimension;
      const QuadratureRule<double, interDim>& quad
        = QuadratureRules<double, interDim>::rule(intersectionNeighborLocal.type(), 2);

      for (size_t i=0; i<quad.size(); i++)
      {

        FieldVector<double,dimworld> globalPos = intersectionGlobal.global(quad[i].position());
        FieldVector<double,dimworld> localPos  = iIt->outside()->geometry().global(intersectionNeighborLocal.global(quad[i].position()));

        if ( (globalPos - localPos).infinity_norm() > 1e-6)
          DUNE_THROW(GridError, "global( intersectionNeighborLocal(global() ) is not the same as intersectionGlobal.global() at " << quad[i].position() << "!");

      }

    }

  }

  // check implementation of hasBoundaryIntersections
  if( hasBoundaryIntersection != eIt->hasBoundaryIntersections() )
  {
    DUNE_THROW(GridError,"Entity::hasBoundaryIntersections implemented wrong! \n");
  }

  // Chech whether integral over the outer normal is zero
  //
  // Note: This is wrong on curved surfaces (take, e.g., the upper half sphere).
  //       Therefore we only issue a warning.
  if( sumNormal.two_norm() > 1e-8 )
  {
    std :: cout << "Warning: Integral over outer normals is nonzero: "
                << sumNormal << std :: endl;
    ++errorState.sumNormalsNonZero;
  }
}

/** \brief Test both IntersectionIterators
 */
template <class GridViewType>
void checkViewIntersectionIterator(const GridViewType& view) {

  using namespace Dune;

  typedef typename GridViewType ::template Codim<0>::Iterator
  ElementIterator;
  ElementIterator eIt    = view.template begin<0>();
  ElementIterator eEndIt = view.template end<0>();

  CheckIntersectionIteratorErrorState errorState;
  for (; eIt!=eEndIt; ++eIt)
    checkIntersectionIterator( view, eIt, errorState );

  if( errorState.sumNormalsNonZero > 0 )
  {
    std :: cerr << "Warning: Integral over outer normals is not always zero."
                << std :: endl;
    std :: cerr << "         This behaviour may be correct for entities with"
                << " nonzero curature." << std :: endl;;
  }
}

template <class GridType>
void checkIntersectionIterator(const GridType& grid, bool skipLevelIntersectionTest = false) {

  using namespace Dune;

  // Loop over all levels
  if(skipLevelIntersectionTest)
  {
    std::cerr<<"WARNING: skip test of LevelIntersectionIterator! \n";
  }
  else
  {
    for (int i=0; i<=grid.maxLevel(); i++)
      checkViewIntersectionIterator(grid.levelView(i));
  }

  // test leaf intersection iterator
  {
    checkViewIntersectionIterator(grid.leafView());
  }

}

#endif
