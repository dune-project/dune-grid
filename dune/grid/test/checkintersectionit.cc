// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_CHECK_INTERSECTIONITERATOR_CC
#define DUNE_CHECK_INTERSECTIONITERATOR_CC

#include <cmath>

#include <dune/geometry/quadraturerules/gaussquadrature.hh>
#include <dune/geometry/genericreferenceelements.hh>

#include "checkgeometry.cc"

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


template< class Grid >
struct EnableIntersectionIteratorReverseCheck
{
  static const bool v = true;
};


// Check that normal and normal2 pointing in the same direction
template< class ctype, int dimworld, class String >
inline void checkParallel ( const Dune::FieldVector< ctype, dimworld > &normal,
                            const Dune::FieldVector< ctype, dimworld > &refNormal,
                            const String & name )
{
  if( (normal.two_norm()*refNormal.two_norm() - normal*refNormal) > 1e-8 )
  {
    std::cerr << "Error: " << name << " does not point in the direction of outer normal." << std::endl;
    std::cerr << "       " << name << " = " << normal << ", outer normal = " << refNormal << std :: endl;
    assert( false );
  }
}


// Check whether the normal is orthogonal to the intersection, i.e.,
// whether (J^-1 * n) = 0. Here J is the jacobian of the intersection
// geometry (geometry) and n is a normal.
template< class ctype, int dimworld, class JacobianInverseTransposed, class String >
inline void checkJIn ( const Dune :: FieldVector< ctype, dimworld > &normal,
                       const JacobianInverseTransposed &jit,
                       const String & name )
{
  const int facedim = JacobianInverseTransposed::cols;
  Dune :: FieldVector< ctype, facedim > x( ctype( 0 ) );
  jit.umtv( normal, x );
  if (x.infinity_norm() > 1e-8)
  {
    const Dune::FieldMatrix< ctype, dimworld, facedim > &mjit = jit;
    std :: cerr << "Error:  (J^-1 * n) != 0." << std :: endl;
    std :: cerr << "       " << name << " = " << normal
                << std :: endl;
    std :: cerr << "       J^-1^T = \n" << mjit << std::endl;
    assert( false );
  }
}


template< class GridView >
void checkIntersectionAssignment ( const GridView &view, const typename GridView::template Codim< 0 >::Iterator &l )
{
  typedef typename GridView::template Codim< 0 >::Iterator Iterator;
  typedef typename GridView::IntersectionIterator IntersectionIterator;

  Iterator l1( l );
  Iterator l2( l );
  ++l2;

  if( l2 != view.template end< 0 >() )
  {
    IntersectionIterator it1 = view.ibegin( *l1 );
    IntersectionIterator it2 = view.ibegin( *l2 );

    // verify prerequisites
    assert( l1 != l2 );
    assert( it1 != it2 );
    assert( it1->inside() == l1 );
    assert( it2->inside() == l2 );

    // check assignment
    it1 = it2;
    assert( it1 == it2 );
    assert( it1->inside() == l2 );
  }
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
  typedef typename IntersectionIterator::Intersection Intersection;

  const GridType& grid = view.grid();
  const bool checkOutside = EnableIntersectionIteratorReverseCheck< GridType >::v;
  const typename GridViewType::IndexSet& indexSet = view.indexSet();

  typedef typename GridType::ctype ctype;

  const int dim      = GridType::dimension;
  // const int dimworld = GridType::dimensionworld;

  typedef typename GridViewType::template Codim< 0 >::Geometry ElementGeometry;
  const ElementGeometry &geoInside = eIt->geometry();
  const GenericReferenceElement< ctype, dim > &refElement = GenericReferenceElements< ctype, dim >::general( eIt->type() );

  typename ElementGeometry::GlobalCoordinate sumNormal( ctype( 0 ) );

  // /////////////////////////////////////////////////////////
  //   Check the types defined by the iterator
  // /////////////////////////////////////////////////////////
  dune_static_assert((is_same<
                        typename Intersection::ctype,
                        typename GridType::ctype>::value),
                     "IntersectionIterator has wrong ctype");

  dune_static_assert((is_same<
                        typename IntersectionIterator::Intersection,
                        typename GridViewType::Intersection>::value),
                     "IntersectionIterator has wrong Intersection type");

  dune_static_assert((static_cast<int>(Intersection::dimension)
                      == static_cast<int>(GridType::dimension)),"IntersectionIterator has wrong dimension");

  dune_static_assert((static_cast<int>(Intersection::dimensionworld)
                      == static_cast<int>(GridType::dimensionworld)),"IntersectionIterator has wrong dimensionworld");

  IntersectionIterator iIt    = view.ibegin(*eIt);
  IntersectionIterator iEndIt = view.iend(*eIt);

  bool hasBoundaryIntersection = false;

  for (;iIt!=iEndIt; ++iIt)
  {
    const int indexInInside  = iIt->indexInInside();

    const GenericReferenceElement< ctype, dim-1 > &refFace = GenericReferenceElements< ctype, dim-1 >::general( iIt->type() );

    // //////////////////////////////////////////////////////////////////////
    //   Compute the integral of the outer normal over the whole element.
    //   This has to be zero.
    // //////////////////////////////////////////////////////////////////////
    const int interDim = Intersection::LocalGeometry::mydimension;

    typedef typename Intersection::Entity EntityType;
    typedef typename EntityType::EntityPointer EntityPointer;

    assert(eIt == iIt->inside());

    // check that boundary id has positive value and that intersection is
    // conform
    if( iIt->boundary() )
    {
      // entity has boundary intersections
      hasBoundaryIntersection = true;

#if !DISABLE_DEPRECATED_METHOD_CHECK
      if( iIt->boundaryId() < 0 )
      {
        DUNE_THROW(GridError, "boundary id has negative value (" << iIt->boundaryId() << ") !");
      }
#endif // #if !DISABLE_DEPRECATED_METHOD_CHECK
      if( ! iIt->conforming() && ! iIt->neighbor() )
      {
        DUNE_THROW(GridError, "Boundary intersection should be conforming!");
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

        if (outsideIIt->neighbor() && outsideIIt->outside() == iIt->inside())
        {
          insideFound = true;
          const int idxInside = outsideIIt->indexInInside();
          const int idxOutside = iIt->indexInOutside();
          if( idxOutside != idxInside )
          {
            std::cerr << "Error: outside()->outside() == inside()"
                      << ", but with incorrect numbering" << std::endl;
            std::cerr << "       (inside().indexInOutside = " << idxOutside
                      << ", outside().indexInInside = " << idxInside << ")." << std::endl;
            assert( false );
          }
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
        derr << "Warning: Skipping reverse intersection iterator check." << std::endl;
        called = true;
      }
    }

    // Check if conforming() methods is compatible with static
    // information on GridView
    if ( GridViewType::conforming && !iIt->conforming())
    {
      DUNE_THROW(GridError, "GridView says conforming but intersection is not conforming!");
    }

    // /////////////////////////////////////////////////////////////
    //   Check the consistency of numberInSelf, numberInNeighbor
    //   and the indices of the subface between.
    // /////////////////////////////////////////////////////////////
    if( iIt->conforming() && iIt->neighbor() && !iIt->boundary() )
    {
      EntityPointer outside = iIt->outside();
      const int indexInOutside = iIt->indexInOutside();

      if( indexSet.subIndex( *eIt, indexInInside, 1 ) != indexSet.subIndex( *outside, indexInOutside, 1 ) )
      {
        std::cerr << "Error: Index of conforming intersection differs when "
                  << "obtained from inside and outside." << std::endl;
        std::cerr << "       inside index = " << indexSet.subIndex( *eIt, indexInInside, 1 )
                  << ", outside index = " << indexSet.subIndex( *outside, indexInOutside, 1 ) << std::endl;
        assert( false );
      }

      const typename GridType::LocalIdSet &localIdSet = grid.localIdSet();
      if( localIdSet.subId( *eIt, indexInInside, 1 ) != localIdSet.subId( *outside, indexInOutside, 1 ) )
      {
        std::cerr << "Error: Local id of conforming intersection differs when "
                  << "obtained from inside and outside." << std::endl;
        std::cerr << "       inside id = " << localIdSet.subId( *eIt, indexInInside, 1 )
                  << ", outside id = " << localIdSet.subId( *outside, indexInOutside, 1 ) << std::endl;
        assert( false );
      }

      const typename GridType::GlobalIdSet &globalIdSet = grid.globalIdSet();
      if( globalIdSet.subId( *eIt, indexInInside, 1 ) != globalIdSet.subId( *outside, indexInOutside, 1 ) )
      {
        std::cerr << "Error: Global id of conforming intersection differs when "
                  << "obtained from inside and outside." << std::endl;
        std::cerr << "       inside id = " << globalIdSet.subId( *eIt, indexInInside, 1 )
                  << ", outside id = " << globalIdSet.subId( *outside, indexInOutside, 1 ) << std::endl;
        assert( false );
      }
    }

    // //////////////////////////////////////////////////////////
    //   Check the geometry returned by geometry()
    // //////////////////////////////////////////////////////////
    typedef typename Intersection::Geometry Geometry;
    typedef typename Geometry :: ctype ctype ;
    const Geometry &geometry = iIt->geometry();

    checkGeometry(geometry);

    if( geometry.type() != iIt->type() )
    {
      std::cerr << "Error: Reference geometry type returned by intersection "
                << "differs from the one returned by the global geometry." << std::endl;
      assert( false );
    }

    // //////////////////////////////////////////////////////////
    //   Check the geometry returned by geometryInInside()
    // //////////////////////////////////////////////////////////

    typedef typename Intersection::LocalGeometry LocalGeometry;
    typedef typename Intersection::Geometry IntersectionGeometry;
    const LocalGeometry &geometryInInside = iIt->geometryInInside();
    checkLocalGeometry( geometryInInside, eIt->type(), "geometryInInside" );

    //  Check the consistency of geometryInInside() and geometry

    if (geometryInInside.corners() != geometry.corners())
      DUNE_THROW(GridError, "Geometry of intersection is inconsistent from left hand side and global view!");

    // Use a quadrature rule as a set of test points
    typedef Dune::GenericGeometry::GaussPoints< ctype > Points;
    typedef Dune::GenericGeometry::GenericQuadratureFactory< interDim, ctype, Points > QuadratureFactory;
    const typename QuadratureFactory::Object &quad = *QuadratureFactory::create( iIt->type(), 3 );
    for (size_t i=0; i<quad.size(); i++)
    {
      const typename LocalGeometry::LocalCoordinate &pt = quad[ i ].position();
      const typename IntersectionGeometry::Jacobian &jit
        = geometry.jacobianInverseTransposed( pt );

      // independently calculate the integration outer normal for the inside element
      typename LocalGeometry::GlobalCoordinate xInside = geometryInInside.global( pt );
      const typename LocalGeometry::GlobalCoordinate &refNormal = refElement.volumeOuterNormal( indexInInside );
      typename IntersectionGeometry::GlobalCoordinate refIntNormal;
      geoInside.jacobianInverseTransposed( xInside ).mv( refNormal, refIntNormal );
      // note: refElement.template mapping< 1 >( indexInInside ) is affine,
      //       hence we may use any point to obtain the integrationElement
      refIntNormal *= geoInside.integrationElement( xInside ) * geometryInInside.integrationElement( pt )
                      / refElement.template mapping< 1 >( indexInInside ).integrationElement( pt );

      // Check outer normal
      // const typename IntersectionGeometry::GlobalCoordinate normal = iIt->outerNormal( pt );
      const typename Intersection::GlobalCoordinate normal = iIt->outerNormal( pt );
      checkParallel( normal, refIntNormal, "outerNormal" );

      // Check normal vector is orthogonal to all vectors connecting
      // the vertices
      for (int c=1; c<geometry.corners(); c++)
      {
        typename IntersectionGeometry::GlobalCoordinate x = geometry.corner( c-1 );
        x -= geometry.corner( c );
        // this only makes sense for non-curved faces
        if( geometry.affine() && x*normal >= 1e3*std::numeric_limits< ctype >::epsilon() )
        {
          std::cerr << "outerNormal not orthogonal to line between corner "
                    << (c-1) << " and corner " << c << "." << std::endl;
          std::cerr << "Note: This is ok for curved faces, though." << std::endl;
          assert( false );
        }
      }

      // Check JacobianInverseTransposed
      // (J^-1 * n) = 0. Here J is the jacobian of the intersection
      // geometry (geometry) and n is a normal.
      checkJIn( normal, jit, "outerNormal" );

      // Check integration outer normal
      typename IntersectionGeometry::GlobalCoordinate intNormal =
        iIt->integrationOuterNormal( pt );
      sumNormal.axpy( quad[ i ].weight(), intNormal );

      const ctype det = geometry.integrationElement( pt );
      if( std :: abs( det - intNormal.two_norm() ) > 1e-8 )
      {
        std :: cerr << "Error: integrationOuterNormal yields wrong length."
                    << std :: endl;
        std :: cerr << "       |integrationOuterNormal| = " << intNormal.two_norm()
                    << ", integrationElement = " << det << std :: endl;
        assert( false );
      }

      checkJIn( intNormal, jit, "integrationOuterNormal" );
      checkParallel ( intNormal, refIntNormal, "integrationOuterNormal");
      if( (intNormal - refIntNormal).two_norm() > 1e-8 )
      {
        std::cerr << "Error: Wrong integration outer normal (" << intNormal
                  << ", should be " << refIntNormal << ")." << std::endl;
        std::cerr << "       Intersection: " << geometry.corner( 0 );
        for( int c = 1; c < geometry.corners(); ++c )
          std::cerr << ", " << geometry.corner( c );
        std::cerr << "." << std::endl;
        assert( false );
      }

      // Check unit outer normal
      // const typename IntersectionGeometry::GlobalCoordinate unitNormal = iIt->unitOuterNormal( pt );
      const typename Intersection::GlobalCoordinate unitNormal = iIt->unitOuterNormal( pt );
      if( std :: abs( ctype( 1 ) - unitNormal.two_norm() ) > 1e-8 )
      {
        std :: cerr << "Error: unitOuterNormal yields wrong length." << std :: endl;
        std :: cerr << "       |unitOuterNormal| = " << unitNormal.two_norm()
                    << std :: endl;
        assert( false );
      }

      const bool isCartesian = Dune :: Capabilities :: isCartesian< GridType > :: v ;
      // check normal for Cartesian grids
      if( isCartesian )
      {
        if( ! geometry.affine() )
          DUNE_THROW( Dune::GridError, "Intersection geometry is not affine, although isCartesian is true");

        // check that normal of Cartesian grid is given in
        // a specific way
        typename Intersection::GlobalCoordinate normal ( 0 );
        normal[ indexInInside / 2 ] = 2 * ( indexInInside % 2 ) - 1;

        if( (normal - unitNormal).infinity_norm() > 1e-8 )
          DUNE_THROW( Dune::GridError, "Normal is not in Cartesian format, although isCartesian is true");
      }

      checkJIn( unitNormal, jit, "unitOuterNormal" );
      checkParallel ( unitNormal, refIntNormal, "unitOuterNormal");

      // check geometryInInside
      typename IntersectionGeometry::GlobalCoordinate globalPos = geometry.global(quad[i].position());
      typename IntersectionGeometry::GlobalCoordinate localPos  = eIt->geometry().global(geometryInInside.global(quad[i].position()));

      if( (globalPos - localPos).infinity_norm() > 1e-6 )
      {
        std::cerr << "Error: inside()->geometry().global( intersection.geometryInInside().global( x ) ) != intersection.geometry().global( x )" << std::endl;

        std::cerr << "       x = " << quad[ i ].position() << std::endl;

        std::cerr << "       interesction.geometry() = " << geometry.corner( 0 );
        for( int i = 1; i < geometry.corners(); ++i )
          std::cerr << " | " << geometry.corner( i );
        std::cerr << std::endl;

        std::cerr << "       intersection.geometryInInside() = " << geometryInInside.corner( 0 );
        for( int i = 1; i < geometryInInside.corners(); ++i )
          std::cerr << " | " << geometryInInside.corner( i );
        std::cerr << std::endl;

        std::cerr << "       inside()->geometry() = " << eIt->geometry().corner( 0 );
        for( int i = 1; i < eIt->geometry().corners(); ++i )
          std::cerr << " | " << eIt->geometry().corner( i );
        std::cerr << std::endl;

        assert( false );
      }
    }
    QuadratureFactory::release( &quad );

    if( (iIt->centerUnitOuterNormal() - iIt->unitOuterNormal( refFace.position( 0, 0 ) )).two_norm() > 1e-8 )
    {
      std::cerr << "Error: centerUnitOuterNormal() does not match unitOuterNormal( "
                << refFace.position( 0, 0 ) << " )." << std::endl;
      std::cerr << "       centerUnitOuterNormal = " << iIt->centerUnitOuterNormal()
                << ", unitOuterNormal( " << refFace.position( 0, 0 ) << " ) = "
                << iIt->unitOuterNormal( refFace.position( 0, 0 ) ) << std::endl;
    }

    // ////////////////////////////////////////////////////////////////
    //   Check the geometry returned by geometryInOutside()
    // ////////////////////////////////////////////////////////////////

    if( iIt->neighbor() )
    {
      const typename Intersection::LocalGeometry &geometryInOutside = iIt->geometryInOutside();
      checkLocalGeometry( geometryInOutside, iIt->outside()->type(), "geometryInOutside" );

      if (geometryInInside.corners() != geometryInOutside.corners())
        DUNE_THROW(GridError, "Geometry of intersection is inconsistent from left and right hand side!");

      if( !iIt->boundary() )
      {
        // (Ab)use a quadrature rule as a set of test points
        const int interDim = Intersection::LocalGeometry::mydimension;
        typedef Dune::GenericGeometry::GaussPoints< ctype > Points;
        typedef Dune::GenericGeometry::GenericQuadratureFactory< interDim, ctype, Points > QuadratureFactory;
        const typename QuadratureFactory::Object &quad = *QuadratureFactory::create( iIt->type(), 3 );
        for (size_t i=0; i<quad.size(); i++)
        {
          typename IntersectionGeometry::GlobalCoordinate globalPos = geometry.global(quad[i].position());
          typename IntersectionGeometry::GlobalCoordinate localPos  = iIt->outside()->geometry().global(geometryInOutside.global(quad[i].position()));

          if ( (globalPos - localPos).infinity_norm() > 1e-6)
            DUNE_THROW(GridError, "global( geometryInOutside(global() ) is not the same as geometry.global() at " << quad[i].position() << "!");

        }
        QuadratureFactory::release( &quad );
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
  if( (sumNormal.two_norm() > 1e-8) && (eIt->partitionType() != Dune::GhostEntity) )
  {
    // std :: cout << "Warning: Integral over outer normals is nonzero: "
    //             << sumNormal << std :: endl;
    ++errorState.sumNormalsNonZero;
  }

  // check assignment operator for IntersectionIterator
  checkIntersectionAssignment( view, eIt );
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
                << " nonzero curvature." << std :: endl;;
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