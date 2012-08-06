// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRID_INTERSECTION_IMP_CC
#define DUNE_ALU2DGRID_INTERSECTION_IMP_CC

#include <stack>
#include <utility>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/alugrid/2d/geometry.hh>
#include <dune/grid/alugrid/2d/entity.hh>
#include <dune/grid/alugrid/2d/grid.hh>

namespace Dune
{

  // Implementation of ALU2dGridIntersectionBase
  // -------------------------------------------

  template< class Grid, class Info >
  inline ALU2dGridIntersectionBase< Grid, Info >
  ::ALU2dGridIntersectionBase ( const Factory &factory, const IntersectionInfo &info )
    : current( info ),
      factory_( factory ),
      localGeomStorage_( LocalGeometryStorageType::instance() )
  {}


  template< class Grid, class Info >
  inline ALU2dGridIntersectionBase< Grid, Info >
  ::ALU2dGridIntersectionBase ( const This &other )
    : current( other.current ),
      factory_( other.factory_ ),
      localGeomStorage_( LocalGeometryStorageType::instance() )
  {}


  template< class Grid, class Info >
  inline const typename ALU2dGridIntersectionBase< Grid, Info >::This &
  ALU2dGridIntersectionBase< Grid, Info >::operator= ( const This &other )
  {
    current = other.current;
    assert( &factory_ == &other.factory_ );

    invalidate();
    return *this;
  }


  template< class Grid, class Info >
  inline void ALU2dGridIntersectionBase< Grid, Info > :: checkValid ()
  {
    if( current.outside() )
    {
      const int index = current.outside()->getIndex();
      if( !this->grid().rankManager().isValid( index, All_Partition ) )
        current.setOutside( 0, -222 );
    }
  }


  template< class Grid, class Info >
  inline int ALU2dGridIntersectionBase< Grid, Info > :: boundaryId() const
  {
    assert( current.inside() );
    // ALUGrid stores negative values, so make 'em positive
    return (current.isBoundary() ? std::abs( current.boundary()->type() ) : 0);
  }

  template< class Grid, class Info >
  inline size_t ALU2dGridIntersectionBase< Grid, Info > :: boundarySegmentIndex() const
  {
    // only call this method on boundary intersections
    assert( current.isBoundary() );
#ifdef ALUGRID_VERTEX_PROJECTION
    return current.boundary()->segmentIndex();
#else
    derr << "Method available in any version of ALUGrid > 1.14 \n";
    return 0;
#endif
  }

  //! return true if intersection is with neighbor on this level
  template< class Grid, class Info >
  inline bool ALU2dGridIntersectionBase< Grid, Info > :: neighbor () const
  {
    return bool( current.outside() );
  }

  //! return EntityPointer to the Entity on the inside of this intersection.
  template< class Grid, class Info >
  inline typename ALU2dGridIntersectionBase< Grid, Info >::EntityPointer
  ALU2dGridIntersectionBase< Grid, Info >::inside() const
  {
    assert( (current.inside() != 0) && (current.index() < current.nFaces()) );
    return EntityPointerImp( factory_, *current.inside(), -1, current.walkLevel() );
  }


  template< class Grid, class Info >
  inline typename ALU2dGridIntersectionBase< Grid, Info >::EntityPointer
  ALU2dGridIntersectionBase< Grid, Info >::outside() const
  {
    assert( current.inside() && current.outside() );
    return EntityPointerImp( factory_, *current.outside(), -1, current.walkLevel() );
  }

  template< class Grid, class Info >
  inline int ALU2dGridIntersectionBase< Grid, Info >::indexInInside () const
  {
    const int i = current.index();
    if( (eltype == ALU2DSPACE triangle) || ((eltype == ALU2DSPACE mixed) && (current.nFaces() == 3)) )
      return 2 - i;
    else
      return ((i^2)>>1) | ((i&1)<<1);
  }


  template< class Grid, class Info >
  inline int ALU2dGridIntersectionBase< Grid, Info >::indexInOutside () const
  {
    const int i = current.opposite();
    if( (eltype == ALU2DSPACE triangle) || ((eltype == ALU2DSPACE mixed) && (current.nFaces() == 3)) )
      return 2 - i;
    else
      return ((i^2)>>1) | ((i&1)<<1);
  }


  template< class Grid, class Info >
  inline int ALU2dGridIntersectionBase< Grid, Info >::twistInInside () const
  {
    return 0;
  }


  template< class Grid, class Info >
  inline int ALU2dGridIntersectionBase< Grid, Info >::twistInOutside () const
  {
    const int i = current.index();
    const int o = current.opposite();
    // The twist is either 0 or 1, depending on the edge numbers.
    // The edge is always twisted with respect to the ALU reference element.
    if( (eltype == ALU2DSPACE triangle) || ((eltype == ALU2DSPACE mixed) && (current.nFaces() == 3)) )
      //return (1 + i + o) % 2;
      return 1 ^ ((i^o) & 1);
    else
      return 1 ^ ((i^o) & 1) ^ ((i^o) >> 1);
  }


  template< class Grid, class Info >
  inline typename ALU2dGridIntersectionBase< Grid, Info >::NormalType
  ALU2dGridIntersectionBase< Grid, Info >::outerNormal ( const LocalCoordinate &local ) const
  {
    assert( (current.inside() != 0) && (current.index() < current.nFaces()) );

    typedef double (&normal_t)[dimensionworld];

    NormalType outerNormal;
    if ( dimensionworld == 2 || current.inside()->numvertices() == 3 ) // current.inside()->affine()
      current.inside()->outernormal( current.index(), (normal_t)(&outerNormal)[0] );
    else
    {
      const ReferenceElement< alu2d_ctype, dimension > &refElement =
        ReferenceElements< alu2d_ctype, dimension >::cube();
      typename LocalGeometry::GlobalCoordinate xInside = geometryInInside().global( local );
      typename LocalGeometry::GlobalCoordinate refNormal = refElement.volumeOuterNormal( indexInInside() );
      inside()->geometry().jacobianInverseTransposed( xInside ).mv( refNormal, outerNormal );
      outerNormal *= inside()->geometry().integrationElement( xInside );
    }
    if( current.useOutside() )
      outerNormal *= 0.5;
    return outerNormal;
  }


  template< class Grid, class Info >
  inline typename ALU2dGridIntersectionBase< Grid, Info >::NormalType
  ALU2dGridIntersectionBase< Grid, Info >::integrationOuterNormal ( const LocalCoordinate &local ) const
  {
    return outerNormal( local );
  }


  template< class Grid, class Info >
  inline typename ALU2dGridIntersectionBase< Grid, Info >::NormalType
  ALU2dGridIntersectionBase< Grid, Info >::unitOuterNormal ( const LocalCoordinate &local ) const
  {
    NormalType unitNormal( outerNormal( local ) );
    unitNormal *= (1.0 / unitNormal.two_norm());
    return unitNormal;
  }


  template< class Grid, class Info >
  inline typename ALU2dGridIntersectionBase< Grid, Info >::NormalType
  ALU2dGridIntersectionBase< Grid, Info >::centerUnitOuterNormal () const
  {
    const ReferenceElement< ctype, dimension-1 > &refElement
      = ReferenceElements< ctype, dimension-1 >::general( type() );
    return unitOuterNormal( refElement.position( 0, 0 ) );
  }


  template< class Grid, class Info >
  inline typename ALU2dGridIntersectionBase< Grid, Info >::LocalGeometry
  ALU2dGridIntersectionBase< Grid, Info >::geometryInInside () const
  {
    const int i = current.index();
    assert( (current.inside() != 0) && (i < current.nFaces()) );

    // only in non-conform situation we use default method
    if( current.useOutside() )
    {
      if( !intersectionSelfLocal_.valid() )
        intersectionSelfLocal_.buildLocalGeom( inside()->geometry(), geometry() );
      assert( intersectionSelfLocal_.valid() );
      return LocalGeometry( intersectionSelfLocal_ );
    }
    else
    {
      // parameters are face and twist
      const int localTwist = (current.nFaces() == 3) ? (i % 2) : (i >> 1)^(i & 1);
      const int twist = (twistInInside() + localTwist) % 2;
      return LocalGeometry( localGeomStorage_.localGeom( i, twist, current.nFaces() ) );
    }
  }

  template< class Grid, class Info >
  inline typename ALU2dGridIntersectionBase< Grid, Info >::LocalGeometry
  ALU2dGridIntersectionBase< Grid, Info >::geometryInOutside () const
  {
    assert( current.inside() && current.outside() );

    // only in non-conform situation we use default method
    if( !conforming() )
    {
      if( !intersectionNeighborLocal_.valid() )
        intersectionNeighborLocal_.buildLocalGeom( outside()->geometry(), geometry() );
      assert( intersectionNeighborLocal_.valid() );
      return LocalGeometry( intersectionNeighborLocal_ );
    }
    else
    {
      // parameters are face and twist
      const int localTwist = (current.nFaces() == 3) ? (current.opposite() % 2) : (current.opposite() >> 1)^(current.opposite() & 1);
      const int twist = (twistInOutside() + localTwist) % 2;
      return LocalGeometry( localGeomStorage_.localGeom( current.opposite(), twist, current.nFaces() ) );
    }
  }

  template< class Grid, class Info >
  inline typename ALU2dGridIntersectionBase< Grid, Info >::Geometry
  ALU2dGridIntersectionBase< Grid, Info >::geometry () const
  {
    assert( current.inside() );

    if( !intersectionGlobal_.valid() )
    {
      if( current.useOutside() )
        intersectionGlobal_.buildGeom( *current.outside(), current.opposite() );
      else
        intersectionGlobal_.buildGeom( *current.inside(), current.index() );
    }

    assert( intersectionGlobal_.valid() );
    return Geometry( intersectionGlobal_ );
  }


  template< class Grid, class Info >
  inline GeometryType ALU2dGridIntersectionBase< Grid, Info >::type () const
  {
    return GeometryType( (eltype == ALU2DSPACE triangle ?
                          GenericGeometry :: SimplexTopology< 1 > :: type :: id :
                          GenericGeometry :: CubeTopology   < 1 > :: type :: id), 1);
  }

  template< class Grid, class Info >
  inline void ALU2dGridIntersectionBase< Grid, Info >::invalidate ()
  {
    intersectionGlobal_.invalidate();
    intersectionSelfLocal_.invalidate();
    intersectionNeighborLocal_.invalidate();
  }



  // Implementation of ALU2dGridLevelIntersectionIterator
  // ----------------------------------------------------

  template< class Grid >
  inline ALU2dGridLevelIntersectionIterator< Grid >
  ::ALU2dGridLevelIntersectionIterator ( const Factory &factory, HElementType *el, int wLevel, bool end )
    : intersection_( IntersectionImpl( factory ) )
  {
    current().setInside( el );
    current().walkLevel_ = wLevel;
    if( !end )
    {
      assert( current().walkLevel() >= 0 );
      assert( current().inside() );

      current().index_ = -1;
      current().setOutside( 0, -1 );

      increment();
    }
    else
      current().done();
  }


  template< class Grid >
  inline ALU2dGridLevelIntersectionIterator< Grid >
  ::ALU2dGridLevelIntersectionIterator ( const This &other )
    : intersection_( other.intersectionImpl() ),
      nbStack_( other.nbStack_ )
  {}


  template< class Grid >
  inline const typename ALU2dGridLevelIntersectionIterator< Grid >::This &
  ALU2dGridLevelIntersectionIterator< Grid >::operator= ( const This &other )
  {
    intersectionImpl() = other.intersectionImpl();
    nbStack_ = other.nbStack_;
    return *this;
  }


  //! increment iterator
  template<class Grid>
  inline void ALU2dGridLevelIntersectionIterator<Grid> :: increment ()
  {
    doIncrement();
  #if ALU2DGRID_PARALLEL
    this->checkValid();
  #endif
  }



  // Implementation of ALU2dGridLeafIntersectionIterator
  // ---------------------------------------------------

  template< class Grid >
  inline ALU2dGridLeafIntersectionIterator< Grid >
  ::ALU2dGridLeafIntersectionIterator ( const Factory &factory, HElementType *el, int wLevel, bool end )
    : intersection_( IntersectionImpl( factory ) )
  {
    current().setInside( el );
    current().walkLevel_ = wLevel;
    if( !end )
    {
      assert( current().walkLevel() >= 0 );
      assert( current().inside() );

      current().index_ = -1;
      current().setOutside( 0, -1 );

      increment();
    }
    else
      current().done();
  }


  template< class Grid >
  inline ALU2dGridLeafIntersectionIterator< Grid >
  ::ALU2dGridLeafIntersectionIterator ( const This &other )
    : intersection_( other.intersectionImpl() ),
      nbStack_( other.nbStack_ )
  {}


  template< class Grid >
  inline const typename ALU2dGridLeafIntersectionIterator< Grid >::This &
  ALU2dGridLeafIntersectionIterator< Grid >::operator= ( const This &other )
  {
    intersectionImpl() = other.intersectionImpl();
    nbStack_ = other.nbStack_;
    return *this;
  }


  //! increment iterator
  template<class Grid>
  inline void ALU2dGridLeafIntersectionIterator<Grid> :: increment ()
  {
    doIncrement();
  #if ALU2DGRID_PARALLEL
    this->checkValid();
  #endif
  }

} // namespace Dune

#endif // #ifndef DUNE_ALU2DGRID_INTERSECTION_IMP_CC
