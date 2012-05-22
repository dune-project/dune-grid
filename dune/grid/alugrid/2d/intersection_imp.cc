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

  //********************************************************************
  //
  //  --ALU2dGridIntersectionBase
  //
  //
  //********************************************************************

  template< class Grid >
  inline ALU2dGridIntersectionBase< Grid >
  ::ALU2dGridIntersectionBase ( const Factory &factory, int wLevel )
    : factory_( factory ),
      localGeomStorage_( LocalGeometryStorageType::instance() ),
      walkLevel_( wLevel )
  {
    this->done( 0 );
  }


  template< class Grid >
  inline ALU2dGridIntersectionBase< Grid >
  ::ALU2dGridIntersectionBase ( const This &other )
    : current( other.current ),
      factory_( other.factory_ ),
      localGeomStorage_( LocalGeometryStorageType::instance() ),
      walkLevel_( other.walkLevel_ )
  {}


  template< class Grid >
  inline const typename ALU2dGridIntersectionBase< Grid >::This &
  ALU2dGridIntersectionBase< Grid >::operator= ( const This &other )
  {
    current = other.current;
    assert( &factory_ == &other.factory_ );
    walkLevel_ = other.walkLevel_;

    invalidate();
    return *this;
  }


#if 0
  template<class Grid>
  inline void
  ALU2dGridIntersectionBase<Grid> ::
  assign(const ALU2dGridIntersectionBase<Grid> & org)
  {
    assert( &factory_ == &org.factory_ );
    walkLevel_ = org.walkLevel_;
    current = org.current;

    // unset geometry information
    invalidate();
  }
#endif

  //! check whether entities are the same or whether iterator is done
  template<class Grid>
  inline bool ALU2dGridIntersectionBase< Grid >
  ::equals ( const ALU2dGridIntersectionBase< Grid > &other ) const
  {
    return ((current.inside() == other.current.inside()) && (current.index_ == other.current.index_));
  }


  //! return level of inside() entitiy
  template<class Grid>
  inline int ALU2dGridIntersectionBase<Grid> :: level () const
  {
    assert( current.inside() );
    return current.inside()->level();
  }


  //! return true if intersection is with boundary
  template<class Grid>
  inline void ALU2dGridIntersectionBase<Grid> :: checkValid ()
  {
    if( current.outside() )
    {
      const int index = current.outside()->getIndex();
      if( !this->grid().rankManager().isValid( index, All_Partition ) )
        current.setOutside( 0, -222 );
    }
  }

  //! return true if intersection is with boundary
  template<class Grid>
  inline bool ALU2dGridIntersectionBase<Grid> :: boundary() const
  {
    return current.isBoundary();
  }

  template<class Grid>
  inline int ALU2dGridIntersectionBase<Grid> :: boundaryId() const
  {
    assert( current.inside() );
    // ALUGrid stores negative values, so make 'em positive
    return (current.isBoundary() ? std::abs( current.boundary()->type() ) : 0);
  }

  template<class Grid>
  inline size_t ALU2dGridIntersectionBase<Grid> :: boundarySegmentIndex() const
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
  template<class Grid>
  inline bool ALU2dGridIntersectionBase<Grid> :: neighbor () const
  {
    return bool( current.outside() );
  }

  //! return EntityPointer to the Entity on the inside of this intersection.
  template< class Grid >
  inline typename ALU2dGridIntersectionBase< Grid >::EntityPointer
  ALU2dGridIntersectionBase< Grid >::inside() const
  {
    assert( (current.inside() != 0) && (current.index_ < current.nFaces()) );
    return EntityPointerImp( factory_, *current.inside(), -1, walkLevel_ );
  }


  template< class Grid >
  inline void ALU2dGridIntersectionBase< Grid >::done ( HElementType *inside )
  {
    current.setInside( inside );
    current.setOutside( 0, -1 );
    current.index_= current.nFaces();
  }


  template< class Grid >
  inline typename ALU2dGridIntersectionBase< Grid >::EntityPointer
  ALU2dGridIntersectionBase< Grid >::outside() const
  {
    assert( current.inside() && current.outside() );
    return EntityPointerImp( factory_, *current.outside(), -1, walkLevel_ );
  }

  template< class Grid >
  inline int ALU2dGridIntersectionBase< Grid >::indexInInside () const
  {
    const int i = current.index_;
    if( (eltype == ALU2DSPACE triangle) || ((eltype == ALU2DSPACE mixed) && (current.nFaces() == 3)) )
      return 2 - i;
    else
      return ((i^2)>>1) | ((i&1)<<1);
  }


  template< class Grid >
  inline int ALU2dGridIntersectionBase< Grid >::indexInOutside () const
  {
    const int i = current.opposite();
    if( (eltype == ALU2DSPACE triangle) || ((eltype == ALU2DSPACE mixed) && (current.nFaces() == 3)) )
      return 2 - i;
    else
      return ((i^2)>>1) | ((i&1)<<1);
  }


  template< class Grid >
  inline int ALU2dGridIntersectionBase< Grid >::twistInInside () const
  {
    return 0;
  }


  template< class Grid >
  inline int ALU2dGridIntersectionBase< Grid >::twistInOutside () const
  {
    const int i = current.index_;
    const int o = current.opposite();
    // The twist is either 0 or 1, depending on the edge numbers.
    // The edge is always twisted with respect to the ALU reference element.
    if( (eltype == ALU2DSPACE triangle) || ((eltype == ALU2DSPACE mixed) && (current.nFaces() == 3)) )
      //return (1 + i + o) % 2;
      return 1 ^ ((i^o) & 1);
    else
      return 1 ^ ((i^o) & 1) ^ ((i^o) >> 1);
  }


  template< class Grid >
  inline typename ALU2dGridIntersectionBase< Grid >::NormalType
  ALU2dGridIntersectionBase< Grid >::outerNormal ( const LocalCoordinate &local ) const
  {
    assert( (current.inside() != 0) && (current.index_ < current.nFaces()) );

    typedef double (&normal_t)[dimensionworld];

    NormalType outerNormal;
    if ( dimensionworld == 2 || current.inside()->numvertices() == 3 ) // current.inside()->affine()
      current.inside()->outernormal( current.index_, (normal_t)(&outerNormal)[0] );
    else
    {
      const GenericReferenceElement< alu2d_ctype, dimension > &refElement =
        GenericReferenceElements< alu2d_ctype, dimension >::cube();
      typename LocalGeometry::GlobalCoordinate xInside = geometryInInside().global( local );
      typename LocalGeometry::GlobalCoordinate refNormal = refElement.volumeOuterNormal( indexInInside() );
      inside()->geometry().jacobianInverseTransposed( xInside ).mv( refNormal, outerNormal );
      outerNormal *= inside()->geometry().integrationElement( xInside );
    }
    if( current.useOutside_ )
      outerNormal *= 0.5;
    return outerNormal;
  }


  template< class Grid >
  inline typename ALU2dGridIntersectionBase< Grid >::NormalType
  ALU2dGridIntersectionBase< Grid >::integrationOuterNormal ( const LocalCoordinate &local ) const
  {
    return outerNormal( local );
  }


  template< class Grid >
  inline typename ALU2dGridIntersectionBase< Grid >::NormalType
  ALU2dGridIntersectionBase< Grid >::unitOuterNormal ( const LocalCoordinate &local ) const
  {
    NormalType unitNormal( outerNormal( local ) );
    unitNormal *= (1.0 / unitNormal.two_norm());
    return unitNormal;
  }


  template< class Grid >
  inline typename ALU2dGridIntersectionBase< Grid >::NormalType
  ALU2dGridIntersectionBase< Grid >::centerUnitOuterNormal () const
  {
    const GenericReferenceElement< ctype, dimension-1 > &refElement
      = GenericReferenceElements< ctype, dimension-1 >::general( type() );
    return unitOuterNormal( refElement.position( 0, 0 ) );
  }


  template< class Grid >
  inline typename ALU2dGridIntersectionBase< Grid >::LocalGeometry
  ALU2dGridIntersectionBase< Grid >::geometryInInside () const
  {
    assert( (current.inside() != 0) && (current.index_ < current.nFaces()) );

    // only in non-conform situation we use default method
    if( current.useOutside_ )
    {
      if( !intersectionSelfLocal_.valid() )
        intersectionSelfLocal_.buildLocalGeom( inside()->geometry(), geometry() );
      assert( intersectionSelfLocal_.valid() );
      return LocalGeometry( intersectionSelfLocal_ );
    }
    else
    {
      // parameters are face and twist
      const int localTwist = (current.nFaces() == 3) ? (current.index_ % 2) : (current.index_>>1)^(current.index_&1);
      const int twist = (twistInInside() + localTwist) % 2;
      return LocalGeometry( localGeomStorage_.localGeom( current.index_, twist, current.nFaces() ) );
    }
  }

  template< class Grid >
  inline typename ALU2dGridIntersectionBase< Grid >::LocalGeometry
  ALU2dGridIntersectionBase< Grid >::geometryInOutside () const
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

  template< class Grid >
  inline typename ALU2dGridIntersectionBase< Grid >::Geometry
  ALU2dGridIntersectionBase< Grid >::geometry () const
  {
    assert( current.inside() );

    if( !intersectionGlobal_.valid() )
    {
      if( current.useOutside_ )
        intersectionGlobal_.buildGeom( *current.outside(), current.opposite() );
      else
        intersectionGlobal_.buildGeom( *current.inside(), current.index_ );
    }

    assert( intersectionGlobal_.valid() );
    return Geometry( intersectionGlobal_ );
  }


  template< class Grid >
  inline GeometryType ALU2dGridIntersectionBase< Grid >::type () const
  {
    return GeometryType( (eltype == ALU2DSPACE triangle ?
                          GenericGeometry :: SimplexTopology< 1 > :: type :: id :
                          GenericGeometry :: CubeTopology   < 1 > :: type :: id), 1);
  }

  template< class Grid >
  inline void ALU2dGridIntersectionBase< Grid >::invalidate ()
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
    : intersection_( IntersectionImpl( factory, wLevel ) )
  {
    if( !end )
    {
      assert( walkLevel() >= 0 );

      current().setInside( el );
      current().index_ = -1;
      current().setOutside( 0, -1 );

      assert( current().inside() );

      increment();
    }
    else
      intersectionImpl().done( el );
  }


#if 0
  template<class Grid>
  inline ALU2dGridLevelIntersectionIterator<Grid> ::
  ALU2dGridLevelIntersectionIterator(const FactoryType& factory, int wLevel) :
    ALU2dGridIntersectionBase<Grid>::ALU2dGridIntersectionBase(factory, wLevel)
  {}
#endif


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


#if 0
  template<class Grid>
  inline void
  ALU2dGridLevelIntersectionIterator<Grid> ::
  assign(const ALU2dGridLevelIntersectionIterator<Grid> & org) {
    ALU2dGridIntersectionBase<Grid>:: ALU2dGridIntersectionBase::assign(org);
    nbStack_ = org.nbStack_;
  }
#endif


  //! increment iterator
  template<class Grid>
  inline void ALU2dGridLevelIntersectionIterator<Grid> :: increment ()
  {
    doIncrement();
  #if ALU2DGRID_PARALLEL
    this->checkValid();
  #endif
  }


  //! reset IntersectionIterator to first neighbour
  template<class Grid>
  template<class EntityType>
  inline void ALU2dGridLevelIntersectionIterator<Grid> ::
  first(const EntityType & en, int wLevel)
  {
    setFirstItem(en.getItem(),wLevel);
  #if ALU2DGRID_PARALLEL
    this->checkValid();
  #endif
  }



  // Implementation of ALU2dGridLeafIntersectionIterator
  // ---------------------------------------------------

  template< class Grid >
  inline ALU2dGridLeafIntersectionIterator< Grid >
  ::ALU2dGridLeafIntersectionIterator ( const Factory &factory, HElementType *el, int wLevel, bool end )
    : intersection_( IntersectionImpl( factory, wLevel ) )
  {
    if( !end )
    {
      assert( walkLevel() >= 0 );

      current().setInside( el );
      current().index_ = -1;
      current().setOutside( 0, -1 );

      assert( current().inside() );

      increment();
    }
    else
      intersectionImpl().done( el );
  }


#if 0
  template<class Grid>
  inline ALU2dGridLeafIntersectionIterator<Grid> ::
  ALU2dGridLeafIntersectionIterator(const FactoryType& factory, int wLevel) :
    ALU2dGridIntersectionBase<Grid>::ALU2dGridIntersectionBase(factory, wLevel)
  {}
#endif


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

#if 0
  //! The copy constructor
  template<class Grid>
  inline void
  ALU2dGridLeafIntersectionIterator<Grid> ::
  assign(const ALU2dGridLeafIntersectionIterator<Grid> & org){
    ALU2dGridIntersectionBase<Grid>::ALU2dGridIntersectionBase::assign(org);
    nbStack_ = org.nbStack_;
  }
#endif


  //! increment iterator
  template<class Grid>
  inline void ALU2dGridLeafIntersectionIterator<Grid> :: increment ()
  {
    doIncrement();
  #if ALU2DGRID_PARALLEL
    this->checkValid();
  #endif
  }

  //! reset IntersectionIterator to first neighbour
  template<class Grid>
  template<class EntityType>
  inline void ALU2dGridLeafIntersectionIterator<Grid> ::
  first(const EntityType & en, int wLevel)
  {
    // if called on non-leaf, just return end iterator
    if( en.isLeaf() )
      setFirstItem( en.getItem(), wLevel );
    else
      done( &en.getItem() );
  }

} //end namespace Dune

#endif // #ifndef DUNE_ALU2DGRID_INTERSECTION_IMP_CC
