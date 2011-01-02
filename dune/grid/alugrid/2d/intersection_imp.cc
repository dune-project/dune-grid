// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRID_INTERSECTION_IMP_CC
#define DUNE_ALU2DGRID_INTERSECTION_IMP_CC

#include <stack>
#include <utility>

#include <dune/grid/common/genericreferenceelements.hh>
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

  //constructor for end iterator
  template<class GridImp>
  inline ALU2dGridIntersectionBase<GridImp> ::
  ALU2dGridIntersectionBase(const FactoryType& factory, int wLevel) :
    intersectionGlobal_(GeometryImp()),
    intersectionSelfLocal_(LocalGeometryImp()),
    intersectionNeighborLocal_(LocalGeometryImp()),
    factory_( factory ),
    localGeomStorage_( LocalGeometryStorageType :: instance() ),
    walkLevel_(wLevel)
  {
    this->done( 0 );
  }


  //! The copy constructor
  template<class GridImp>
  inline ALU2dGridIntersectionBase<GridImp> ::
  ALU2dGridIntersectionBase(const ALU2dGridIntersectionBase<GridImp> & org) :
    current(org.current),
    intersectionGlobal_(GeometryImp()),
    intersectionSelfLocal_(LocalGeometryImp()),
    intersectionNeighborLocal_(LocalGeometryImp()),
    factory_( org.factory_ ),
    localGeomStorage_( LocalGeometryStorageType :: instance() ),
    walkLevel_(org.walkLevel_)
  {}

  template<class GridImp>
  inline void
  ALU2dGridIntersectionBase<GridImp> ::
  assign(const ALU2dGridIntersectionBase<GridImp> & org)
  {
    assert( &factory_ == &org.factory_ );
    walkLevel_ = org.walkLevel_;
    current = org.current;

    // unset geometry information
    unsetUp2Date();
  }

  //! check whether entities are the same or whether iterator is done
  template<class GridImp>
  inline bool ALU2dGridIntersectionBase< GridImp >
  ::equals ( const ALU2dGridIntersectionBase< GridImp > &other ) const
  {
    return ((current.inside() == other.current.inside()) && (current.index_ == other.current.index_));
  }


  //! return level of inside() entitiy
  template<class GridImp>
  inline int ALU2dGridIntersectionBase<GridImp> :: level () const
  {
    assert( current.inside() );
    return current.inside()->level();
  }


  //! reset IntersectionIterator to first neighbour
  template<class GridImp>
  inline void ALU2dGridIntersectionBase<GridImp> ::
  first ( const EntityImp &en, int wLevel )
  {
    setFirstItem(en.getItem(),wLevel);
#if ALU2DGRID_PARALLEL
    checkValid();
#endif
  }


  //! return true if intersection is with boundary
  template<class GridImp>
  inline void ALU2dGridIntersectionBase<GridImp> :: checkValid ()
  {
    if( current.outside() )
    {
      const int index = current.outside()->getIndex();
      if( !this->grid().rankManager().isValid( index, All_Partition ) )
        current.setOutside( 0, -222 );
    }
  }

  //! return true if intersection is with boundary
  template<class GridImp>
  inline bool ALU2dGridIntersectionBase<GridImp> :: boundary() const
  {
    return current.isBoundary();
  }

  template<class GridImp>
  inline int ALU2dGridIntersectionBase<GridImp> :: boundaryId() const
  {
    assert( current.inside() );
    // ALUGrid stores negative values, so make 'em positive
    return (current.isBoundary() ? std::abs( current.boundary()->type() ) : 0);
  }

  template<class GridImp>
  inline size_t ALU2dGridIntersectionBase<GridImp> :: boundarySegmentIndex() const
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
  template<class GridImp>
  inline bool ALU2dGridIntersectionBase<GridImp> :: neighbor () const
  {
    return bool( current.outside() );
  }

  //! return EntityPointer to the Entity on the inside of this intersection.
  template< class GridImp >
  inline typename ALU2dGridIntersectionBase< GridImp >::EntityPointer
  ALU2dGridIntersectionBase< GridImp >::inside() const
  {
    assert( (current.inside() != 0) && (current.index_ < current.nFaces()) );
    return EntityPointerImp( factory_, *current.inside(), -1, walkLevel_ );
  }

  template< class GridImp >
  inline void ALU2dGridIntersectionBase<GridImp>::done ( const HElementType *inside )
  {
    current.setInside( const_cast< HElementType * >( inside ) );
    current.setOutside( 0, -1 );
    current.index_= current.nFaces();
  }

  //! return EntityPointer to the Entity on the outside of this intersection.
  template< class GridImp >
  inline typename ALU2dGridIntersectionBase< GridImp >::EntityPointer
  ALU2dGridIntersectionBase< GridImp >::outside() const
  {
    assert( current.inside() && current.outside() );
    return EntityPointerImp( factory_, *current.outside(), -1, walkLevel_ );
  }

  //! local number of codim 1 entity in self where intersection is contained in
  template<class GridImp>
  inline int ALU2dGridIntersectionBase<GridImp>::indexInInside () const
  {
    const int i = current.index_;
    if( (eltype == ALU2DSPACE triangle) || ((eltype == ALU2DSPACE mixed) && (current.nFaces() == 3)) )
      return 2 - i;
    else
      return ((i^2)>>1) | ((i&1)<<1);
  }

  //! local number of codim 1 entity in neighbor where intersection is contained in
  template<class GridImp>
  inline int ALU2dGridIntersectionBase<GridImp>::indexInOutside () const
  {
    const int i = current.opposite();
    if( (eltype == ALU2DSPACE triangle) || ((eltype == ALU2DSPACE mixed) && (current.nFaces() == 3)) )
      return 2 - i;
    else
      return ((i^2)>>1) | ((i&1)<<1);
  }

  template< class GridImp >
  inline int ALU2dGridIntersectionBase< GridImp >::twistInInside () const
  {
    return 0;
  }

  template< class GridImp >
  inline int ALU2dGridIntersectionBase< GridImp >::twistInOutside () const
  {
    // twist is either 0 or 1 depending on the edge numbers
    if( (eltype == ALU2DSPACE triangle) || ((eltype == ALU2DSPACE mixed) && (current.nFaces() == 3)) )
      return (1 + current.index_ + current.opposite()) % 2;
    else
      return (current.index_ + current.opposite()) % 2;
  }


#ifndef ALUGRID_SURFACE_2D
  template< class GridImp >
  inline typename ALU2dGridIntersectionBase<GridImp>::NormalType
  ALU2dGridIntersectionBase< GridImp >::outerNormal ( const LocalCoordinate &local ) const
  {
    assert( (current.inside() != 0) && (current.index_ < current.nFaces()) );

    const double (&u)[ dimworld ] = current.inside()->vertex( current.index_+1 )->coord();
    const double (&v)[ dimworld ] = current.inside()->vertex( current.index_+2 )->coord();

    NormalType outerNormal;
    outerNormal[ 0 ] = v[ 1 ] - u[ 1 ];
    outerNormal[ 1 ] = u[ 0 ] - v[ 0 ];

    if( current.useOutside_ )
      outerNormal *= 0.5;
    return outerNormal;
  }
#endif // #ifndef ALUGRID_SURFACE_2D

#ifdef ALUGRID_SURFACE_2D
  template< class GridImp >
  inline typename ALU2dGridIntersectionBase<GridImp>::NormalType
  ALU2dGridIntersectionBase< GridImp >::outerNormal ( const LocalCoordinate &local ) const
  {
    assert( (current.inside() != 0) && (current.index_ < current.nFaces()) );

    typedef double (&normal_t)[dimworld];

    NormalType outerNormal;
    if ( dimworld == 2 || current.inside()->numvertices() == 3 ) // current.inside()->affine()
      current.inside()->outernormal( current.index_, (normal_t)(&outerNormal)[0] );
    else
    {
      const GenericReferenceElement< alu2d_ctype, dim > &refElement =
        GenericReferenceElements< alu2d_ctype, dim >::cube();
      typename LocalGeometry::GlobalCoordinate xInside = geometryInInside().global( local );
      typename LocalGeometry::GlobalCoordinate refNormal = refElement.volumeOuterNormal( indexInInside() );
      inside()->geometry().jacobianInverseTransposed( xInside ).mv( refNormal, outerNormal );
      outerNormal *= inside()->geometry().integrationElement( xInside );
    }
    if( current.useOutside_ )
      outerNormal *= 0.5;
    return outerNormal;
  }
#endif // #ifdef ALUGRID_SURFACE_2D


  template< class GridImp >
  inline typename ALU2dGridIntersectionBase<GridImp>::NormalType
  ALU2dGridIntersectionBase< GridImp >::integrationOuterNormal ( const LocalCoordinate &local ) const
  {
    return outerNormal( local );
  }

  template< class GridImp >
  inline typename ALU2dGridIntersectionBase<GridImp>::NormalType
  ALU2dGridIntersectionBase< GridImp >::unitOuterNormal ( const LocalCoordinate &local ) const
  {
    NormalType unitNormal( outerNormal( local ) );
    unitNormal *= (1.0 / unitNormal.two_norm());
    return unitNormal;
  }


  template<class GridImp>
  inline const typename ALU2dGridIntersectionBase<GridImp>::LocalGeometry&
  ALU2dGridIntersectionBase< GridImp >::geometryInInside () const
  {
    assert( (current.inside() != 0) && (current.index_ < current.nFaces()) );

    // only in non-conform situation we use default method
    if( current.useOutside_ )
    {
      if( !GridImp::getRealImplementation( intersectionSelfLocal_ ).up2Date() )
        GridImp::getRealImplementation( intersectionSelfLocal_ ).buildLocalGeom( inside()->geometry(), geometry() );
      assert( GridImp::getRealImplementation( intersectionSelfLocal_ ).up2Date() );
      return intersectionSelfLocal_;
    }
    else
    {
      // parameters are face and twist
      const int localTwist = (current.nFaces() == 3) ? (current.index_ % 2) : (current.index_>>1)^(current.index_&1);
      const int twist = (twistInInside() + localTwist) % 2;
      return localGeomStorage_.localGeom( current.index_, twist, current.nFaces() );
    }
  }

  template<class GridImp>
  inline const typename ALU2dGridIntersectionBase<GridImp>::LocalGeometry&
  ALU2dGridIntersectionBase< GridImp >::geometryInOutside () const
  {
    assert( current.inside() && current.outside() );

    // only in non-conform situation we use default method
    if( (current.nFaces() != 3) || !conforming() )
    {
      if( !GridImp::getRealImplementation( intersectionNeighborLocal_ ).up2Date() )
        GridImp::getRealImplementation( intersectionNeighborLocal_ ).buildLocalGeom( outside()->geometry(), geometry() );
      assert( GridImp::getRealImplementation( intersectionNeighborLocal_ ).up2Date());
      return intersectionNeighborLocal_;
    }
    else
    {
      // parameters are face and twist
      const int localTwist = (current.nFaces() == 3) ? (current.opposite() % 2) : (current.opposite() >> 1)^(current.opposite() & 1);
      const int twist = (twistInOutside() + localTwist) % 2;
      return localGeomStorage_.localGeom( current.opposite(), twist, current.nFaces()  );
    }
  }

  template<class GridImp>
  inline const typename ALU2dGridIntersectionBase<GridImp>::Geometry&
  ALU2dGridIntersectionBase<GridImp>::geometry () const
  {
    assert( current.inside() );

    if( !GridImp::getRealImplementation( intersectionGlobal_ ).up2Date() )
    {
      if( current.useOutside_ )
        GridImp::getRealImplementation( intersectionGlobal_ ).buildGeom( *current.outside(), current.opposite() );
      else
        GridImp::getRealImplementation( intersectionGlobal_ ).buildGeom( *current.inside(), current.index_ );
    }

    assert( GridImp::getRealImplementation( intersectionGlobal_ ).up2Date() );
    return intersectionGlobal_;
  }


  template< class GridImp >
  inline GeometryType ALU2dGridIntersectionBase< GridImp >::type () const
  {
    return GeometryType( (eltype == ALU2DSPACE triangle ?
                          GenericGeometry :: SimplexTopology< 1 > :: type :: id :
                          GenericGeometry :: CubeTopology   < 1 > :: type :: id), 1);
  }

  template< class GridImp >
  inline void ALU2dGridIntersectionBase< GridImp >:: unsetUp2Date ()
  {
    GridImp::getRealImplementation( intersectionGlobal_ ).unsetUp2Date();
    GridImp::getRealImplementation( intersectionSelfLocal_ ).unsetUp2Date();
    GridImp::getRealImplementation( intersectionNeighborLocal_ ).unsetUp2Date();
  }


  //********************************************************************
  //
  //  --ALU2dGridLevelIntersectionIterator
  //  --LevelIntersectionIterator
  //
  //********************************************************************

  //! Constructor
  template<class GridImp>
  inline ALU2dGridLevelIntersectionIterator<GridImp> ::
  ALU2dGridLevelIntersectionIterator(const FactoryType& factory, const HElementType* el, int wLevel, bool end)
    : ALU2dGridIntersectionBase<GridImp>::ALU2dGridIntersectionBase( factory, wLevel )
  {
    if (!end)
    {
      assert(this->walkLevel_ >= 0);
      setFirstItem(*el,wLevel);
    }
    else
      this->done( el );
  }

  template<class GridImp>
  inline ALU2dGridLevelIntersectionIterator<GridImp> ::
  ALU2dGridLevelIntersectionIterator(const FactoryType& factory, int wLevel) :
    ALU2dGridIntersectionBase<GridImp>::ALU2dGridIntersectionBase(factory, wLevel)
  {}

  //! The copy constructor
  template<class GridImp>
  inline ALU2dGridLevelIntersectionIterator<GridImp> ::
  ALU2dGridLevelIntersectionIterator(const ALU2dGridLevelIntersectionIterator<GridImp> & org) :
    ALU2dGridIntersectionBase<GridImp>:: ALU2dGridIntersectionBase(org)
  {
    nbStack_ = org.nbStack_;
  }

  //! The copy constructor
  template<class GridImp>
  inline void
  ALU2dGridLevelIntersectionIterator<GridImp> ::
  assign(const ALU2dGridLevelIntersectionIterator<GridImp> & org) {
    ALU2dGridIntersectionBase<GridImp>:: ALU2dGridIntersectionBase::assign(org);
    nbStack_ = org.nbStack_;
  }


  //! increment iterator
  template<class GridImp>
  inline void ALU2dGridLevelIntersectionIterator<GridImp> :: increment ()
  {
    doIncrement();
  #if ALU2DGRID_PARALLEL
    this->checkValid();
  #endif
  }


  //! reset IntersectionIterator to first neighbour
  template<class GridImp>
  template<class EntityType>
  inline void ALU2dGridLevelIntersectionIterator<GridImp> ::
  first(const EntityType & en, int wLevel)
  {
    setFirstItem(en.getItem(),wLevel);
  #if ALU2DGRID_PARALLEL
    this->checkValid();
  #endif
  }

  //********************************************************************
  //
  //  --ALU2dGridLeafIntersectionIterator
  //  --LeafIntersectionIterator
  //
  //********************************************************************


  // Constructor
  template<class GridImp>
  inline ALU2dGridLeafIntersectionIterator<GridImp> ::
  ALU2dGridLeafIntersectionIterator(const FactoryType& factory, const HElementType* el, int wLevel, bool end)
    : ALU2dGridIntersectionBase<GridImp>::ALU2dGridIntersectionBase( factory, wLevel )
  {
    if (!end)
    {
      assert(this->walkLevel_ >= 0);
      setFirstItem(*el,wLevel);
    }
    else
      done( el );
  }

  template<class GridImp>
  inline ALU2dGridLeafIntersectionIterator<GridImp> ::
  ALU2dGridLeafIntersectionIterator(const FactoryType& factory, int wLevel) :
    ALU2dGridIntersectionBase<GridImp>::ALU2dGridIntersectionBase(factory, wLevel)
  {}


  //! The copy constructor
  template<class GridImp>
  inline ALU2dGridLeafIntersectionIterator<GridImp> ::
  ALU2dGridLeafIntersectionIterator(const ALU2dGridLeafIntersectionIterator<GridImp> & org)
    :  ALU2dGridIntersectionBase<GridImp>::ALU2dGridIntersectionBase(org)
      ,  nbStack_(org.nbStack_)
  {}

  //! The copy constructor
  template<class GridImp>
  inline void
  ALU2dGridLeafIntersectionIterator<GridImp> ::
  assign(const ALU2dGridLeafIntersectionIterator<GridImp> & org){
    ALU2dGridIntersectionBase<GridImp>::ALU2dGridIntersectionBase::assign(org);
    nbStack_ = org.nbStack_;
  }


  //! increment iterator
  template<class GridImp>
  inline void ALU2dGridLeafIntersectionIterator<GridImp> :: increment ()
  {
    doIncrement();
  #if ALU2DGRID_PARALLEL
    this->checkValid();
  #endif
  }

  //! reset IntersectionIterator to first neighbour
  template<class GridImp>
  template<class EntityType>
  inline void ALU2dGridLeafIntersectionIterator<GridImp> ::
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
