// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRID_INTERSECTION_IMP_CC
#define DUNE_ALU2DGRID_INTERSECTION_IMP_CC

#include <stack>
#include <utility>

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

  //! Constructor
  template<class GridImp>
  inline ALU2dGridIntersectionBase<GridImp> ::
  ALU2dGridIntersectionBase(const GridImp & grid, const HElementType* el, int wLevel, bool end) :
    intersectionGlobal_(GeometryImp()),
    intersectionSelfLocal_(LocalGeometryImp()),
    intersectionNeighborLocal_(LocalGeometryImp()),
    grid_(grid),
    localGeomStorage_( LocalGeometryStorageType :: instance() ),
    walkLevel_(wLevel)
  {
    if (!end)
    {
      assert(walkLevel_ >= 0);
      setFirstItem(*el,wLevel);
    }
    else
      done( el );
  }

  //constructor for end iterator
  template<class GridImp>
  inline ALU2dGridIntersectionBase<GridImp> ::
  ALU2dGridIntersectionBase(const GridImp & grid, int wLevel) :
    intersectionGlobal_(GeometryImp()),
    intersectionSelfLocal_(LocalGeometryImp()),
    intersectionNeighborLocal_(LocalGeometryImp()),
    grid_(grid),
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
    grid_(org.grid_),
    localGeomStorage_( LocalGeometryStorageType :: instance() ),
    walkLevel_(org.walkLevel_)
  {}

  template<class GridImp>
  inline void
  ALU2dGridIntersectionBase<GridImp> ::
  assign(const ALU2dGridIntersectionBase<GridImp> & org)
  {
    assert( &grid_ == &org.grid_);
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
      if( !this->grid_.rankManager().isValid( index, All_Partition ) )
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
    return EntityPointerImp( grid_, *current.inside(), -1, walkLevel_ );
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
    return EntityPointerImp( grid_, *current.outside(), -1, walkLevel_ );
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
    return (1 + current.index_ + current.opposite()) % 2;
  }

  template< class GridImp >
  inline typename ALU2dGridIntersectionBase<GridImp>::NormalType
  ALU2dGridIntersectionBase< GridImp >::outerNormal ( const LocalCoordinate &local ) const
  {
    assert( (current.inside() != 0) && (current.index_ < current.nFaces()) );

    typedef double (&normal_t)[dimworld];

    NormalType outerNormal;
    current.inside()->outernormal( current.index_, (normal_t)(&outerNormal)[0] );
    if( current.useOutside_ )
      outerNormal *= 0.5;
    return outerNormal;
  }

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
      if( ! this->grid_.getRealImplementation(intersectionSelfLocal_).up2Date() )
      {
        this->grid_.getRealImplementation(intersectionSelfLocal_).
        buildLocalGeom( inside()->geometry() , geometry() );
      }
      assert(this->grid_.getRealImplementation(intersectionSelfLocal_).up2Date());
      return intersectionSelfLocal_;
    }
    else
    {
      const int twist = (current.nFaces() == 3) ? (current.index_ % 2) : (current.index_>>1)^(current.index_&1);
      // parameters are face and twist
      return localGeomStorage_.localGeom(  current.index_, twist, current.nFaces() );
    }
  }

  template<class GridImp>
  inline const typename ALU2dGridIntersectionBase<GridImp>::LocalGeometry&
  ALU2dGridIntersectionBase< GridImp >::geometryInOutside () const
  {
    assert( current.inside() && current.outside() );

    // only in non-conform situation we use default method
    //if( ! this->conforming() )
    //if( this->current.useOutside_ )
    {
      if( ! this->grid_.getRealImplementation(intersectionNeighborLocal_).up2Date() )
      {
        // we don't know here wether we have non-conform or conform situation on the neighbor
        this->grid_.getRealImplementation(intersectionNeighborLocal_).
        buildLocalGeom( outside()->geometry(), geometry() );
      }
      assert(this->grid_.getRealImplementation(intersectionNeighborLocal_).up2Date());
      return intersectionNeighborLocal_;
    }
    /*
       else
       {
       // parameters are face and twist
       const int twist = (current.nFaces() == 3) ? (current.opposite_ % 2) : (current.opposite_>>1)^(current.opposite_&1);
       return localGeomStorage_.localGeom( current.opposite_, 1-twist, current.nFaces()  );
       }
     */
  }

  template<class GridImp>
  inline const typename ALU2dGridIntersectionBase<GridImp>::Geometry&
  ALU2dGridIntersectionBase<GridImp>::geometry () const
  {
    assert( current.inside() );

    if( ! this->grid_.getRealImplementation(intersectionGlobal_).up2Date() )
    {
      if( current.useOutside_ )
        this->grid_.getRealImplementation(intersectionGlobal_).
        buildGeom(*current.outside(), current.opposite());
      else
        this->grid_.getRealImplementation(intersectionGlobal_).
        buildGeom(*current.inside(), current.index_);
    }

    assert(this->grid_.getRealImplementation(intersectionGlobal_).up2Date());
    return intersectionGlobal_;
  }


  template< class GridImp >
  inline GeometryType ALU2dGridIntersectionBase< GridImp >::type () const
  {
    return GeometryType(
             (eltype == ALU2DSPACE triangle ? GeometryType::simplex : GeometryType::cube),
             1 );
  }

  template< class GridImp >
  inline void ALU2dGridIntersectionBase< GridImp >:: unsetUp2Date ()
  {
    grid_.getRealImplementation(intersectionGlobal_).unsetUp2Date();
    grid_.getRealImplementation(intersectionSelfLocal_).unsetUp2Date();
    grid_.getRealImplementation(intersectionNeighborLocal_).unsetUp2Date();
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
  ALU2dGridLevelIntersectionIterator(const GridImp & grid, const HElementType* el, int wLevel, bool end) :
    ALU2dGridIntersectionBase<GridImp>::ALU2dGridIntersectionBase(grid, el, wLevel, end)
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
  ALU2dGridLevelIntersectionIterator(const GridImp & grid, int wLevel) :
    ALU2dGridIntersectionBase<GridImp>::ALU2dGridIntersectionBase(grid, wLevel)
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


  template<class GridImp>
  inline int ALU2dGridLevelIntersectionIterator<GridImp> ::
  getOppositeInFather(const int nrInChild, const int nrOfChild) const
  {
    int ret = (nrInChild==0) ? (2-nrOfChild) :
              ((nrInChild-nrOfChild==2 || nrInChild-nrOfChild==0) ? -1 : 0);
    assert(ret >= -1 && ret < 3);
    return ret;
  }

  template<class GridImp>
  inline int ALU2dGridLevelIntersectionIterator<GridImp> ::
  getOppositeInChild(const int nrInFather, const int nrOfChild) const
  {
    int ret = (nrInFather==0) ? (nrOfChild+1) :
              ((nrInFather-nrOfChild==1) ? -1 : 0);
    assert( ret >= -1 && ret < 3 );
    return ret;
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

  //! increment iterator
  template<class GridImp>
  inline void ALU2dGridLevelIntersectionIterator<GridImp> :: doIncrement ()
  {
    assert( current.index_ < current.nFaces() );

    this->unsetUp2Date();

    if( nbStack_.empty() )
    {
      ++current.index_;
      if( current.index_ >= current.nFaces() )
      {
        assert( current.index_ == current.nFaces() );
        return;
      }

      addNeighboursToStack();
      // if more then one element in stack we have non-conform intersection
      current.useOutside_ = (nbStack_.size() > 1);

      if( nbStack_.empty() )
        return current.setOutside( 0, -1 );
    }

    setupIntersection();

    assert( !current.outside() || (current.outside()->level() == walkLevel_) );
  }


  template< class GridImp >
  inline void
  ALU2dGridLevelIntersectionIterator< GridImp >::addNeighboursToStack ()
  {
    assert( current.index_ < current.nFaces() );

    ThinelementType *neighbor = current.inside()->neighbour( current.index_ );
    assert( neighbor );

    IntersectionInfo info;
    if( neighbor->thinis( ThinelementType::bndel_like ) )
    {
      HBndElType *bndel = (HBndElType *)neighbor;
      if( bndel->type() != HBndElType::periodic )
        return;

      PeriodicBndElType *bndnb = ((PeriodicBndElType *)bndel)->periodic_nb;
      assert( bndnb && bndnb->neighbour( 0 ) && bndnb->neighbour( 0 )->thinis( ThinelementType::element_like ) );
      info.first = (HElementType *)bndnb->neighbour( 0 );
      info.second = bndnb->opposite( 0 );
    }
    else
    {
      assert( neighbor->thinis( ThinelementType::element_like ) );
      info.first = (HElementType *)neighbor;
      info.second = current.inside()->opposite( current.index_ );
    }
    assert( info.first );

    while( info.first->level() > walkLevel_ )
    {
      info.second = getOppositeInFather( info.second, info.first->childNr() );
      assert( (info.second >= 0) && (info.second < current.nFaces()) );
      info.first = info.first->father();
    }

    if( info.first->level() >= walkLevel_ )
    {
      nbStack_.push( info );
      return;
    }

    // why should we go up, here?
    while( info.first && (info.first->level() < walkLevel_ - 1) )
    {
      info.second = getOppositeInFather( info.second, info.first->childNr() );
      assert( (info.second >= 0) && (info.second < current.nFaces()) );
      info.first = info.first->father();
    }

    if( info.first )
    {
      assert( info.first->level() == walkLevel_ - 1 );

      const int opposite = info.second;
      for( info.first = info.first->down(); info.first; info.first = info.first->next() )
      {
        info.second = getOppositeInChild( opposite, info.first->childNr() );
        if( info.second != -1 )
          nbStack_.push( info );
      }
    }
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

  //! reset IntersectionIterator to first neighbour
  template<class GridImp>
  inline void ALU2dGridLevelIntersectionIterator<GridImp> :: setFirstItem(const HElementType & elem, int wLevel)
  {
    // empty stack first
    while( ! nbStack_.empty() )
      nbStack_.pop();

    current.setInside( const_cast< HElementType * >( &elem ) );
    current.index_ = -1;
    current.setOutside( 0, -1 );

    walkLevel_ = wLevel;

    assert( current.inside() );

    increment();
  }


  template< class GridImp >
  inline void ALU2dGridLevelIntersectionIterator< GridImp >::setupIntersection ()
  {
    assert( !nbStack_.empty() );

    IntersectionInfo &info = nbStack_.top();
    current.setOutside( info.first, info.second );
    nbStack_.pop();
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
  ALU2dGridLeafIntersectionIterator(const GridImp & grid, const HElementType* el, int wLevel, bool end) :
    ALU2dGridIntersectionBase<GridImp>::ALU2dGridIntersectionBase(grid, el, wLevel, end)
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
  ALU2dGridLeafIntersectionIterator(const GridImp & grid, int wLevel) :
    ALU2dGridIntersectionBase<GridImp>::ALU2dGridIntersectionBase(grid, wLevel)
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

  //! increment iterator
  template<class GridImp>
  inline void ALU2dGridLeafIntersectionIterator<GridImp> :: doIncrement ()
  {
    assert( current.index_ < current.nFaces() );

    this->unsetUp2Date();

    // do we still have neighbours?
    if( !nbStack_.empty() )
      return setupIntersection();

    ++current.index_;
    if( current.index_ >= current.nFaces())
    {
      assert( current.index_ == current.nFaces() );
      return;
    }

    ThinelementType *neighbor = current.inside()->neighbour( current.index_ );
    assert( neighbor );

    if( neighbor->thinis( ThinelementType::bndel_like ) )
    {
      HBndElType *bndel = (HBndElType *)neighbor;
      if( bndel->type() != HBndElType::periodic )
      {
        current.useOutside_ = false;
        return current.setOutside( 0, -1 );
      }

      PeriodicBndElType *bndnb = ((PeriodicBndElType *)bndel)->periodic_nb;
      assert( bndnb && bndnb->neighbour( 0 ) && bndnb->neighbour( 0 )->thinis( ThinelementType::element_like ) );
      current.useOutside_ = !bndnb->leaf();
      if( current.useOutside_ )
      {
        IntersectionInfo info;

        // insert left intersection
        HBndElType *left = bndnb->down();
        assert( left && left->leaf() );
        assert( left->neighbour( 0 ) && left->neighbour( 0 )->thinis( ThinelementType::element_like ) );
        info.first = (HElementType *)left->neighbour( 0 );
        info.second = left->opposite( 0 );
        nbStack_.push( info );

        HBndElType *right = left->next();
        assert( right && right->leaf() );
        assert( right->neighbour( 0 ) && right->neighbour( 0 )->thinis( ThinelementType::element_like ) );
        info.first = (HElementType *)right->neighbour( 0 );
        info.second = right->opposite( 0 );
        nbStack_.push( info );

        setupIntersection();
      }
      else
        current.setOutside( (HElementType *)bndnb->neighbour( 0 ), bndnb->opposite( 0 ) );
    }
    else
    {
      current.useOutside_ = current.inside()->hasHangingNode( current.index_ );
      const int opposite = current.inside()->opposite( current.index_ );
      if( current.useOutside_ )
      {
        IntersectionInfo info;

        // insert left intersection
        ThinelementType *left = current.inside()->getLeftIntersection( current.index_ );
        assert( left && left->thinis( ThinelementType::element_like ) );
        info.first = (HElementType *)left;   // neighbor
        info.second = opposite;              // opposite vertex
        assert( info.first->leaf() );
        nbStack_.push( info );

        // insert right intersection
        ThinelementType *right = current.inside()->getRightIntersection( current.index_ );
        assert( right && right->thinis( ThinelementType::element_like ) );
        info.first = (HElementType *)right;  // neighbor
        info.second = opposite;              // opposite vertex
        assert( info.first->leaf() );
        nbStack_.push( info );

        setupIntersection();
      }
      else
        current.setOutside( (HElementType *)current.inside()->neighbour( current.index_ ), opposite );
    }
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

  //! reset IntersectionIterator to first neighbour
  template<class GridImp>
  inline void ALU2dGridLeafIntersectionIterator<GridImp> :: setFirstItem(const HElementType & elem, int wLevel)
  {
    // empty stack first
    while( ! nbStack_.empty() )
      nbStack_.pop();

    current.setInside( const_cast< HElementType * >( &elem ) );
    current.index_ = -1;
    assert( current.inside() );
    walkLevel_ = wLevel;
    increment();
  }


  template< class GridImp >
  inline void ALU2dGridLeafIntersectionIterator< GridImp >::setupIntersection ()
  {
    assert( !nbStack_.empty() );

    IntersectionInfo &info = nbStack_.top();
    current.setOutside( info.first, info.second );
    nbStack_.pop();
  }

} //end namespace Dune

#endif // #ifndef DUNE_ALU2DGRID_INTERSECTION_IMP_CC
