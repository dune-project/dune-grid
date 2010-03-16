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
    nFaces_(3),
    walkLevel_(wLevel)
  {
    if (!end)
    {
      assert(walkLevel_ >= 0);
      this->setFirstItem(*el,wLevel);
    }
    else
      this->done();
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
    nFaces_(3),
    walkLevel_(wLevel)
  {
    this->done();
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
    nFaces_(org.nFaces_),
    walkLevel_(org.walkLevel_)
  {}

  template<class GridImp>
  inline void
  ALU2dGridIntersectionBase<GridImp> ::
  assign(const ALU2dGridIntersectionBase<GridImp> & org)
  {
    assert( &grid_ == &org.grid_);
    nFaces_    = org.nFaces_;
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
    return ((current.item_ == other.current.item_) && (current.index_ == other.current.index_));
  }


  //! return level of inside() entitiy
  template<class GridImp>
  inline int ALU2dGridIntersectionBase<GridImp> :: level () const
  {
    assert( this->current.item_ );
    return this->current.item_->level();
  }


  //! reset IntersectionIterator to first neighbour
  template<class GridImp>
  template<class EntityType>
  inline void ALU2dGridIntersectionBase<GridImp> ::
  first(const EntityType & en, int wLevel)
  {
    this->setFirstItem(en.getItem(),wLevel);
  #if ALU2DGRID_PARALLEL
    this->checkValid();
  #endif

  }

  //! reset IntersectionIterator to first neighbour
  template< class GridImp >
  inline void ALU2dGridIntersectionBase< GridImp >::setFirstItem(const HElementType & elem, int wLevel)
  {
    current.item_ = const_cast< HElementType * >( &elem );
    assert( current.item_ != 0 );
    walkLevel_ = wLevel;
    current.index_ = 0;
    current.opposite_ = current.item_->opposite( current.index_ );

    unsetUp2Date();
  }


  //! return true if intersection is with boundary
  template<class GridImp>
  inline void ALU2dGridIntersectionBase<GridImp> :: checkValid ()
  {
    // in case of non-valid neighbor set boundary
    if ( this->current.neigh_ &&
         ! this->grid_.rankManager().isValid( this->current.neigh_->getIndex(), All_Partition) )
    {
      this->current.neigh_      = 0;
      this->current.isBoundary_ = true ;
      // this detects a ghost face
      this->current.opposite_   = -222;
    }
  }

  //! return true if intersection is with boundary
  template<class GridImp>
  inline bool ALU2dGridIntersectionBase<GridImp> :: boundary() const
  {
    assert( (this->current.isBoundary_) ? (this->current.neigh_ == 0 ) : 1);
    // make sure that isBoundary is only true on physical boundary
    // in parallel this is not true
    assert( (this->current.opposite_ != -222) ?
            (this->current.isBoundary_ == (this->current.item_->nbbnd(this->current.index_) != 0)) : 1 );
    return this->current.isBoundary_;
  }

  template<class GridImp>
  inline int ALU2dGridIntersectionBase<GridImp> :: boundaryId() const
  {
    int isBoundaryType=0;
    assert(this->current.item_);
    // return boundary id in case of boundary
    // note ALUGrid stores negative values
    if( this->boundary() )
    {
      assert( (this->current.opposite_ != -222) ?
              (this->current.item_->nbbnd(this->current.index_) != 0) : 1);
      isBoundaryType =
  #if ALU2DGRID_PARALLEL
        // if we have ghost face, return 222
        (this->current.opposite_ == -222) ? 222 :
  #endif
        std::abs(this->current.item_->nbbnd(this->current.index_)->type());
    }

    assert( isBoundaryType >= 0 );
    return isBoundaryType;
  }

  template<class GridImp>
  inline size_t ALU2dGridIntersectionBase<GridImp> :: boundarySegmentIndex() const
  {
    // only call this method on boundary intersections
    assert( boundary() );
  #ifdef ALUGRID_VERTEX_PROJECTION
    return this->current.item_->nbbnd(this->current.index_)->segmentIndex();
  #else
    derr << "Method available in any version of ALUGrid > 1.14 \n";
    return 0;
  #endif
  }

  //! return true if intersection is with neighbor on this level
  template<class GridImp>
  inline bool ALU2dGridIntersectionBase<GridImp> :: neighbor () const {
    return (this->current.neigh_ && !boundary() );
  }

  //! return EntityPointer to the Entity on the inside of this intersection.
  template<class GridImp>
  inline typename ALU2dGridIntersectionBase<GridImp> :: EntityPointer
  ALU2dGridIntersectionBase<GridImp> :: inside() const {
    assert(this->current.item_);
    return EntityPointerImp(grid_, *(this->current.item_), -1, walkLevel_);
  }

  template< class GridImp >
  inline void ALU2dGridIntersectionBase<GridImp>::done ()
  {
    current.item_ = 0;
    current.neigh_ = 0;
    current.index_= nFaces_;
  }

  //! return EntityPointer to the Entity on the outside of this intersection.
  template<class GridImp>
  inline typename ALU2dGridIntersectionBase<GridImp> :: EntityPointer
  ALU2dGridIntersectionBase<GridImp> :: outside() const {
    assert(!this->boundary());
    assert( this->current.neigh_ );
    return EntityPointerImp(grid_, *(this->current.neigh_),-1, walkLevel_);
  }

  //! local number of codim 1 entity in self where intersection is contained in
  template<class GridImp>
  inline int ALU2dGridIntersectionBase<GridImp>::indexInInside () const
  {
    return 2 - this->current.index_;
  }

  //! local number of codim 1 entity in neighbor where intersection is contained in
  template<class GridImp>
  inline int ALU2dGridIntersectionBase<GridImp>::indexInOutside () const
  {
    return 2 - this->current.opposite_;
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
    return (1 + this->current.index_ + this->current.opposite_) % 2;
  }

  template<class GridImp>
  inline typename ALU2dGridIntersectionBase<GridImp>::NormalType &
  ALU2dGridIntersectionBase<GridImp> :: outerNormal (const FieldVector<alu2d_ctype, dim-1>& local) const
  {
    assert(this->current.item_ != 0);
    typedef double (&normal_t)[dimworld];

    this->current.item_->outernormal(  this->current.index_,  ((normal_t) (&outerNormal_)[0]) );
    if( this->current.useOutside_ )
      outerNormal_ *= 0.5;
    return outerNormal_;
  }

  template<class GridImp>
  inline typename ALU2dGridIntersectionBase<GridImp>::NormalType &
  ALU2dGridIntersectionBase<GridImp> :: integrationOuterNormal (const FieldVector<alu2d_ctype, dim-1>& local) const {
    return this->outerNormal(local);
  }

  template<class GridImp>
  inline typename ALU2dGridIntersectionBase<GridImp>::NormalType
  ALU2dGridIntersectionBase<GridImp> :: unitOuterNormal (const FieldVector<alu2d_ctype, dim-1>& local) const
  {
    NormalType unitNormal( this->outerNormal(local) );
    unitNormal *= (1.0/unitNormal.two_norm());
    return unitNormal;
  }


  template<class GridImp>
  inline const typename ALU2dGridIntersectionBase<GridImp>::LocalGeometry&
  ALU2dGridIntersectionBase< GridImp >::geometryInInside () const
  {
    assert(this->current.item_ != 0);

    // only in non-conform situation we use default method
    if( this->current.useOutside_ )
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
      // parameters are face and twist
      return localGeomStorage_.localGeom(  this->current.index_,
                                           (this->current.index_ % 2) );
    }
  }

  template<class GridImp>
  inline const typename ALU2dGridIntersectionBase<GridImp>::LocalGeometry&
  ALU2dGridIntersectionBase< GridImp >::geometryInOutside () const
  {
    assert( this->current.item_  );
    assert( this->current.neigh_ );

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
       return localGeomStorage_.localGeom(  this->current.opposite_,
                                          1 - (this->current.index_ % 2) );
       }
     */
  }

  template<class GridImp>
  inline const typename ALU2dGridIntersectionBase<GridImp>::Geometry&
  ALU2dGridIntersectionBase<GridImp>::geometry () const
  {
    assert(this->current.item_ != 0);

    if( ! this->grid_.getRealImplementation(intersectionGlobal_).up2Date() )
    {
      if( this->current.useOutside_ )
        this->grid_.getRealImplementation(intersectionGlobal_).
        buildGeom(*(this->current.neigh_), this->current.opposite_);
      else
        this->grid_.getRealImplementation(intersectionGlobal_).
        buildGeom(*(this->current.item_), this->current.index_);
    }

    assert(this->grid_.getRealImplementation(intersectionGlobal_).up2Date());
    return intersectionGlobal_;
  }


  template< class GridImp >
  inline GeometryType ALU2dGridIntersectionBase< GridImp >::type () const
  {
    return GeometryType( GeometryType::simplex, dimension-1 );
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
    {
      this-> done();
    }
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
    neighbourStack_ = org.neighbourStack_;
  }

  //! The copy constructor
  template<class GridImp>
  inline void
  ALU2dGridLevelIntersectionIterator<GridImp> ::
  assign(const ALU2dGridLevelIntersectionIterator<GridImp> & org) {
    ALU2dGridIntersectionBase<GridImp>:: ALU2dGridIntersectionBase::assign(org);
    neighbourStack_ = org.neighbourStack_;
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
    assert( current.index_ < nFaces_ );
#if 0
    // this cannot happen
    if( current.index_ >= nFaces_ )
    {
      this->done();
      return;
    }
#endif

    this->unsetUp2Date();

    if( neighbourStack_.empty() )
    {
      ++current.index_;
      if( current.index_ >= nFaces_ )
      {
        this->done();
        return;
      }

      addNeighboursToStack();
      // if more then one element in stack we have non-conform intersection
      this->current.useOutside_ = (neighbourStack_.size() > 1);

      if( neighbourStack_.empty() )
      {
        current.neigh_ = 0;
        current.opposite_ = current.item_->opposite( current.index_ );
        // we have only boundary if intersection is on domain boundary
        current.isBoundary_  = (current.item_->nbbnd( current.index_ ) != 0);
        return;
      }
    }

    setupIntersection();

    assert( !current.isBoundary_ );
    assert( !current.neigh_ || (current.neigh_->level() == walkLevel_) );
  }


  template< class GridImp >
  inline void
  ALU2dGridLevelIntersectionIterator< GridImp >::addNeighboursToStack ()
  {
    assert( current.index_ < nFaces_ );

    IntersectionInfo dummy;
    dummy.first = current.item_->nbel( current.index_ );
    if( dummy.first == 0 )
      return;

    dummy.second.first = current.item_->opposite( current.index_ );
    dummy.second.second = false;

    while( dummy.first->level() > walkLevel_ )
    {
      dummy.second.first = getOppositeInFather( dummy.second.first, dummy.first->childNr() );
      assert( (dummy.second.first >= 0) && (dummy.second.first) < nFaces_); // parsing error???
      dummy.first = dummy.first->father();
    }

    if( dummy.first->level() < walkLevel_ )
    {
      while( dummy.first != 0 )
      {
        const int level = dummy.first->level();
        if( level >= walkLevel_ - 1 )
          break;

        dummy.second.first = getOppositeInFather( dummy.second.first, dummy.first->childNr() );
        assert( (dummy.second.first >= 0) && (dummy.second.first) < nFaces_); // parsing error???
        dummy.first = dummy.first->father();
      }

      if( dummy.first != 0 )
      {
        assert( dummy.first->level() == walkLevel_ - 1 );

        for( HElementType *tmp = dummy.first->down(); tmp != 0; tmp = tmp->next() )
        {
          int tmpOpposite = getOppositeInChild( dummy.second.first, tmp->childNr() );
          if( tmpOpposite != -1 )
            neighbourStack_.push( IntersectionInfo( tmp, std::pair< int, bool >( tmpOpposite, false ) ) );
        }
      }
    }
    else
      neighbourStack_.push( dummy );
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
    while( ! neighbourStack_.empty() )
    {
      neighbourStack_.pop();
    }

    this->current.item_ = const_cast<HElementType *> (&elem);
    this->current.index_ = -1;
    this->current.neigh_ = 0;
    this->current.opposite_= -1;


    this->walkLevel_ = wLevel;

    assert(this->current. item_ );

    increment();
  }


  template< class GridImp >
  inline void ALU2dGridLevelIntersectionIterator< GridImp >::setupIntersection ()
  {
    assert( !neighbourStack_.empty() );

    IntersectionInfo &info = neighbourStack_.top();
    current.neigh_ = info.first;
    current.opposite_ = info.second.first;
    current.isBoundary_ = info.second.second;
    neighbourStack_.pop();
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
    {
      this->done();
    }
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
    assert( current.index_ < nFaces_ );
#if 0
    // this cannot happen
    if( current.index_ >= nFaces_ )
    {
      this->done();
      return ;
    }
#endif

    this->unsetUp2Date();

    // non conform case and we still have neighbours
    if( current.item_->hashvtx( current.index_ ) && !nbStack_.empty() )
    {
      setupIntersection();
      return;
    }

    ++current.index_;
    if( current.index_ >= nFaces_ )
    {
      this->done();
      return;
    }

    current.useOutside_ = current.item_->hasHangingNode( current.index_ );
    // non-conform case
    if( current.useOutside_ )
    {
      // stack should be empty here
      assert( nbStack_.empty() );

      IntersectionInfo dummy;

      // insert left intersection
      dummy.first = current.item_->getLeftIntersection( current.index_ );  // neighbor
      dummy.second.first = current.item_->opposite( current.index_ );      // opposite vertex
      dummy.second.second = (dummy.first ==0);                             // is boundary
      nbStack_.push( dummy );

      // insert right intersection
      dummy.first = current.item_->getRightIntersection( current.index_ ); // neighbor
      dummy.second.first = current.item_->opposite( current.index_ );      // opposite vertex
      dummy.second.second = (dummy.first == 0);                            // is boundary
      nbStack_.push( dummy );

      setupIntersection();
    }
    //conform case
    else
    {
      current.neigh_ = current.item_->nbel( current.index_ );
      current.isBoundary_ = (current.neigh_ == 0);
      current.opposite_ = (current.isBoundary_ ? -1 : current.item_->opposite( current.index_ ));
    }

#if 0
    // this cannot happen
    if( current.index_ >= nFaces_ )
      this->done();
#endif
  }


  //! reset IntersectionIterator to first neighbour
  template<class GridImp>
  template<class EntityType>
  inline void ALU2dGridLeafIntersectionIterator<GridImp> ::
  first(const EntityType & en, int wLevel)
  {
    // if called on non-leaf, just return end iterator
    if(! en.isLeaf())
    {
      this->done();
      return ;
    }

    setFirstItem(en.getItem(),wLevel);
  }

  //! reset IntersectionIterator to first neighbour
  template<class GridImp>
  inline void ALU2dGridLeafIntersectionIterator<GridImp> :: setFirstItem(const HElementType & elem, int wLevel)
  {
    // empty stack first
    while( ! nbStack_.empty() )
    {
      nbStack_.pop();
    }

    this->current.item_ = const_cast<HElementType *> (&elem);
    this->current.index_ = -1;
    assert(this->current.item_ );
    this->walkLevel_ = wLevel;
    increment();
  }


  template< class GridImp >
  inline void ALU2dGridLeafIntersectionIterator< GridImp >::setupIntersection ()
  {
    assert( !nbStack_.empty() );

    IntersectionInfo &info = nbStack_.top();
    current.neigh_ = static_cast< HElementType * >( info.first );
    current.opposite_ = info.second.first;
    current.isBoundary_ = info.second.second;
    nbStack_.pop();
  }

} //end namespace Dune

#endif // #ifndef DUNE_ALU2DGRID_INTERSECTION_IMP_CC
