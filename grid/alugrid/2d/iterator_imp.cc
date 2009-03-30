// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "geometry.hh"
#include "entity.hh"
#include "grid.hh"
//#include "faceutility.hh"
#include <stack>
#include <utility>

#ifndef DUNE_ALU2DGRID_ITERATOR_IMP_CC
#define DUNE_ALU2DGRID_ITERATOR_IMP_CC


namespace Dune {


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
    intersectionSelfLocal_(GeometryImp()),
    intersectionNeighborLocal_(GeometryImp()),
    grid_(grid),
    nFaces_(3),
    walkLevel_(wLevel),
    generatedGlobalGeometry_(false),
    generatedLocalGeometries_(false),
    done_(end)
  {
    if (!end)
    {
      assert(walkLevel_ >= 0);
      this->setFirstItem(*el,wLevel);
    }
    else
    {
      this->done();
    }
  }

  //constructor for end iterator
  template<class GridImp>
  inline ALU2dGridIntersectionBase<GridImp> ::
  ALU2dGridIntersectionBase(const GridImp & grid, int wLevel) :
    intersectionGlobal_(GeometryImp()),
    intersectionSelfLocal_(GeometryImp()),
    intersectionNeighborLocal_(GeometryImp()),
    grid_(grid),
    nFaces_(3),
    walkLevel_(wLevel),
    generatedGlobalGeometry_(false),
    generatedLocalGeometries_(false),
    done_(true)
  {
    this->done();
  }


  //! The copy constructor
  template<class GridImp>
  inline ALU2dGridIntersectionBase<GridImp> ::
  ALU2dGridIntersectionBase(const ALU2dGridIntersectionBase<GridImp> & org) :
    current(org.current),
    intersectionGlobal_(GeometryImp()),
    intersectionSelfLocal_(GeometryImp()),
    intersectionNeighborLocal_(GeometryImp()),
    grid_(org.grid_),
    nFaces_(org.nFaces_),
    walkLevel_(org.walkLevel_),
    generatedGlobalGeometry_(false),
    generatedLocalGeometries_(false),
    done_(org.done_)
  {}

  template<class GridImp>
  inline void
  ALU2dGridIntersectionBase<GridImp> ::
  assign(const ALU2dGridIntersectionBase<GridImp> & org)
  {
    assert( &grid_ == &org.grid_);
    nFaces_    = org.nFaces_;
    walkLevel_ = org.walkLevel_;
    generatedGlobalGeometry_ = false;
    generatedLocalGeometries_ = false;
    done_ = org.done_;
    current = org.current;

    // unset geometry information
    this->grid_.getRealImplementation(intersectionGlobal_).unsetUp2Date();
    this->grid_.getRealImplementation(intersectionSelfLocal_).unsetUp2Date();
    this->grid_.getRealImplementation(intersectionNeighborLocal_).unsetUp2Date();
  }

  //! check whether entities are the same or whether iterator is done
  template<class GridImp>
  inline bool ALU2dGridIntersectionBase<GridImp> ::
  equals (const ALU2dGridIntersectionBase<GridImp> & i) const
  {
    return ((this->current.item_ == i.current.item_) &&
            (this->done_ == i.done_));
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
  template<class GridImp>
  inline void ALU2dGridIntersectionBase<GridImp> :: setFirstItem(const HElementType & elem, int wLevel)
  {
    this->current.item_ = const_cast<HElementType *> (&elem);
    assert( this->current.item_ );
    walkLevel_ = wLevel;
    done_ = false;
    this->current.index_ = 0;
    this->current.opposite_ = this->current.item_->opposite(this->current.index_);

    this->grid_.getRealImplementation(intersectionGlobal_).unsetUp2Date();
    this->grid_.getRealImplementation(intersectionSelfLocal_).unsetUp2Date();
    this->grid_.getRealImplementation(intersectionNeighborLocal_).unsetUp2Date();
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

  template<class GridImp>
  inline void ALU2dGridIntersectionBase<GridImp> :: done() {
    done_ = true;
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
  inline int ALU2dGridIntersectionBase<GridImp> :: numberInInside () const {
    return this->current.index_;
  }

  //! local number of codim 1 entity in neighbor where intersection is contained in
  template<class GridImp>
  inline int ALU2dGridIntersectionBase<GridImp> :: numberInOutside () const {
    return this->current.opposite_;
  }

  template< class GridImp >
  inline int ALU2dGridIntersectionBase< GridImp >::twistInSelf () const
  {
    //return 0;
    return 0;
  }

  template< class GridImp >
  inline int ALU2dGridIntersectionBase< GridImp >::twistInNeighbor () const
  {
    //return 1;
    return 1;
  }

  template<class GridImp>
  inline typename ALU2dGridIntersectionBase<GridImp>::NormalType &
  ALU2dGridIntersectionBase<GridImp> :: outerNormal (const FieldVector<alu2d_ctype, dim-1>& local) const
  {
    assert(this->current.item_ != 0);
    typedef double (&normal_t)[2];

    if (this->current.isNotConform_ )
    {
      this->current.neigh_->outernormal( numberInOutside(), ((normal_t) (&outerNormal_)[0]) );
      outerNormal_ *= -1.0;
      return outerNormal_;
    }

    this->current.item_->outernormal(this->current.index_,  ((normal_t) (&outerNormal_)[0]) );
    return outerNormal_;
  }

  template<class GridImp>
  inline typename ALU2dGridIntersectionBase<GridImp>::NormalType &
  ALU2dGridIntersectionBase<GridImp> :: integrationOuterNormal (const FieldVector<alu2d_ctype, dim-1>& local) const {
    return this->outerNormal(local);
  }

  template<class GridImp>
  inline typename ALU2dGridIntersectionBase<GridImp>::NormalType &
  ALU2dGridIntersectionBase<GridImp> :: unitOuterNormal (const FieldVector<alu2d_ctype, dim-1>& local) const {
    unitOuterNormal_ = this->outerNormal(local);
    unitOuterNormal_ *= (1.0/unitOuterNormal_.two_norm());
    return unitOuterNormal_;
  }


  template<class GridImp>
  inline const typename ALU2dGridIntersectionBase<GridImp>::LocalGeometry&
  ALU2dGridIntersectionBase< GridImp >::geometryInInside () const
  {
    assert(this->current.item_ != 0);

    if( ! this->grid_.getRealImplementation(intersectionSelfLocal_).up2Date() )
    {
      // only in non-conform situation we use default method
      if( this->current.isNotConform_ )
      {
        this->grid_.getRealImplementation(intersectionSelfLocal_).
        builtLocalGeom( inside()->geometry() , geometry() );
      }
      else
      {
        this->grid_.getRealImplementation(intersectionSelfLocal_).
        builtLocalGeom( numberInInside() , 0 );
      }
    }

    assert(this->grid_.getRealImplementation(intersectionSelfLocal_).up2Date());
    return intersectionSelfLocal_;
  }

  template<class GridImp>
  inline const typename ALU2dGridIntersectionBase<GridImp>::LocalGeometry&
  ALU2dGridIntersectionBase< GridImp >::geometryInOutside () const
  {
    assert(this->current.item_ != 0);
    assert(this->current.neigh_ != 0);

    if( ! this->grid_.getRealImplementation(intersectionNeighborLocal_).up2Date() )
    {
      // we don't know here wether we have non-conform or conform situation on the neighbor
      this->grid_.getRealImplementation(intersectionNeighborLocal_).
      builtLocalGeom( outside()->geometry(), geometry() );
    }
    assert(this->grid_.getRealImplementation(intersectionNeighborLocal_).up2Date());
    return intersectionNeighborLocal_;
  }

  template<class GridImp>
  inline const typename ALU2dGridIntersectionBase<GridImp>::Geometry&
  ALU2dGridIntersectionBase<GridImp>::geometry () const
  {
    assert(this->current.item_ != 0);

    if( ! this->grid_.getRealImplementation(intersectionGlobal_).up2Date() )
    {
      if( this->current.isNotConform_ )
        this->grid_.getRealImplementation(intersectionGlobal_).builtGeom(*(this->current.neigh_), this->current.opposite_);
      else
        this->grid_.getRealImplementation(intersectionGlobal_).builtGeom(*(this->current.item_), this->current.index_);
    }

    assert(this->grid_.getRealImplementation(intersectionGlobal_).up2Date());
    return intersectionGlobal_;
  }


  template< class GridImp >
  inline GeometryType ALU2dGridIntersectionBase< GridImp >::type () const
  {
    return GeometryType( GeometryType::simplex, dimension-1 );
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
    if (this->current.index_ >= this->nFaces_) {
      this->done();
      return;
    }

    this->grid_.getRealImplementation(this->intersectionGlobal_).unsetUp2Date();
    this->grid_.getRealImplementation(this->intersectionSelfLocal_).unsetUp2Date();
    this->grid_.getRealImplementation(this->intersectionNeighborLocal_).unsetUp2Date();

    if (neighbourStack_.empty())
    {
      ++this->current.index_;
      if (this->current.index_ >= this->nFaces_)
      {
        this->done();
        return ;
      }

      if(this->current.item_->hasHangingNode(this->current.index_))
        this->current.isNotConform_ = true;
      else
        this->current.isNotConform_ = false;

      addNeighboursToStack();
      if (neighbourStack_.empty())
      {
        this->current.neigh_ = 0;
        this->current.opposite_ = this->current.item_->opposite(this->current.index_);
        // we have only boundary if intersection is on domain boundary
        this->current.isBoundary_  = (this->current.item_->nbbnd(this->current.index_) != 0);
        return;
      }
    }

    IntersectionInfo& info    = neighbourStack_.top();
    this->current.neigh_      = info.first;
    this->current.opposite_   = info.second.first;
    this->current.isBoundary_ = info.second.second;
    assert( this->current.isBoundary_ == false);
    neighbourStack_.pop();

    assert( (this->current.neigh_) ? (this->current.neigh_->level()==this->walkLevel_) : 1);
    return;
  }

  template<class GridImp>
  inline void ALU2dGridLevelIntersectionIterator<GridImp> ::
  addNeighboursToStack ()
  {
    assert (this->current.index_ < 3);
    IntersectionInfo dummy;

    dummy.first = this->current.item_->nbel(this->current.index_);
    if(dummy.first == 0)
    {
      return ;
    }

    dummy.second.first = this->current.item_->opposite(this->current.index_);
    dummy.second.second = false;

    if (dummy.first->level() == this->walkLevel_) {
      neighbourStack_.push(dummy);
    }
    else if (dummy.first->level() > this->walkLevel_)
    {
      while (dummy.first->level() > this->walkLevel_)
      {
        dummy.second.first =
          getOppositeInFather(dummy.second.first, dummy.first->childNr());
        assert(dummy.second.first >= 0 && dummy.second.first < 3);
        dummy.first = dummy.first->father();
      }
      assert(dummy.first->level()==this->walkLevel_);
      neighbourStack_.push(dummy);
    }
    else
    {
      while ( dummy.first )
      {
        const int lev = dummy.first->level();
        if(lev >= this->walkLevel_ - 1) break;

        dummy.second.first =
          getOppositeInFather(dummy.second.first, dummy.first->childNr() );
        assert(dummy.second.first >= 0 && dummy.second.first < 3);

        // get next
        dummy.first = dummy.first->father();
      }

      if (dummy.first)
      {
        assert(dummy.first->level() == this->walkLevel_ - 1);

        HElementType * tmp = dummy.first->down();
        while ( tmp )
        {
          int tmpOpposite = getOppositeInChild(dummy.second.first, tmp->childNr() );
          if (tmpOpposite != -1)
          {
            IntersectionInfo dummy2(tmp, std::pair<int, bool> (tmpOpposite, false));
            neighbourStack_.push(dummy2);
          }
          tmp = tmp->next();
        }
      }
    }
    // if more then one element in stack we have non-conform intersection
    this->current.isNotConform_ = (neighbourStack_.size() > 1);
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
    this->current.index_ = 0;
    this->current.neigh_ = elem.nbel(0);
    this->current.opposite_= elem.opposite(0);

    assert(this->current.opposite_ >= 0 && this->current.opposite_ < 3);
    assert(this->current. item_ );

    this->walkLevel_ = wLevel;
    this->done_ = false;

    if(this->current.item_->hasHangingNode(0))
      this->current.isNotConform_ = true;
    else
      this->current.isNotConform_ = false;

    if(!this->current.neigh_ )
    {
      this->current.isBoundary_ = true ;
      this->current.opposite_ = -1;
    }
    else {
      this->current.isBoundary_ = false;
      if (this->current.neigh_->level() != this->walkLevel_)
      {
        // index is increased in increment again
        this->current.index_ = -1;
        this->current.neigh_ = 0;
        increment();
      }
    }
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
    if (this->current.index_ >= this->nFaces_)
    {
      this->done();
      return ;
    }

    this->grid_.getRealImplementation(this->intersectionGlobal_).unsetUp2Date();
    this->grid_.getRealImplementation(this->intersectionSelfLocal_).unsetUp2Date();
    this->grid_.getRealImplementation(this->intersectionNeighborLocal_).unsetUp2Date();

    // non conform case and we still have neighbours
    if(this->current.item_->hashvtx(this->current.index_) && !nbStack_.empty())
    {
      IntersectionInfo& info = nbStack_.top();
      this->current.neigh_      = static_cast<HElementType *>(info.first);
      this->current.opposite_   = info.second.first;
      this->current.isBoundary_ = info.second.second;
      nbStack_.pop();
      return;
    }

    ++this->current.index_;
    if (this->current.index_ >= this->nFaces_)
    {
      this->done();
      return;
    }

    // non-conform case
    if (this->current.item_->hasHangingNode(this->current.index_))
    {
      this->current.isNotConform_ = true;
      // stack should be empty here
      assert( nbStack_.empty() );

      IntersectionInfo dummy;

      // insert left intersection
      dummy.first = this->current.item_->getLeftIntersection(this->current.index_);
      dummy.second.first = this->current.item_->opposite(this->current.index_);
      dummy.second.second = (dummy.first==0) ? true : false ;
      nbStack_.push(dummy);

      // insert right intersection
      dummy.first = this->current.item_->getRightIntersection(this->current.index_);
      dummy.second.first = this->current.item_->opposite(this->current.index_);
      dummy.second.second = (dummy.first==0) ? true : false ;
      nbStack_.push(dummy);

      // nbStack should contain some elements here
      assert( ! nbStack_.empty() );

      IntersectionInfo& info    = nbStack_.top();
      this->current.neigh_      = static_cast<HElementType *>(info.first);
      this->current.opposite_   = info.second.first;
      this->current.isBoundary_ = info.second.second;
      nbStack_.pop();
    }
    //conform case
    else
    {
      this->current.isNotConform_ = false;
      this->current.neigh_ = this->current.item_->nbel(this->current.index_);
      if (this->current.neigh_ == 0)
      {
        this->current.isBoundary_ = true;
        this->current.opposite_ = -1;
      }
      else
      {
        this->current.isBoundary_ = false;
        this->current.opposite_ = this->current.item_->opposite(this->current.index_);
      }
    }

    if (this->current.index_ >= this->nFaces_)
      this->done();

    return;
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
    this->done_ = false;
    increment();
  }


  //********************************************************************
  //
  //  --ALU2dGridLeafIterator
  //  --LeafIterator
  //
  //********************************************************************

  //! constructor
  template<int cdim, PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLeafIterator<cdim, pitype, GridImp> ::
  ALU2dGridLeafIterator(const GridImp & grid, bool end) :
    EntityPointerType (grid),
    endIter_(end),
    level_(-1),
    elem_(0),
    iter_(),
    marker_(grid.getLeafMarker())
  {
    if(!end)
    {
      // update marker Vector
      if( (cdim == 2) && (! marker_.up2Date()) ) marker_.update(grid);

      iter_ = IteratorType(grid.myGrid());
      iter_->first();

      if((!iter_->done()))
      {
        elem_ = &(iter_->getitem());

#if ALU2DGRID_PARALLEL
        bool valid = (cdim == 0) ?
                     (this->grid_.rankManager().isValid( elem_->getIndex(), pitype )) :
                     marker_.isValidVertex( elem_->getIndex() );
#else
        bool valid = true;
#endif

        // if we found valid item, call update
        if( valid )
        {
          this->updateEntityPointer(elem_, -1,
                                    GetLevel<ElementType,LeafMarkerVectorType,cdim>::level(*elem_,marker_));
        }
        else
        {
          increment();
        }
      }
    }
    else
    {
      endIter_ = true;
      this->done();
    }
  }

  //! copy constructor
  template<int cdim, PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLeafIterator<cdim, pitype, GridImp> ::
  ALU2dGridLeafIterator(const ALU2dGridLeafIterator<cdim,pitype,GridImp> & org)
    : EntityPointerType (org)
      , endIter_( org.endIter_ )
      , level_( org.level_ )
      , elem_(org.elem_)
      , iter_ ( org.iter_ )
      , marker_(org.marker_)
  {}

  //! assignment
  template<int cdim, PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLeafIterator<cdim, pitype, GridImp> &
  ALU2dGridLeafIterator<cdim, pitype, GridImp> ::
  operator = (const ThisType & org)
  {
    EntityPointerType :: operator = (org);
    endIter_ =  org.endIter_ ;
    level_   =  org.level_;
    elem_    =  org.elem_;
    iter_    =  org.iter_;

    // there is only one reference, so we don't copy
    assert(&marker_ == &org.marker_);
    return *this;
  }


  //! prefix increment
  template<int cdim, PartitionIteratorType pitype, class GridImp>
  inline void ALU2dGridLeafIterator<cdim, pitype, GridImp> :: increment ()
  {
    if(endIter_) return ;

    // go to next item
    iter_->next();

    if(iter_->done()) {
      endIter_ = true;
      this->done();
      return ;
    }

    // get element pointer
    elem_ = &(iter_->getitem());

#if ALU2DGRID_PARALLEL
    bool valid = (cdim == 0) ?
                 (this->grid_.rankManager().isValid( elem_->getIndex(), pitype )) :
                 marker_.isValidVertex( elem_->getIndex() );

    while ( ! valid )
    {
      // go to next item
      iter_->next();

      if(iter_->done())
      {
        endIter_ = true;
        this->done();
        return ;
      }

      // get element pointer
      elem_ = &(iter_->getitem());
      valid = (cdim == 0) ?
              (this->grid_.rankManager().isValid( elem_->getIndex(), pitype )) :
              marker_.isValidVertex( elem_->getIndex() );
    }
#endif

    this->updateEntityPointer(elem_, -1,
                              GetLevel<ElementType,LeafMarkerVectorType,cdim>::level(*elem_,marker_));
  }

  //********************************************************************
  //
  //  -- ALU2dGridLeafIterator
  //  -- LeafIterator
  //  -- specialized for codim=1
  //
  //********************************************************************

  //! constructor
  template<PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLeafIterator<1, pitype, GridImp> ::
  ALU2dGridLeafIterator(const GridImp & grid, bool end) :
    EntityPointerType (grid),
    endIter_(end),
    level_(-1),
    face_(0),
    elem_(0),
    iter_(),
    marker_(grid.getLeafMarker())
  {
    if(!end)
    {
      // update marker Vector
      if( ! marker_.up2Date() ) marker_.update(grid);

      iter_ = IteratorType(grid.myGrid());
      iter_->first();
      if((!iter_->done()))
      {
        elem_ = &(iter_->getitem());
        this->updateEntityPointer(elem_, face_, elem_->level());
        increment();
      }
    }
    else
    {
      endIter_ = true;
      this->done();
    }
  }

  //! copy constructor
  template<PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLeafIterator<1, pitype, GridImp> ::
  ALU2dGridLeafIterator(const ALU2dGridLeafIterator<1,pitype,GridImp> & org)
    : EntityPointerType (org)
      , endIter_( org.endIter_ )
      , level_( org.level_ )
      , face_(org.face_)
      , elem_(org.elem_)
      , iter_ ( org.iter_ )
      , marker_ (org.marker_)
  {}

  //! assignment
  template<PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLeafIterator<1, pitype, GridImp> &
  ALU2dGridLeafIterator<1, pitype, GridImp> ::
  operator = (const ThisType & org)
  {
    EntityPointerType :: operator = (org);
    endIter_ =  org.endIter_ ;
    level_   =  org.level_;
    face_    =  org.face_;
    elem_    =  org.elem_;
    iter_    =  org.iter_;

    // there is only one reference, so we don't copy
    assert(&marker_ == &org.marker_);
    return *this;
  }


  //! prefix increment
  template<PartitionIteratorType pitype, class GridImp>
  inline void ALU2dGridLeafIterator<1, pitype, GridImp> :: increment() {

    IteratorType & iter = iter_;
    if(iter->done()) {
      face_= -1;
      return ;
    }

    int goNext = goNextElement();
    if (goNext) {
      assert(face_==3);

      // go next element
      iter->next();
      if(iter->done())
      {
        endIter_=true;
        face_= -1;
        this->done();
        return ;
      }
      // get new element
      elem_ = &(iter->getitem());
#if ALU2DGRID_PARALLEL
      while ( ! this->grid_.rankManager().isValid( elem_->getIndex(), pitype ) )
      {
        // go to next item
        iter->next();
        if(iter->done())
        {
          endIter_=true;
          face_= -1;
          this->done();
          return ;
        }
        elem_ = &(iter->getitem());
      }
#endif

      face_=0;
      // whatever this is good for
      this->updateEntityPointer(elem_, face_, elem_->level());
      increment();
    }
    else
    {
      if(iter->done())
      {
        endIter_=true;
        face_= -1;
        this->done();
        return ;
      }
      elem_ = &(iter->getitem());
      this->updateEntityPointer(elem_, face_, elem_->level());
      ++face_;
    }
    return;
  }

  template<PartitionIteratorType pitype, class GridImp>
  inline int ALU2dGridLeafIterator<1, pitype, GridImp> :: goNextElement() {
    assert(face_>=0);
    int elIdx = this->item_->getIndex();

    while (face_ < 3) {
      int idx = this->item_->edge_idx(face_);
      // check if face is visited on this element
      if(marker_.isOnElement(elIdx,idx,1))
        return 0;
      else
        ++face_;
    }
    return 1;
  }


  //********************************************************************
  //
  //  --ALU2dLevelLeafIterator
  //  --LevelIterator, specialized for cd=0
  //
  //********************************************************************

  //! constructor for cd=0
  template<PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLevelIterator<0, pitype, GridImp> ::
  ALU2dGridLevelIterator(const GridImp & grid, int level, bool end) :
    EntityPointerType (grid),
    endIter_(end),
    level_(level),
    iter_()
  {
    if(!end)
    {
      iter_ = IteratorType(grid.myGrid(), level_);
      iter_->first();
      if((!iter_->done()))
      {
        item_ = &(iter_->getitem());
        this->updateEntityPointer(item_, -1 , level_);
#if ALU2DGRID_PARALLEL
        if ( ! this->grid_.rankManager().isValid( item_->getIndex(), pitype ) )
        {
          increment();
        }
#endif
      }

    }
    else
    {
      endIter_ = true;
      this->done();
    }
  }

  //! copy constructor for cd=1
  template<PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLevelIterator<0, pitype, GridImp> ::
  ALU2dGridLevelIterator(const ALU2dGridLevelIterator<0,pitype,GridImp> & org)
    : EntityPointerType (org)
      , endIter_( org.endIter_ )
      , level_( org.level_ )
      , item_(org.item_)
      , iter_ (org.iter_)
  {}

  //! assignment
  template<PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLevelIterator<0, pitype, GridImp> &
  ALU2dGridLevelIterator<0, pitype, GridImp> ::
  operator = (const ThisType & org)
  {
    EntityPointerType :: operator = (org);
    endIter_ =  org.endIter_ ;
    level_   =  org.level_;
    item_    =  org.item_;
    iter_    =  org.iter_;
    return *this;
  }

  //! prefix increment
  template<PartitionIteratorType pitype, class GridImp>
  inline void ALU2dGridLevelIterator<0, pitype, GridImp> :: increment () {

    if(endIter_) return ;

    IteratorType & iter = iter_;

    iter->next();
    if(iter->done()) {
      endIter_ = true;
      this->done();
      return ;
    }
    item_ = &iter->getitem();

#if ALU2DGRID_PARALLEL
    while ( ! this->grid_.rankManager().isValid( item_->getIndex(), pitype ) )
    {
      // go to next item
      iter_->next();

      if(iter_->done()) {
        endIter_ = true;
        this->done();
        return ;
      }

      // get element pointer
      item_ = &(iter_->getitem());
    }
#endif

    this->updateEntityPointer(item_, -1 , level_);
    return;
  }


  //********************************************************************
  //
  //  --ALU2dLevelLeafIterator
  //  --LevelIterator, specialized for cd=1
  //
  //********************************************************************

  //! constructor for cd=1
  template<PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLevelIterator<1, pitype, GridImp> ::
  ALU2dGridLevelIterator(const GridImp & grid, int level, bool end) :
    EntityPointerType (grid),
    endIter_(end),
    level_(level),
    myFace_(0),
    iter_(),
    marker_(& grid.getMarkerVector(level))
  {
    if(!end)
    {
      // update marker Vector if necessary
      if( ! marker().up2Date() ) marker().update(grid,level_);

      iter_ = IteratorType(grid.myGrid(), level_);
      iter_->first();

      if((!iter_->done()))
      {
        elem_ = &(iter_->getitem());
        this->updateEntityPointer(elem_, myFace_, level_);
        increment();
      }
    }
    else
    {
      endIter_ = true;
      this->done();
    }
  }

  //! copy constructor for cd=1
  template<PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLevelIterator<1, pitype, GridImp> ::
  ALU2dGridLevelIterator(const ALU2dGridLevelIterator<1,pitype,GridImp> & org)
    : EntityPointerType (org)
      , endIter_( org.endIter_ )
      , level_( org.level_ )
      , myFace_(org.myFace_)
      , item_(org.item_)
      , elem_(org.elem_)
      , iter_ ( org.iter_ )
      , marker_(org.marker_)
  {}

  //! assignment
  template<PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLevelIterator<1, pitype, GridImp> &
  ALU2dGridLevelIterator<1, pitype, GridImp> ::
  operator = (const ThisType & org)
  {
    EntityPointerType :: operator = (org);
    endIter_ = org.endIter_ ;
    level_   = org.level_;
    myFace_  = org.myFace_;
    item_    = org.item_;
    elem_    = org.elem_;
    iter_    = org.iter_;
    marker_  = org.marker_;

    assert(marker_ == org.marker_);
    return *this;
  }

  template<PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLevelIterator<1, pitype, GridImp> ::
  ~ALU2dGridLevelIterator()
  {}

  //! prefix increment
  template<PartitionIteratorType pitype, class GridImp>
  inline void ALU2dGridLevelIterator<1, pitype, GridImp> :: increment ()
  {
    // if already end iter, return
    if(endIter_) return ;

    IteratorType & iter = iter_;
    assert(myFace_>=0);

    int goNext = 1;
    item_ = &iter->getitem();
    int elIdx = item_->getIndex();

    while (myFace_ < 3) {
      int idx = item_->edge_idx(myFace_);
      // check if face is visited on this element
      if( marker().isOnElement(elIdx,idx,1) )
      {
        goNext = 0;
        break;
      }
      ++myFace_;
    }

    if (goNext)
    {
      assert(myFace_==3);

      iter->next();
      if(iter->done())
      {
        endIter_ = true;
        myFace_= 0;
        this->done();
        return ;
      }

      // get new element
      item_ = &iter->getitem();
#if ALU2DGRID_PARALLEL
      while ( ! this->grid_.rankManager().isValid( item_->getIndex(), pitype ) )
      {
        iter->next();
        if(iter->done())
        {
          endIter_ = true;
          myFace_= 0;
          this->done();
          return ;
        }
        item_ = &iter->getitem();
      }
#endif

      myFace_= 0;
      // whatever this does
      this->updateEntityPointer(item_, myFace_, level_);
      increment();
      return;
    }

    if(iter->done())
    {
      endIter_ = true;
      myFace_= 0; // set face to non valid value
      this->done();
      return ;
    }
    item_ = &iter->getitem();
    this->updateEntityPointer(item_, myFace_, level_);
    ++myFace_;
  }


  //********************************************************************
  //
  //  --ALU2dLevelLeafIterator
  //  --LevelIterator, specialized for cd=2
  //
  //********************************************************************

  //! constructor for cd=2
  template<PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLevelIterator<2, pitype, GridImp> ::
  ALU2dGridLevelIterator(const GridImp & grid, int level, bool end) :
    EntityPointerType (grid),
    endIter_(end),
    level_(level),
    myFace_(0),
    iter_(),
    marker_(& grid.getMarkerVector(level))
  {
    if(!end)
    {
      // update marker Vector if necessary
      if( ! marker().up2Date() ) marker().update(grid,level_);

      iter_ = IteratorType(grid.myGrid(), level_);
      iter_->first();

      if((!iter_->done()))
      {
        item_ = &iter_->getitem();
#if ALU2DGRID_PARALLEL
        while ( ! this->grid_.rankManager().isValid( item_->getIndex(), pitype ) )
        {
          iter_->next();
          if(iter_->done())
          {
            endIter_ = true;
            myFace_= 0;
            this->done();
            return ;
          }
          item_ = &iter_->getitem();
        }
#endif
        vertex_ = item_->getVertex(myFace_);
        this->updateEntityPointer(vertex_, myFace_, level_);
        increment();
      }
    }
    else
    {
      endIter_ = true;
      this->done();
    }
  }

  //! copy constructor for cd=2
  template<PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLevelIterator<2, pitype, GridImp> ::
  ALU2dGridLevelIterator(const ALU2dGridLevelIterator<2,pitype,GridImp> & org)
    : EntityPointerType (org)
      , endIter_( org.endIter_ )
      , level_( org.level_ )
      , myFace_(org.myFace_)
      , item_(org.item_)
      , vertex_(org.vertex_)
      , iter_ ( org.iter_ )
      , marker_(org.marker_)
  {}

  //! assignment
  template<PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLevelIterator<2, pitype, GridImp> &
  ALU2dGridLevelIterator<2, pitype, GridImp> ::
  operator = (const ThisType & org)
  {
    EntityPointerType :: operator = (org);
    endIter_ = org.endIter_ ;
    level_   = org.level_;
    myFace_  = org.myFace_;
    item_    = org.item_;
    vertex_  = org.vertex_;
    iter_    = org.iter_;
    marker_  = org.marker_;

    assert(marker_ == org.marker_);
    return *this;
  }


  template<PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLevelIterator<2, pitype, GridImp> ::
  ~ALU2dGridLevelIterator()
  {}

  //! prefix increment
  template<PartitionIteratorType pitype, class GridImp>
  inline void ALU2dGridLevelIterator<2, pitype, GridImp> :: increment () {

    if(endIter_)
      return ;

    IteratorType & iter = iter_;
    assert(myFace_>=0);

    int goNext = 1;
    item_ = &iter->getitem();
    int elIdx = item_->getIndex();

    while (myFace_ < 3) {
      vertex_ = item_->getVertex(myFace_);
      int idx = vertex_->getIndex();

      // check if face is visited on this element
      if( marker().isOnElement(elIdx,idx,2) )
      {
        goNext = 0;
        break;
      }
      ++myFace_;
    }

    if (goNext) {
      assert(myFace_==3);

      iter->next();
      if(iter->done())
      {
        endIter_ = true;
        myFace_= 0;
        this->done();
        return ;
      }
      item_ = &iter->getitem();

#if ALU2DGRID_PARALLEL
      while ( ! this->grid_.rankManager().isValid( item_->getIndex(), pitype ) )
      {
        iter->next();
        if(iter->done())
        {
          endIter_ = true;
          myFace_= 0;
          this->done();
          return ;
        }
        item_ = &iter->getitem();
      }
#endif

      myFace_ = 0;
      vertex_ = item_->getVertex(myFace_);
      this->updateEntityPointer(vertex_, myFace_, level_);
      increment();
      return;
    }

    if(iter->done()) {
      endIter_ = true;
      myFace_= 0;
      this->done();
      return ;
    }
    item_ = &iter->getitem();
    vertex_ = item_->getVertex(myFace_);
    this->updateEntityPointer(vertex_, myFace_, level_);
    ++myFace_;
  }


  //********************************************************************
  //
  //  --ALU2dGridHierarchicIterator
  //  --HierarchicIterator
  //
  //********************************************************************


  //! the normal Constructor
  template<class GridImp>
  inline ALU2dGridHierarchicIterator<GridImp> ::
  ALU2dGridHierarchicIterator(const GridImp &grid, const HElementType & elem, int maxlevel, bool end) :
    ALU2dGridEntityPointer<0,GridImp> (grid)
    , elem_(&elem)
    , maxlevel_(maxlevel)
    , endIter_(end) {

    if (!end)
    {
      HElementType * item = const_cast<HElementType *> (elem_->down());
      if(item)
      {
        // we have children and they lie in the disired level range
        if(item->level() <= maxlevel_)
        {
          this->updateEntityPointer( item );
        }
        else
        { // otherwise do nothing
          this->done();
        }
      }
      else
      {
        this->done();
      }
    }
  }


  //! the normal Constructor
  template<class GridImp>
  inline ALU2dGridHierarchicIterator<GridImp> ::
  ALU2dGridHierarchicIterator(const ALU2dGridHierarchicIterator<GridImp> &org) :
    ALU2dGridEntityPointer<0,GridImp> (org)
    , elem_ (org.elem_)
    , maxlevel_(org.maxlevel_)
    , endIter_(org.endIter_)  {}


  template <class GridImp>
  inline typename ALU2dGridHierarchicIterator<GridImp>::HElementType * ALU2dGridHierarchicIterator<GridImp>::
  goNextElement(HElementType * oldelem )
  {
    // strategy is:
    // - go down as far as possible and then over all children
    // - then go to father and next and down again

    HElementType * nextelem = oldelem->down();
    if(nextelem)
    {
      if(nextelem->level() <= maxlevel_)
        return nextelem;
    }

    nextelem = oldelem->next();
    if(nextelem)
    {
      if(nextelem->level() <= maxlevel_)
        return nextelem;
    }

    nextelem = oldelem->father();
    if(nextelem == elem_) return 0;

    while( !nextelem->next() )
    {
      nextelem = nextelem->father();
      if(nextelem == elem_) return 0;
    }

    if(nextelem) nextelem = nextelem->next();

    return nextelem;
  }


  //! increment, go to next entity
  template<class GridImp>
  inline void ALU2dGridHierarchicIterator<GridImp> :: increment() {

    assert(this->item_ != 0);

    HElementType * nextItem = goNextElement( this->item_ );
    if(!nextItem)
    {
      this->done();
      return ;
    }

    this->updateEntityPointer(nextItem);
    return ;
  }

} //end namespace Dune
#endif
