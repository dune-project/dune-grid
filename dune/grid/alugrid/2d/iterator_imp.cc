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


namespace Dune
{

  //********************************************************************
  //
  //  --ALU2dGridLeafIterator
  //  --LeafIterator
  //
  //********************************************************************

  //! constructor
  template<int cdim, PartitionIteratorType pitype, class GridImp>
  inline ALU2dGridLeafIterator<cdim, pitype, GridImp> ::
  ALU2dGridLeafIterator(const FactoryType& factory, bool end) :
    EntityPointerType ( factory ),
    endIter_(end),
    level_(-1),
    elem_(0),
    iter_(),
    marker_( factory.grid().getLeafMarker())
  {
    if(!end)
    {
      const GridImp& grid = factory.grid() ;
      // update marker Vector
      if( (cdim == 2) && (! marker_.valid()) ) marker_.update(grid);

      iter_ = IteratorType(grid.myGrid());
      iter_->first();

      if((!iter_->done()))
      {
        elem_ = &(iter_->getitem());

#if ALU2DGRID_PARALLEL
        const bool valid = (cdim == 0) ?
                           (this->grid().rankManager().isValid( elem_->getIndex(), pitype )) :
                           marker_.isValidVertex( elem_->getIndex() );
#else
        const bool valid = true;
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
                 (this->grid().rankManager().isValid( elem_->getIndex(), pitype )) :
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
              (this->grid().rankManager().isValid( elem_->getIndex(), pitype )) :
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
  ALU2dGridLeafIterator(const FactoryType& factory, bool end) :
    EntityPointerType ( factory ),
    endIter_(end),
    level_(-1),
    face_(0),
    elem_(0),
    iter_(),
    marker_( factory.grid().getLeafMarker())
  {
    if(!end)
    {
      const GridImp& grid = factory.grid() ;
      // update marker Vector
      if( ! marker_.valid() ) marker_.update(grid);

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
      assert( (eltype == ALU2DSPACE triangle && face_==3) ||
              (eltype == ALU2DSPACE quadrilateral && face_==4) );

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
      while ( ! this->grid().rankManager().isValid( elem_->getIndex(), pitype ) )
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
  inline int ALU2dGridLeafIterator<1, pitype, GridImp> :: goNextElement()
  {
    assert(face_>=0);
    const ElementType* item = this->seed_.item();
    assert( item );
    int elIdx = item->getIndex();

    while (face_ < item->numfaces() ) {
      int idx = item->edge_idx(face_);
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
  ALU2dGridLevelIterator(const FactoryType& factory, int level, bool end) :
    EntityPointerType ( factory ),
    endIter_(end),
    level_(level),
    iter_()
  {
    if(!end)
    {
      const GridImp& grid = factory.grid() ;

      iter_ = IteratorType(grid.myGrid(), level_);
      iter_->first();
      if((!iter_->done()))
      {
        item_ = &(iter_->getitem());
        this->updateEntityPointer(item_, -1 , level_);
#if ALU2DGRID_PARALLEL
        if ( ! this->grid().rankManager().isValid( item_->getIndex(), pitype ) )
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
    while ( ! this->grid().rankManager().isValid( item_->getIndex(), pitype ) )
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
  ALU2dGridLevelIterator(const FactoryType& factory, int level, bool end) :
    EntityPointerType ( factory ),
    endIter_(end),
    level_(level),
    myFace_(0),
    iter_(),
    marker_(& factory.grid().getMarkerVector(level) )
  {
    if(!end)
    {
      const GridImp& grid = factory.grid() ;

      // update marker Vector if necessary
      if( ! marker().valid() ) marker().update( grid, level_);

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

    while (myFace_ < item_->numfaces() ) {
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
      assert( (eltype == ALU2DSPACE triangle && myFace_==3) ||
              (eltype == ALU2DSPACE quadrilateral && myFace_==4) );

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
      while ( ! this->grid().rankManager().isValid( item_->getIndex(), pitype ) )
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
  ALU2dGridLevelIterator(const FactoryType& factory, int level, bool end) :
    EntityPointerType ( factory ),
    endIter_(end),
    level_(level),
    myFace_(0),
    iter_(),
    marker_(& factory.grid().getMarkerVector(level))
  {
    if(!end)
    {
      const GridImp& grid = factory.grid() ;

      // update marker Vector if necessary
      if( ! marker().valid() ) marker().update(grid,level_);

      iter_ = IteratorType(grid.myGrid(), level_);
      iter_->first();

      if((!iter_->done()))
      {
        item_ = &iter_->getitem();
#if ALU2DGRID_PARALLEL
        while ( ! this->grid().rankManager().isValid( item_->getIndex(), pitype ) )
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

    while (myFace_ < item_->numfaces() ) {
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
      assert( (eltype == ALU2DSPACE triangle && myFace_==3) ||
              (eltype == ALU2DSPACE quadrilateral && myFace_==4) );

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
      while ( ! this->grid().rankManager().isValid( item_->getIndex(), pitype ) )
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
  ALU2dGridHierarchicIterator(const FactoryType& factory, const HElementType & elem, int maxlevel, bool end)
    : ALU2dGridEntityPointer<0,GridImp>( factory )
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

    assert( this->seed_.item() != 0);

    HElementType * nextItem = goNextElement( this->seed_.item() );
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
