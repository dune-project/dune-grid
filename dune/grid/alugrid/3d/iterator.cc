// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_ITERATOR_CC
#define DUNE_ALUGRID_ITERATOR_CC

#include "alu3dinclude.hh"
#include "iterator.hh"

namespace Dune {

  /*************************************************************************
  #       ######  #    #  ######  #          #     #####  ######  #####
  #       #       #    #  #       #          #       #    #       #    #
  #       #####   #    #  #####   #          #       #    #####   #    #
  #       #       #    #  #       #          #       #    #       #####
  #       #        #  #   #       #          #       #    #       #   #
  ######  ######    ##    ######  ######     #       #    ######  #    #
  *************************************************************************/
  //--LevelIterator
  // Constructor for begin iterator
  template<int codim, PartitionIteratorType pitype, class GridImp >
  alu_inline ALU3dGridLevelIterator<codim,pitype,GridImp> ::
  ALU3dGridLevelIterator(const GridImp & grid, int level, bool )
    : ALU3dGridEntityPointer<codim,GridImp> (grid,level)
      , level_(level)
      , iter_ (0)
  {
    iter_  = new IteratorType ( this->grid_ , level_, grid.nlinks() );
    assert( iter_ );
    this->firstItem(this->grid_,*this,level_);
  }

  // Constructor for end iterator
  template<int codim, PartitionIteratorType pitype, class GridImp >
  alu_inline ALU3dGridLevelIterator<codim,pitype,GridImp> ::
  ALU3dGridLevelIterator(const GridImp & grid, int level)
    : ALU3dGridEntityPointer<codim,GridImp> (grid ,level)
      , level_(level)
      , iter_ (0)
  {
    this->done();
  }

  template<int codim, PartitionIteratorType pitype, class GridImp >
  alu_inline ALU3dGridLevelIterator<codim,pitype,GridImp> ::
  ALU3dGridLevelIterator(const ALU3dGridLevelIterator<codim,pitype,GridImp> & org )
    : ALU3dGridEntityPointer<codim,GridImp> ( org.grid_ , org.level_ )
      , level_( org.level_ )
      , iter_(0)
  {
    assign(org);
  }

  template<int codim, PartitionIteratorType pitype, class GridImp>
  alu_inline ALU3dGridLevelIterator<codim, pitype, GridImp> ::
  ~ALU3dGridLevelIterator ()
  {
    removeIter();
  }

  template<int codim, PartitionIteratorType pitype, class GridImp>
  alu_inline void ALU3dGridLevelIterator<codim, pitype, GridImp> ::
  removeIter ()
  {
    this->done();
    if(iter_)
    {
      delete iter_;
      iter_ = 0;
    }
  }

  template<int codim, PartitionIteratorType pitype, class GridImp>
  alu_inline void ALU3dGridLevelIterator<codim, pitype, GridImp> ::
  assign(const ThisType & org)
  {
    assert( iter_ == 0 );
    ALU3dGridEntityPointer <codim,GridImp> :: clone (org);
    level_ = org.level_;
    if( org.iter_ )
    {
      iter_ = new IteratorType ( *(org.iter_) );
      assert( iter_ );
      if(!(iter_->done()))
      {
        this->setItem(this->grid_, *this, *iter_, level_ );
        assert( this->equals(org) );
      }
    }
    else
    {
      this->done();
    }
  }
  template<int codim, PartitionIteratorType pitype, class GridImp>
  alu_inline ALU3dGridLevelIterator<codim, pitype, GridImp> &
  ALU3dGridLevelIterator<codim, pitype, GridImp> ::
  operator = (const ThisType & org)
  {
    removeIter();
    assign(org);
    return *this;
  }

  template<int codim, PartitionIteratorType pitype, class GridImp >
  alu_inline void ALU3dGridLevelIterator<codim,pitype,GridImp> :: increment ()
  {
    this->incrementIterator(this->grid_,*this,level_);
    return ;
  }

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  alu_inline typename ALU3dGridLevelIterator<cdim, pitype, GridImp> :: Entity &
  ALU3dGridLevelIterator<cdim, pitype, GridImp> :: dereference () const
  {
#ifndef NDEBUG
    const ALU3dGridLevelIterator<cdim, pitype, GridImp> endIterator (this->grid_,level_);
    // assert that iterator not equals end iterator
    assert( ! this->equals(endIterator) );
#endif

    // don't dereference empty entity pointer
    assert( this->item_ );
    assert( this->entity_ );
    assert( this->item_ == & this->entityImp().getItem() );
    return (*this->entity_);
  }

  //*******************************************************************
  //
  //  LEAFITERATOR
  //
  //--LeafIterator
  //*******************************************************************
  // constructor for end iterators
  template<int cdim, PartitionIteratorType pitype, class GridImp>
  alu_inline ALU3dGridLeafIterator<cdim, pitype, GridImp> ::
  ALU3dGridLeafIterator( const GridImp &grid, int level )
    : ALU3dGridEntityPointer <cdim,GridImp> ( grid, level )
      , iter_ (0)
      , walkLevel_(level)
  {
    this->done();
  }

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  alu_inline ALU3dGridLeafIterator<cdim, pitype, GridImp> ::
  ALU3dGridLeafIterator(const GridImp &grid, int level ,
                        bool isBegin)
    : ALU3dGridEntityPointer <cdim,GridImp> ( grid, level )
      , iter_ (0)
      , walkLevel_(level)
  {
    // create interior iterator
    iter_ = new IteratorType ( this->grid_ , level , grid.nlinks() );
    assert( iter_ );
    // -1 to identify as leaf iterator
    this->firstItem(this->grid_,*this,-1);
  }

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  alu_inline ALU3dGridLeafIterator<cdim, pitype, GridImp> ::
  ALU3dGridLeafIterator(const ThisType & org)
    : ALU3dGridEntityPointer <cdim,GridImp> ( org.grid_, org.walkLevel_ )
      , iter_(0)
      , walkLevel_(org.walkLevel_)
  {
    // assign iterator without cloning entity pointer again
    assign(org);
  }

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  alu_inline ALU3dGridLeafIterator<cdim, pitype, GridImp> ::
  ~ALU3dGridLeafIterator()
  {
    removeIter();
  }

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  alu_inline void ALU3dGridLeafIterator<cdim, pitype, GridImp> ::
  removeIter ()
  {
    this->done();
    if(iter_)
    {
      delete iter_;
      iter_ = 0;
    }
  }

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  alu_inline ALU3dGridLeafIterator<cdim, pitype, GridImp> &
  ALU3dGridLeafIterator<cdim, pitype, GridImp> ::
  operator = (const ThisType & org)
  {
    removeIter();
    assign(org);
    return *this;
  }

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  alu_inline void ALU3dGridLeafIterator<cdim, pitype, GridImp> ::
  assign (const ThisType & org)
  {
    assert( iter_ == 0 );
    ALU3dGridEntityPointer <cdim,GridImp> :: clone (org);

    if( org.iter_ )
    {
      assert( !org.iter_->done() );
      iter_ = new IteratorType ( *(org.iter_) );
      assert( iter_ );

      if( !(iter_->done() ))
      {
        assert( !iter_->done());
        assert( !org.iter_->done() );
        // -1 to identify leaf iterator
        this->setItem(this->grid_,*this, *iter_,-1);
        assert( this->equals(org) );
      }
    }
    else
    {
      this->done();
    }

    walkLevel_ = org.walkLevel_;
  }

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  alu_inline void ALU3dGridLeafIterator<cdim, pitype, GridImp> :: increment ()
  {
    // -1 to identify leaf iterator
    this->incrementIterator(this->grid_,*this,-1);
    return ;
  }

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  alu_inline typename ALU3dGridLeafIterator<cdim, pitype, GridImp> :: Entity &
  ALU3dGridLeafIterator<cdim, pitype, GridImp> :: dereference () const
  {
#ifndef NDEBUG
    const ALU3dGridLeafIterator<cdim, pitype, GridImp> endIterator (this->grid_, this->grid_.maxLevel());
    // assert that iterator not equals end iterator
    assert( ! this->equals(endIterator) );
#endif

    // don't dereference empty entity pointer
    assert( this->item_ );
    assert( this->entity_ );
    assert( this->item_ == & this->entityImp().getItem() );
    return (*this->entity_);
  }


  /************************************************************************************
  #     #
  #     #     #    ######  #####      #     #####  ######  #####
  #     #     #    #       #    #     #       #    #       #    #
  #######     #    #####   #    #     #       #    #####   #    #
  #     #     #    #       #####      #       #    #       #####
  #     #     #    #       #   #      #       #    #       #   #
  #     #     #    ######  #    #     #       #    ######  #    #
  ************************************************************************************/
  // --HierarchicIterator
  template <class GridImp>
  alu_inline ALU3dGridHierarchicIterator<GridImp> ::
  ALU3dGridHierarchicIterator(const GridImp & grid ,
                              const HElementType & elem, int maxlevel ,bool end)
    : ALU3dGridEntityPointer<0,GridImp> ( grid, maxlevel )
      , elem_(&elem)
      , ghost_( 0 )
      , nextGhost_( 0 )
      , maxlevel_(maxlevel)
  {
    if (!end)
    {
      HElementType * item =
        const_cast<HElementType *> (elem.down());
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

#ifdef ALU3DGRID_PARALLEL
  template <class GridImp>
  alu_inline ALU3dGridHierarchicIterator<GridImp> ::
  ALU3dGridHierarchicIterator(const GridImp & grid ,
                              const HBndSegType& ghost,
                              int maxlevel,
                              bool end)
    : ALU3dGridEntityPointer<0,GridImp> ( grid, maxlevel )
      , elem_( 0 )
      , ghost_( &ghost )
      , nextGhost_( 0 )
      , maxlevel_(maxlevel)
  {
    if( ! end )
    {
      // lock entity pointer
      this->locked_ = true ;
      nextGhost_ = const_cast<HBndSegType *> (ghost.down());

      // we have children and they lie in the disired level range
      if( nextGhost_ && nextGhost_->ghostLevel() <= maxlevel_)
      {
        this->updateGhostPointer( *nextGhost_ );
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
#endif

  template <class GridImp>
  alu_inline ALU3dGridHierarchicIterator<GridImp> ::
  ALU3dGridHierarchicIterator(const ThisType& org)
    : ALU3dGridEntityPointer<0,GridImp> ( org.grid_, org.maxlevel_ )
  {
    assign( org );
  }

  template <class GridImp>
  alu_inline ALU3dGridHierarchicIterator<GridImp> &
  ALU3dGridHierarchicIterator<GridImp> ::
  operator = (const ThisType& org)
  {
    assign( org );
    return *this;
  }

  template <class GridImp>
  alu_inline void
  ALU3dGridHierarchicIterator<GridImp> ::
  assign(const ThisType& org)
  {
    // copy my data
    elem_      = org.elem_;
#ifdef ALU3DGRID_PARALLEL
    ghost_     = org.ghost_;
    nextGhost_ = org.nextGhost_;
#endif
    maxlevel_  = org.maxlevel_;

    // copy entity pointer
    // this method will probably free entity
    ALU3dGridEntityPointer<0,GridImp> :: clone(org);
  }

  template <class GridImp>
  alu_inline int
  ALU3dGridHierarchicIterator<GridImp>::
  getLevel(const HBndSegType* face) const
  {
    // return ghost level
    assert( face );
    return face->ghostLevel();
  }

  template <class GridImp>
  alu_inline int
  ALU3dGridHierarchicIterator<GridImp>::
  getLevel(const HElementType * item) const
  {
    // return normal level
    assert( item );
    return item->level();
  }
  template <class GridImp>
  template <class HItemType>
  alu_inline HItemType*
  ALU3dGridHierarchicIterator<GridImp>::
  goNextElement(const HItemType* startElem, HItemType * oldelem )
  {
    // strategy is:
    // - go down as far as possible and then over all children
    // - then go to father and next and down again

    HItemType * nextelem = oldelem->down();
    if(nextelem)
    {
      // use getLevel method
      if( getLevel(nextelem) <= maxlevel_)
        return nextelem;
    }

    nextelem = oldelem->next();
    if(nextelem)
    {
      // use getLevel method
      if( getLevel(nextelem) <= maxlevel_)
        return nextelem;
    }

    nextelem = oldelem->up();
    if(nextelem == startElem) return 0;

    while( !nextelem->next() )
    {
      nextelem = nextelem->up();
      if(nextelem == startElem) return 0;
    }

    if(nextelem) nextelem = nextelem->next();

    return nextelem;
  }

  template <class GridImp>
  alu_inline void ALU3dGridHierarchicIterator<GridImp> :: increment ()
  {
    assert(this->item_ != 0);

#ifdef ALU3DGRID_PARALLEL
    if( ghost_ )
    {
      assert( nextGhost_ );
      nextGhost_ = goNextElement( ghost_, nextGhost_ );
      if( ! nextGhost_ )
      {
        this->done();
        return ;
      }

      this->updateGhostPointer( *nextGhost_ );
    }
    else
#endif
    {
      HElementType * nextItem = goNextElement( elem_, this->item_ );
      if( ! nextItem)
      {
        this->done();
        return ;
      }

      this->updateEntityPointer(nextItem);
    }
    return ;
  }

  template <class GridImp>
  alu_inline typename ALU3dGridHierarchicIterator<GridImp> :: Entity &
  ALU3dGridHierarchicIterator<GridImp> :: dereference () const
  {
    // don't dereference empty entity pointer
    assert( this->item_ );
    assert( this->entity_ );
    assert( this->item_ == & this->entityImp().getItem() );
    return (*this->entity_);
  }

#if COMPILE_ALUGRID_LIB
  template class ALU3dGridLeafIterator< 0, All_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLeafIterator< 0, All_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLeafIterator< 0, Interior_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLeafIterator< 0, Interior_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLeafIterator< 0, InteriorBorder_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLeafIterator< 0, InteriorBorder_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLeafIterator< 0, Ghost_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLeafIterator< 0, Ghost_Partition, const ALU3dGrid<3,3,tetra> >;

  template class ALU3dGridLeafIterator< 1, All_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLeafIterator< 1, All_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLeafIterator< 1, Interior_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLeafIterator< 1, Interior_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLeafIterator< 1, InteriorBorder_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLeafIterator< 1, InteriorBorder_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLeafIterator< 1, Ghost_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLeafIterator< 1, Ghost_Partition, const ALU3dGrid<3,3,tetra> >;

  template class ALU3dGridLeafIterator< 2, All_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLeafIterator< 2, All_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLeafIterator< 2, Interior_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLeafIterator< 2, Interior_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLeafIterator< 2, InteriorBorder_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLeafIterator< 2, InteriorBorder_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLeafIterator< 2, Ghost_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLeafIterator< 2, Ghost_Partition, const ALU3dGrid<3,3,tetra> >;

  template class ALU3dGridLeafIterator< 3, All_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLeafIterator< 3, All_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLeafIterator< 3, Interior_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLeafIterator< 3, Interior_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLeafIterator< 3, InteriorBorder_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLeafIterator< 3, InteriorBorder_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLeafIterator< 3, Ghost_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLeafIterator< 3, Ghost_Partition, const ALU3dGrid<3,3,tetra> >;

  template class ALU3dGridLevelIterator< 0, All_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLevelIterator< 0, All_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLevelIterator< 0, Interior_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLevelIterator< 0, Interior_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLevelIterator< 0, InteriorBorder_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLevelIterator< 0, InteriorBorder_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLevelIterator< 0, Ghost_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLevelIterator< 0, Ghost_Partition, const ALU3dGrid<3,3,tetra> >;

  template class ALU3dGridLevelIterator< 1, All_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLevelIterator< 1, All_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLevelIterator< 1, Interior_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLevelIterator< 1, Interior_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLevelIterator< 1, InteriorBorder_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLevelIterator< 1, InteriorBorder_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLevelIterator< 1, Ghost_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLevelIterator< 1, Ghost_Partition, const ALU3dGrid<3,3,tetra> >;

  template class ALU3dGridLevelIterator< 2, All_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLevelIterator< 2, All_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLevelIterator< 2, Interior_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLevelIterator< 2, Interior_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLevelIterator< 2, InteriorBorder_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLevelIterator< 2, InteriorBorder_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLevelIterator< 2, Ghost_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLevelIterator< 2, Ghost_Partition, const ALU3dGrid<3,3,tetra> >;

  template class ALU3dGridLevelIterator< 3, All_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLevelIterator< 3, All_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLevelIterator< 3, Interior_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLevelIterator< 3, Interior_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLevelIterator< 3, InteriorBorder_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLevelIterator< 3, InteriorBorder_Partition, const ALU3dGrid<3,3,tetra> >;
  template class ALU3dGridLevelIterator< 3, Ghost_Partition, const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridLevelIterator< 3, Ghost_Partition, const ALU3dGrid<3,3,tetra> >;

  template class ALU3dGridHierarchicIterator< const ALU3dGrid<3,3,hexa> >;
  template class ALU3dGridHierarchicIterator< const ALU3dGrid<3,3,tetra> >;
#endif

} // end namespace Dune
#endif // DUNE_ALUGRID_ITERATOR_IMP_CC
