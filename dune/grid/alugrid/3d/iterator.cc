// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_ITERATOR_CC
#define DUNE_ALUGRID_ITERATOR_CC

#if COMPILE_ALUGRID_INLINE == 0
#include <config.h>
#endif

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
  ALU3dGridLevelIterator(const FactoryType& factory, int level, bool )
    : ALU3dGridEntityPointer<codim,GridImp> (factory,level)
      , level_(level)
      , iter_ (0)
  {
    const GridImp& grid = factory.grid();
    iter_  = new IteratorType ( grid, level_, grid.nlinks() );
    assert( iter_ );
    this->firstItem( grid, *this, level_);
  }

  // Constructor for end iterator
  template<int codim, PartitionIteratorType pitype, class GridImp >
  alu_inline ALU3dGridLevelIterator<codim,pitype,GridImp> ::
  ALU3dGridLevelIterator(const FactoryType& factory, int level)
    : ALU3dGridEntityPointer<codim,GridImp> (factory ,level)
      , level_(level)
      , iter_ (0)
  {
    this->done();
  }

  template<int codim, PartitionIteratorType pitype, class GridImp >
  alu_inline ALU3dGridLevelIterator<codim,pitype,GridImp> ::
  ALU3dGridLevelIterator(const ALU3dGridLevelIterator<codim,pitype,GridImp> & org )
    : ALU3dGridEntityPointer<codim,GridImp> ( org.factory_ , org.level_ )
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
        this->setItem( this->grid(), *this, *iter_, level_ );
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
    this->incrementIterator(this->grid(),*this,level_);
    return ;
  }

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  alu_inline typename ALU3dGridLevelIterator<cdim, pitype, GridImp> :: Entity &
  ALU3dGridLevelIterator<cdim, pitype, GridImp> :: dereference () const
  {
#ifndef NDEBUG
    const ALU3dGridLevelIterator<cdim, pitype, GridImp> endIterator ( this->factory_,level_);
    // assert that iterator not equals end iterator
    assert( ! this->equals(endIterator) );
#endif

    // don't dereference empty entity pointer
    assert( this->seed_.item() );
    assert( this->entity_ );
    assert( this->seed_.item() == & this->entityImp().getItem() );
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
  ALU3dGridLeafIterator( const FactoryType& factory, int level )
    : ALU3dGridEntityPointer <cdim,GridImp> ( factory, level )
      , iter_ (0)
      , walkLevel_(level)
  {
    this->done();
  }

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  alu_inline ALU3dGridLeafIterator<cdim, pitype, GridImp> ::
  ALU3dGridLeafIterator(const FactoryType& factory, int level ,
                        bool isBegin)
    : ALU3dGridEntityPointer <cdim,GridImp> ( factory, level )
      , iter_ (0)
      , walkLevel_(level)
  {
    const GridImp& grid = factory.grid();
    // create interior iterator
    iter_ = new IteratorType ( grid , level , grid.nlinks() );
    assert( iter_ );
    // -1 to identify as leaf iterator
    this->firstItem(grid,*this,-1);
  }

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  alu_inline ALU3dGridLeafIterator<cdim, pitype, GridImp> ::
  ALU3dGridLeafIterator(const ThisType & org)
    : ALU3dGridEntityPointer <cdim,GridImp> ( org.factory_, org.walkLevel_ )
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
        this->setItem( this->grid(),*this, *iter_,-1);
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
    this->incrementIterator(this->grid(), *this,-1);
    return ;
  }

  template<int cdim, PartitionIteratorType pitype, class GridImp>
  alu_inline typename ALU3dGridLeafIterator<cdim, pitype, GridImp> :: Entity &
  ALU3dGridLeafIterator<cdim, pitype, GridImp> :: dereference () const
  {
#ifndef NDEBUG
    const ALU3dGridLeafIterator<cdim, pitype, GridImp> endIterator (this->factory_, this->grid().maxLevel());
    // assert that iterator not equals end iterator
    assert( ! this->equals(endIterator) );
#endif

    // don't dereference empty entity pointer
    assert( this->seed_.item() );
    assert( this->entity_ );
    assert( this->seed_.item() == & this->entityImp().getItem() );
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
  ALU3dGridHierarchicIterator(const FactoryType& factory,
                              const HElementType & elem, int maxlevel ,bool end)
    : ALU3dGridEntityPointer<0,GridImp> ( factory, maxlevel )
      , elem_(&elem)
      , ghostElem_( )
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

  template <class GridImp>
  alu_inline ALU3dGridHierarchicIterator<GridImp> ::
  ALU3dGridHierarchicIterator(const FactoryType& factory,
                              const HBndSegType& ghost,
                              int maxlevel,
                              bool end)
    : ALU3dGridEntityPointer<0,GridImp> ( factory, maxlevel )
      , elem_( 0 )
      , ghostElem_( ghost )
      , maxlevel_(maxlevel)
  {
    if( ! end )
    {
      // lock entity pointer
      this->locked_ = true ;
      ghostElem_ = const_cast<HBndSegType *> (ghost.down());

      // we have children and they lie in the disired level range
      if( ghostElem_ != 0 && ghostElem_->ghostLevel() <= maxlevel_)
      {
        this->updateGhostPointer( *ghostElem_ );
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

  template <class GridImp>
  alu_inline ALU3dGridHierarchicIterator<GridImp> ::
  ALU3dGridHierarchicIterator(const ThisType& org)
    : ALU3dGridEntityPointer<0,GridImp> ( org.factory_, org.maxlevel_ )
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
    ghostElem_ = org.ghostElem_;
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
    assert(this->seed_.item() != 0);

    if( ghostElem_.valid() )
    {
      ghostElem_ = goNextElement( ghostElem_.ghost(), ghostElem_.nextGhost() );
      if( ! ghostElem_ )
      {
        this->done();
        return ;
      }

      this->updateGhostPointer( *ghostElem_ );
    }
    else
    {
      HElementType * nextItem = goNextElement( elem_, this->seed_.item() );
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
    assert( this->seed_.item() );
    assert( this->entity_ );
    assert( this->seed_.item() == & this->entityImp().getItem() );
    return (*this->entity_);
  }

#if COMPILE_ALUGRID_LIB
  // Instantiation
  template class ALU3dGrid< hexa, No_Comm >;
  template class ALU3dGrid< tetra, No_Comm >;

  // Instantiation without MPI
  template class ALU3dGridLeafIterator< 0, All_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 0, All_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLeafIterator< 0, Interior_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 0, Interior_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLeafIterator< 0, InteriorBorder_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 0, InteriorBorder_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLeafIterator< 0, Overlap_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 0, Overlap_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLeafIterator< 0, OverlapFront_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 0, OverlapFront_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLeafIterator< 0, Ghost_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 0, Ghost_Partition, const ALU3dGrid< tetra, No_Comm > >;

  template class ALU3dGridLeafIterator< 1, All_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 1, All_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLeafIterator< 1, Interior_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 1, Interior_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLeafIterator< 1, InteriorBorder_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 1, InteriorBorder_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLeafIterator< 1, Overlap_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 1, Overlap_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLeafIterator< 1, OverlapFront_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 1, OverlapFront_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLeafIterator< 1, Ghost_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 1, Ghost_Partition, const ALU3dGrid< tetra, No_Comm > >;

  template class ALU3dGridLeafIterator< 2, All_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 2, All_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLeafIterator< 2, Interior_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 2, Interior_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLeafIterator< 2, InteriorBorder_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 2, InteriorBorder_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLeafIterator< 2, Overlap_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 2, Overlap_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLeafIterator< 2, OverlapFront_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 2, OverlapFront_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLeafIterator< 2, Ghost_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 2, Ghost_Partition, const ALU3dGrid< tetra, No_Comm > >;

  template class ALU3dGridLeafIterator< 3, All_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 3, All_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLeafIterator< 3, Interior_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 3, Interior_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLeafIterator< 3, InteriorBorder_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 3, InteriorBorder_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLeafIterator< 3, Overlap_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 3, Overlap_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLeafIterator< 3, OverlapFront_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 3, OverlapFront_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLeafIterator< 3, Ghost_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLeafIterator< 3, Ghost_Partition, const ALU3dGrid< tetra, No_Comm > >;

  template class ALU3dGridLevelIterator< 0, All_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 0, All_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLevelIterator< 0, Interior_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 0, Interior_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLevelIterator< 0, InteriorBorder_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 0, InteriorBorder_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLevelIterator< 0, Overlap_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 0, Overlap_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLevelIterator< 0, OverlapFront_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 0, OverlapFront_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLevelIterator< 0, Ghost_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 0, Ghost_Partition, const ALU3dGrid< tetra, No_Comm > >;

  template class ALU3dGridLevelIterator< 1, All_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 1, All_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLevelIterator< 1, Interior_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 1, Interior_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLevelIterator< 1, InteriorBorder_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 1, InteriorBorder_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLevelIterator< 1, Overlap_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 1, Overlap_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLevelIterator< 1, OverlapFront_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 1, OverlapFront_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLevelIterator< 1, Ghost_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 1, Ghost_Partition, const ALU3dGrid< tetra, No_Comm > >;

  template class ALU3dGridLevelIterator< 2, All_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 2, All_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLevelIterator< 2, Interior_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 2, Interior_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLevelIterator< 2, InteriorBorder_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 2, InteriorBorder_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLevelIterator< 2, Overlap_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 2, Overlap_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLevelIterator< 2, OverlapFront_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 2, OverlapFront_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLevelIterator< 2, Ghost_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 2, Ghost_Partition, const ALU3dGrid< tetra, No_Comm > >;

  template class ALU3dGridLevelIterator< 3, All_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 3, All_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLevelIterator< 3, Interior_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 3, Interior_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLevelIterator< 3, InteriorBorder_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 3, InteriorBorder_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLevelIterator< 3, Overlap_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 3, Overlap_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLevelIterator< 3, OverlapFront_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 3, OverlapFront_Partition, const ALU3dGrid< tetra, No_Comm > >;
  template class ALU3dGridLevelIterator< 3, Ghost_Partition, const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridLevelIterator< 3, Ghost_Partition, const ALU3dGrid< tetra, No_Comm > >;

  template class ALU3dGridHierarchicIterator< const ALU3dGrid< hexa, No_Comm > >;
  template class ALU3dGridHierarchicIterator< const ALU3dGrid< tetra, No_Comm > >;

#if ALU3DGRID_PARALLEL
  // Instantiation
  template class ALU3dGrid< hexa, MPI_Comm >;
  template class ALU3dGrid< tetra, MPI_Comm >;

  // Instantiation with MPI
  template class ALU3dGridLeafIterator< 0, All_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 0, All_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 0, Interior_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 0, Interior_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 0, InteriorBorder_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 0, InteriorBorder_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 0, Overlap_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 0, Overlap_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 0, OverlapFront_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 0, OverlapFront_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 0, Ghost_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 0, Ghost_Partition, const ALU3dGrid< tetra, MPI_Comm > >;

  template class ALU3dGridLeafIterator< 1, All_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 1, All_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 1, Interior_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 1, Interior_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 1, InteriorBorder_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 1, InteriorBorder_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 1, Overlap_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 1, Overlap_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 1, OverlapFront_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 1, OverlapFront_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 1, Ghost_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 1, Ghost_Partition, const ALU3dGrid< tetra, MPI_Comm > >;

  template class ALU3dGridLeafIterator< 2, All_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 2, All_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 2, Interior_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 2, Interior_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 2, InteriorBorder_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 2, InteriorBorder_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 2, Overlap_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 2, Overlap_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 2, OverlapFront_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 2, OverlapFront_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 2, Ghost_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 2, Ghost_Partition, const ALU3dGrid< tetra, MPI_Comm > >;

  template class ALU3dGridLeafIterator< 3, All_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 3, All_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 3, Interior_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 3, Interior_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 3, InteriorBorder_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 3, InteriorBorder_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 3, Overlap_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 3, Overlap_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 3, OverlapFront_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 3, OverlapFront_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 3, Ghost_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLeafIterator< 3, Ghost_Partition, const ALU3dGrid< tetra, MPI_Comm > >;

  template class ALU3dGridLevelIterator< 0, All_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 0, All_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 0, Interior_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 0, Interior_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 0, InteriorBorder_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 0, InteriorBorder_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 0, Overlap_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 0, Overlap_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 0, OverlapFront_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 0, OverlapFront_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 0, Ghost_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 0, Ghost_Partition, const ALU3dGrid< tetra, MPI_Comm > >;

  template class ALU3dGridLevelIterator< 1, All_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 1, All_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 1, Interior_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 1, Interior_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 1, InteriorBorder_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 1, InteriorBorder_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 1, Overlap_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 1, Overlap_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 1, OverlapFront_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 1, OverlapFront_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 1, Ghost_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 1, Ghost_Partition, const ALU3dGrid< tetra, MPI_Comm > >;

  template class ALU3dGridLevelIterator< 2, All_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 2, All_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 2, Interior_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 2, Interior_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 2, InteriorBorder_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 2, InteriorBorder_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 2, Overlap_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 2, Overlap_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 2, OverlapFront_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 2, OverlapFront_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 2, Ghost_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 2, Ghost_Partition, const ALU3dGrid< tetra, MPI_Comm > >;

  template class ALU3dGridLevelIterator< 3, All_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 3, All_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 3, Interior_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 3, Interior_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 3, InteriorBorder_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 3, InteriorBorder_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 3, Overlap_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 3, Overlap_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 3, OverlapFront_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 3, OverlapFront_Partition, const ALU3dGrid< tetra, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 3, Ghost_Partition, const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridLevelIterator< 3, Ghost_Partition, const ALU3dGrid< tetra, MPI_Comm > >;

  template class ALU3dGridHierarchicIterator< const ALU3dGrid< hexa, MPI_Comm > >;
  template class ALU3dGridHierarchicIterator< const ALU3dGrid< tetra, MPI_Comm > >;
#endif // end ALU3DGRID_PARALLEL

#endif // end COMPILE_ALUGRID_LIB

} // end namespace Dune
#endif // DUNE_ALUGRID_ITERATOR_IMP_CC
