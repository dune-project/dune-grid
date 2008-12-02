// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_HIERARCHICITERATOR_CC
#define DUNE_ALBERTA_HIERARCHICITERATOR_CC

#include <dune/grid/albertagrid/hierarchiciterator.hh>

namespace Dune
{

  template< class GridImp >
  inline void AlbertaGridHierarchicIterator<GridImp>::
  makeIterator()
  {
    virtualEntity_.setTraverseStack(0);
    virtualEntity_.setElInfo(0,0,0,0);
  }

  template< class GridImp >
  inline AlbertaGridHierarchicIterator<GridImp>::
  AlbertaGridHierarchicIterator(const GridImp & grid,
                                int actLevel,
                                int maxLevel)
    : AlbertaGridEntityPointer<0,GridImp> (grid,actLevel,true,true)
      , startLevel_(actLevel)
      , level_ (actLevel)
      , maxlevel_ (maxLevel)
      , virtualEntity_( this->entityImp() )
      , end_ (true)
  {
    makeIterator();
  }

  template< class GridImp >
  inline AlbertaGridHierarchicIterator<GridImp>::
  AlbertaGridHierarchicIterator(const GridImp & grid,
                                ALBERTA TRAVERSE_STACK *travStack,int actLevel, int maxLevel, bool leafIt )
    : AlbertaGridEntityPointer<0,GridImp> (grid,actLevel,leafIt,false)
      , startLevel_(actLevel)
      , level_ (actLevel)
      , maxlevel_ ( maxLevel)
      , virtualEntity_( this->entityImp() )
      , end_ (false)
  {
    if(travStack)
    {
      // get new ALBERTA TRAVERSE STACK
      manageStack_.create();

      ALBERTA TRAVERSE_STACK *stack = manageStack_.getStack();

      // cut old traverse stack, kepp only actual element
      ALBERTA copyTraverseStack(stack, travStack );

      // set new traverse level
      if(maxlevel_ < 0)
      {
        std::cout << "WARNING: maxlevel < 0 in AlbertaGridHierarchicIterator! \n";
        // this means, we go until leaf level
        stack->traverse_fill_flag = CALL_LEAF_EL | stack->traverse_fill_flag;
        // exact here has to stand Grid->maxlevel, but is ok anyway
        maxlevel_ = this->grid_.maxLevel();
      }

      // set new traverse level
      stack->traverse_level = maxlevel_;

      virtualEntity_.setTraverseStack(stack);

      // Hier kann ein beliebiges Geometry uebergeben werden,
      // da jedes AlbertGeometry einen Zeiger auf das Macroelement
      // enthaelt.
      ALBERTA EL_INFO * elInfo = firstChild(stack);

      virtualEntity_.setElInfo(elInfo);
    }
    else
    {
      std::cout << "Warning: travStack == NULL in HierarchicIterator(travStack,travLevel) \n";
      makeIterator();
    }
  }

  template< class GridImp >
  inline AlbertaGridHierarchicIterator<GridImp>::
  AlbertaGridHierarchicIterator(const AlbertaGridHierarchicIterator<GridImp> & org)
    : AlbertaGridEntityPointer<0,GridImp> (org.grid_,org.level(),true, org.end_ )
      , startLevel_( org.startLevel_ )
      , level_ ( org.level_ )
      , maxlevel_ ( org.maxlevel_ )
      , virtualEntity_( this->entityImp() )
  {
    if( org.virtualEntity_.getElInfo() )
    {
      // get new ALBERTA TRAVERSE STACK
      manageStack_.create();
      ALBERTA TRAVERSE_STACK *stack = manageStack_.getStack();
      // cut old traverse stack, kepp only actual element
      ALBERTA copyTraverseStack(stack, org.manageStack_.getStack() );

      virtualEntity_.setTraverseStack( stack );
      /// get the actual used enInfo
      ALBERTA EL_INFO * elInfo = stack->elinfo_stack+stack->stack_used;

      virtualEntity_.setElInfo( elInfo );
    }
    else
      this->done();
  }

  template< class GridImp >
  inline AlbertaGridHierarchicIterator<GridImp> &
  AlbertaGridHierarchicIterator<GridImp>::
  operator = (const AlbertaGridHierarchicIterator<GridImp> & org)
  {
    const_cast<int &> (startLevel_) = org.startLevel_;
    level_    = org.level_;
    maxlevel_ = org.maxlevel_;

    if(org.manageStack_.stackExists())
    {
      // full copy of stack
      manageStack_.create();
      ALBERTA TRAVERSE_STACK * stack = manageStack_.getStack();
      const ALBERTA TRAVERSE_STACK * orgStack = org.manageStack_.getStack();
      ALBERTA copyTraverseStack( stack , orgStack );
    }

    if( org.virtualEntity_.getElInfo() )
      virtualEntity_.setEntity( org.virtualEntity_ );
    else
      this->done();
    return *this;
  }

  template< class GridImp >
  inline void AlbertaGridHierarchicIterator< GridImp >::increment()
  {
    ALBERTA EL_INFO * nextinfo = recursiveTraverse(manageStack_.getStack());

    if(!nextinfo)
    {
      this->done();
      return;
    }

    virtualEntity_.setElInfo( nextinfo );
    return ;
  }

  template< class GridImp >
  inline ALBERTA EL_INFO *
  AlbertaGridHierarchicIterator<GridImp>::
  firstChild (ALBERTA TRAVERSE_STACK * stack)
  {
    assert(stack);
    assert(stack->elinfo_stack);

    // stack_used ist the actual element
    stack->stack_used = startLevel_+1;

    // info_stack is 0, we want to visit both children
    stack->info_stack[stack->stack_used] = 0;

    ALBERTA EL * el = stack->elinfo_stack[stack->stack_used].el;

    // go down next child
    if(el->child[0] && (stack->traverse_level >
                        (stack->elinfo_stack+stack->stack_used)->level) )
    {
      if(stack->stack_used >= stack->stack_size - 1)
        ALBERTA enlargeTraverseStack(stack);

      int i = stack->info_stack[stack->stack_used];
      el = el->child[i];
      stack->info_stack[stack->stack_used]++;

      // new: go down maxlevel, but fake the elements
      level_++;
      this->grid_.fillElInfo(i, level_, stack->elinfo_stack+stack->stack_used,
                             stack->elinfo_stack+stack->stack_used+1 ,true);

      stack->stack_used++;
      stack->info_stack[stack->stack_used] = 0;
      return (stack->elinfo_stack + stack->stack_used);
    }
    else
    {
      return 0;
    }
  }

  template< class GridImp >
  inline ALBERTA EL_INFO *
  AlbertaGridHierarchicIterator<GridImp>::
  recursiveTraverse(ALBERTA TRAVERSE_STACK * stack)
  {
    // see function
    // static EL_INFO *traverse_leaf_el(TRAVERSE_STACK *stack)
    // Common/traverse_nr_common.cc, line 392
    ALBERTA EL * el=0;

    if(!stack->elinfo_stack)
    {
      /* somethin' wrong */
      return 0;
    }
    else
    {
      // go up until we can go down again
      el = stack->elinfo_stack[stack->stack_used].el;

      // stack->stack_used is actual element in stack
      // stack->info_stack[stack->stack_used] >= 2
      //    means the two children has been visited
      while((stack->stack_used-startLevel_ > 0) &&
            ((stack->info_stack[stack->stack_used] >= 2)
             || (el->child[0] == 0)
             || ( stack->traverse_level <=
                  (stack->elinfo_stack+stack->stack_used)->level)) )
      {
        stack->stack_used--;
        el = stack->elinfo_stack[stack->stack_used].el;
        level_ = stack->elinfo_stack[stack->stack_used].level;
      }

      // goto next father is done by other iterator and not our problem
      if(stack->stack_used-startLevel_ < 1)
      {
        return 0;
      }
    }

    // go down next child
    if(el->child[0] && (stack->traverse_level >
                        (stack->elinfo_stack+stack->stack_used)->level) )
    {
      if(stack->stack_used >= stack->stack_size - 1)
        ALBERTA enlargeTraverseStack(stack);

      int i = stack->info_stack[stack->stack_used];
      el = el->child[i];
      stack->info_stack[stack->stack_used]++;

      // new: go down maxlevel, but fake the elements
      level_++;
      this->grid_.fillElInfo(i, level_, stack->elinfo_stack+stack->stack_used,
                             stack->elinfo_stack+stack->stack_used+1 ,true);

      stack->stack_used++;
      stack->info_stack[stack->stack_used] = 0;
    }
    else
    {
      return 0;
    }

    return (stack->elinfo_stack + stack->stack_used);
  }  // recursive traverse over all childs

}

#endif
