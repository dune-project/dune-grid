// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU3DGRIDMEMORY_HH
#define DUNE_ALU3DGRIDMEMORY_HH

#include <cassert>
#include <cstdlib>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <dune/common/finitestack.hh>

namespace Dune {

  //! organize the memory management for entitys used by the NeighborIterator
  template <class Object>
  class ALUMemoryProvider
  {
    enum { maxStackObjects = 256 };
    typedef FiniteStack< Object* , maxStackObjects > StackType ;
#ifdef _OPENMP
    std::vector< StackType  > objStackVec_;
#else
    StackType objStack_ ;
#endif

    typedef ALUMemoryProvider < Object > MyType;

    StackType& objStack()
    {
#ifdef _OPENMP
      assert( (int) objStackVec_.size() > omp_get_thread_num() );
      return objStackVec_[ omp_get_thread_num() ];
#else
      return objStack_ ;
#endif
    }

  public:
    typedef Object ObjectType;

    //! delete all objects stored in stack
    ALUMemoryProvider()
#ifdef _OPENMP
      : objStackVec_( omp_get_max_threads() )
#endif
    {}

    //! call deleteEntity
    ~ALUMemoryProvider ();

    //! i.e. return pointer to Entity
    template <class GridType>
    ObjectType * getObject(const GridType &grid, int level);

    //! i.e. return pointer to Entity
    template <class GridType, class EntityImp>
    inline ObjectType * getEntityObject(const GridType &grid, int level , EntityImp * fakePtr )
    {
      if( objStack().empty() )
      {
        return ( new ObjectType(EntityImp(grid,level) ));
      }
      else
      {
        return stackObject();
      }
    }

    //! i.e. return pointer to Entity
    ObjectType * getObjectCopy(const ObjectType & org);

    //! free, move element to stack, returns NULL
    void freeObject (ObjectType * obj);

  protected:
    inline ObjectType * stackObject()
    {
      assert( ! objStack().empty() );
      // finite stack does also return object on pop
      return objStack().pop();
    }

  private:
    //! do not copy pointers
    ALUMemoryProvider(const ALUMemoryProvider<Object> & org);
  };


  //************************************************************************
  //
  //  ALUMemoryProvider implementation
  //
  //************************************************************************
  template <class Object> template <class GridType>
  inline typename ALUMemoryProvider<Object>::ObjectType *
  ALUMemoryProvider<Object>::getObject
    (const GridType &grid, int level )
  {
    if( objStack().empty() )
    {
      return ( new Object (grid,level) );
    }
    else
    {
      return stackObject();
    }
  }

  template <class Object>
  inline typename ALUMemoryProvider<Object>::ObjectType *
  ALUMemoryProvider<Object>::getObjectCopy
    (const ObjectType & org )
  {
    if( objStack().empty() )
    {
      return ( new Object (org) );
    }
    else
    {
      return stackObject();
    }
  }

  template <class Object>
  inline ALUMemoryProvider<Object>::~ALUMemoryProvider()
  {
#ifdef _OPENMP
    for(size_t i = 0; i< objStackVec_.size(); ++i)
    {
      StackType& objStk = objStackVec_[ i ];
#else
    {
      StackType& objStk = objStack_;
#endif
      while ( ! objStk.empty() )
      {
        ObjectType * obj = objStk.pop();
        delete obj;
      }
    }
  }

  template <class Object>
  inline void ALUMemoryProvider<Object>::freeObject(Object * obj)
  {
    StackType& stk = objStack();
    if( stk.full() )
      delete obj;
    else
      stk.push( obj );
  }

#undef USE_FINITE_STACK

} // end namespace Dune

#endif
