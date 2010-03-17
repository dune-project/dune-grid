// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU3DGRIDMEMORY_HH
#define DUNE_ALU3DGRIDMEMORY_HH

#include <cassert>
#include <cstdlib>

#define USE_FINITE_STACK

#ifdef USE_FINITE_STACK
#include <dune/common/finitestack.hh>
#else
#include <stack>
#endif

namespace Dune {

  //! organize the memory management for entitys used by the NeighborIterator
  template <class Object>
  class ALUMemoryProvider
  {
#ifdef USE_FINITE_STACK
    enum { maxStackObjects = 256 };
    FiniteStack< Object* , maxStackObjects > objStack_;
#else
    std::stack < Object* > objStack_;
#endif

    typedef ALUMemoryProvider < Object > MyType;

  public:
    typedef Object ObjectType;

    //! delete all objects stored in stack
    ALUMemoryProvider() {};

    //! call deleteEntity
    ~ALUMemoryProvider ();

    //! i.e. return pointer to Entity
    template <class GridType>
    ObjectType * getObject(const GridType &grid, int level);

    //! i.e. return pointer to Entity
    template <class GridType, class EntityImp>
    inline ObjectType * getEntityObject(const GridType &grid, int level , EntityImp * fakePtr )
    {
      if( objStack_.empty() )
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
      assert( ! objStack_.empty() );
#ifdef USE_FINITE_STACK
      // finite stack does also return object on pop
      return objStack_.pop();
#else
      // std stack does not
      ObjectType * obj = objStack_.top();
      objStack_.pop();
      return obj;
#endif
    }

  private:
    //! do not copy pointers
    ALUMemoryProvider(const ALUMemoryProvider<Object> & org) {}
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
    if( objStack_.empty() )
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
    if( objStack_.empty() )
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
    while ( ! objStack_.empty() )
    {
      ObjectType * obj = stackObject();
      delete obj;
    }
  }

  template <class Object>
  inline void ALUMemoryProvider<Object>::freeObject(Object * obj)
  {
#ifdef USE_FINITE_STACK
    if( objStack_.full() )
      delete obj;
    else
#endif
    objStack_.push( obj );
  }

#undef USE_FINITE_STACK

} // end namespace Dune

#endif
