// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU3DGRIDMEMORY_HH
#define DUNE_ALU3DGRIDMEMORY_HH

#include <cassert>
#include <cstdlib>
#include <vector>

#include <dune/common/finitestack.hh>

namespace Dune {

  //! organize the memory management for entitys used by the NeighborIterator
  template <class Object>
  class ALUMemoryProvider
  {
    enum { maxStackObjects = 256 };
    typedef FiniteStack< Object* , maxStackObjects > StackType ;

    StackType objStack_ ;

    typedef ALUMemoryProvider < Object > MyType;

    StackType& objStack() { return objStack_ ; }

  public:
    typedef Object ObjectType;

    //!default constructor
    ALUMemoryProvider() {}

    //! do not copy pointers
    ALUMemoryProvider(const ALUMemoryProvider<Object> & org)
      : objStack_( org.objStack_ )
    {}

    //! call deleteEntity
    ~ALUMemoryProvider ();

    //! i.e. return pointer to Entity
    template <class FactoryType>
    ObjectType * getObject(const FactoryType &factory, int level);

    //! i.e. return pointer to Entity
    template <class FactoryType, class EntityImp>
    inline ObjectType * getEntityObject(const FactoryType& factory, int level , EntityImp * fakePtr )
    {
      if( objStack().empty() )
      {
        return ( new ObjectType(EntityImp(factory,level) ));
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

  };


  //************************************************************************
  //
  //  ALUMemoryProvider implementation
  //
  //************************************************************************
  template <class Object> template <class FactoryType>
  inline typename ALUMemoryProvider<Object>::ObjectType *
  ALUMemoryProvider<Object>::getObject
    (const FactoryType &factory, int level )
  {
    if( objStack().empty() )
    {
      return ( new Object (factory, level) );
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
    StackType& objStk = objStack_;
    while ( ! objStk.empty() )
    {
      ObjectType * obj = objStk.pop();
      delete obj;
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
