// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_HIERARCHICITERATOR_HH
#define DUNE_ALBERTA_HIERARCHICITERATOR_HH

#include <dune/grid/common/hierarchiciterator.hh>

#include <dune/grid/albertagrid/entitypointer.hh>

namespace Dune
{

  //**********************************************************************
  //
  // --AlbertaGridHierarchicIterator
  // --HierarchicIterator
  /*!
     Mesh entities of codimension 0 ("elements") allow to visit all entities of
     codimension 0 obtained through nested, hierarchic refinement of the entity.
     Iteration over this set of entities is provided by the HIerarchicIterator,
     starting from a given entity.
     This is redundant but important for memory efficient implementations of unstru
     hierarchically refined meshes.
   */

  template<class GridImp>
  class AlbertaGridHierarchicIterator
    : public AlbertaGridEntityPointer< 0, GridImp >
  {
    typedef AlbertaGridHierarchicIterator< GridImp > This;
    typedef AlbertaGridEntityPointer< 0, GridImp > Base;

  public:
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef typename GridImp::ctype ctype;

    typedef MakeableInterfaceObject< Entity > EntityObject;
    typedef typename EntityObject::ImplementationType EntityImp;

  public:
    //! the normal Constructor
    AlbertaGridHierarchicIterator ( const GridImp &grid,
                                    const Alberta::ElementInfo &elementInfo,
                                    int actLevel, int maxLevel, bool leafIt );

    //! the default Constructor
    AlbertaGridHierarchicIterator ( const GridImp &grid, int actLevel, int maxLevel );

    //! copy onstructor
    AlbertaGridHierarchicIterator ( const This &other );
    //! assignment operator
    This &operator= ( const This &other );

    //! increment
    void increment();

  protected:
    using Base::entityImp;

  private:
    void increment ( Alberta::ElementInfo elementInfo );

    //! make empty HierarchicIterator
    void makeIterator();

    int startLevel_;

    //! the actual Level of this Hierarichic Iterator
    int level_;

    //! max level to go down
    int maxlevel_;
  };

}

#endif
