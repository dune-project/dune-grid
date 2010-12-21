// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_HIERARCHICITERATOR_HH
#define DUNE_GRID_HIERARCHICITERATOR_HH

#include <dune/grid/common/entityiterator.hh>

namespace Dune
{

  /**
     @brief Enables iteration over all codim zero entities
     in a subtree
     See also the documentation of Dune::EntityPointer.

     Mesh entities of codimension 0 ("elements") allow to visit all
     entities of codimension 0 obtained through nested, hierarchic
     refinement of the entity.  Iteration over this set of entities is
     provided by the HierarchicIterator, starting from a given entity.
     This is redundant but important for memory efficient
     implementations of unstructured hierarchically refined meshes.

     @ingroup GIEntityPointer
   */
  template<class GridImp, template<class> class HierarchicIteratorImp>
  class HierarchicIterator
    : public EntityIterator< 0, GridImp, HierarchicIteratorImp< GridImp > >
  {
    typedef EntityIterator< 0, GridImp, HierarchicIteratorImp< GridImp > > Base;

  public:
    /**
       @brief Preincrement operator.

       @note Forwarded to LevelIteratorImp.increment()
     */
    HierarchicIterator& operator++()
    {
      ++static_cast< Base & >( *this );
      return *this;
    }

    //===========================================================
    /** @name Implementor interface
     */
    //@{
    //===========================================================

    /** @brief copy constructor from HierarchicIteratorImp
     */
    HierarchicIterator (const HierarchicIteratorImp<const GridImp> & i)
      : Base( i )
    {}
    //@}
  };

}

#endif // DUNE_GRID_HIERARCHICITERATOR_HH
