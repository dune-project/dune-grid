// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_ENTITYITERATOR_HH
#define DUNE_GRID_ENTITYITERATOR_HH

#include <cstddef>
#include <iterator>

#include <dune/grid/common/entitypointer.hh>

namespace Dune
{

  /** \class EntityIterator
   *  \brief interface class for an iterator over grid entities
   *  \ingroup GIEntityPointer
   *
   *  An entity iterator is an iterator over a subset of entities within a
   *  hierarchical grid. It is an extension of the Dune::EntityPointer
   *  interface.
   *
   *  Examples of entity iterators are:
   *  - iterators over the leaf level (LeafGridView::Iterator)
   *  - iterators over a grid level (LevelGridView::Iterator)
   *  - iterators over the children of an entity (Grid::HierarchicIterator)
   *  .
   *
   *  See also the documentation of Dune::EntityPointer.
   *
   *  \tparam  codim        codimension of entities this iterator walks over
   *  \tparam  Grid         type of the grid implementation
   *  \tparam  IteratorImp  type of the iterator implementation
   */
  template< int codim, class Grid, class IteratorImp >
  class EntityIterator
    : public EntityPointer< Grid, IteratorImp >
  {
    typedef EntityPointer< Grid, IteratorImp > Base;

  protected:
    using Base::realIterator;

  public:
    typedef typename Grid::template Codim< codim >::Entity Entity;

    /** \brief prefix increment operator */
    EntityIterator &operator++ ()
    {
      realIterator.increment();
      return *this;
    }

    /** \name Implementor's interface
     *  \{
     */

    /** \brief copy constructor from implementaton */
    EntityIterator ( const IteratorImp &imp )
      : Base( imp )
    {}

    /** \} */
  };

} // namespace Dune

namespace std
{

  template< int codim, class Grid, class IteratorImp >
  struct iterator_traits< Dune::EntityIterator< codim, Grid, IteratorImp > >
  {
    typedef ptrdiff_t difference_type;
    typedef const typename Dune::EntityIterator< codim, Grid, IteratorImp >::Entity value_type;
    typedef value_type *pointer;
    typedef value_type &reference;
    typedef forward_iterator_tag iterator_category;
  };

} // namespace std

#endif // #ifndef DUNE_GRID_ENTITYITERATOR_HH
