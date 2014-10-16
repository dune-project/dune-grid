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

    /** \brief postfix increment operator */
    EntityIterator operator++ (int)
    {
      EntityIterator tmp(*this);
      realIterator.increment();
      return tmp;
    }

    // The dereferencing operators are overridden here to avoid calling
    // the deprecated versions int the EntityPointer facade.

    /** \brief Dereferencing operator. */
    const Entity& operator*() const
    {
      return this->realIterator.dereference();
    }

    /** \brief Pointer operator. */
    const Entity * operator->() const
    {
      return & this->realIterator.dereference();
    }

    /** \brief Checks for equality.
            Only works for EntityPointers and iterators on the same grid.
            Due to the conversion operators one can compare
            all kinds of iterators and EntityPointer.

    bool operator== ( const EntityIterator< codim, Grid, ItImp > &rhs ) const
    {
      return this->equals( rhs );
    }

    /** \brief Checks for inequality.
            Only works for EntityPointers and iterators on the same grid.
            Due to the conversion operators one can compare
            all kinds of iterators and EntityPointer.
     */
    template< class ItImp >
    bool operator!= ( const EntityIterator< codim, Grid, ItImp > &rhs ) const
    {
      return !this->equals( rhs );
    }


    /** \brief Compares an EntityIterator with an Entity for equality.
     *
     * \deprecated This method only exists for backwards compatibility during the 2.4
     *             release cycle and will be removed after dune-grid-2.4 is released.
     */
    DUNE_DEPRECATED_MSG("EntityPointer is deprecated and will be removed after the release of dune-grid-2.4. \
Instead, you can copy and store entities directly now. You probably stumbled across this warning because you \
compared the return value of Entity::father(), Entity::subEntity(), Intersection::inside() or \
Intersection::outside() with an iterator. To fix the problem, derefence the iterator before comparing.")
    bool operator==(const Entity& rhs) const
    {
      return (**this) == rhs;
    }

    /** \brief Compares an EntityIterator with an Entity for inequality.
     *
     * \deprecated This method only exists for backwards compatibility during the 2.4
     *             release cycle and will be removed after dune-grid-2.4 is released.
     */
    DUNE_DEPRECATED_MSG("EntityPointer is deprecated and will be removed after the release of dune-grid-2.4. \
Instead, you can copy and store entities directly now. You probably stumbled across this warning because you \
compared the return value of Entity::father(), Entity::subEntity(), Intersection::inside() or \
Intersection::outside() with an iterator. To fix the problem, derefence the iterator before comparing.")
    bool operator!=(const Entity& rhs) const
    {
      return (**this) != rhs;
    }


    /** \name Implementor's interface
     *  \{
     */

    /** \brief default construct (undefined) iterator */
    EntityIterator ( )
    {}

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
    typedef const typename IteratorImp::Entity value_type;
    typedef value_type *pointer;
    typedef value_type &reference;
    typedef forward_iterator_tag iterator_category;
  };

} // namespace std

#endif // #ifndef DUNE_GRID_ENTITYITERATOR_HH
