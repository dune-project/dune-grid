// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_ENTITYITERATOR_HH
#define DUNE_GRID_COMMON_ENTITYITERATOR_HH

#include <cstddef>
#include <iterator>

namespace Dune
{

  /** \class EntityIterator
   *  \brief interface class for an iterator over grid entities
   *
   *  An entity iterator is an iterator over a subset of entities within a
   *  hierarchical grid.
   *
   *  Examples of entity iterators are:
   *  - iterators over the leaf level (LeafGridView::Iterator)
   *  - iterators over a grid level (LevelGridView::Iterator)
   *  - iterators over the children of an entity (Grid::HierarchicIterator)
   *  .
   *
   *  \tparam  codim        codimension of entities this iterator walks over
   *  \tparam  Grid         type of the grid implementation
   *  \tparam  IteratorImp  type of the iterator implementation
   */
  template< int codim, class Grid, class IteratorImp >
  class EntityIterator
  {
  protected:
    IteratorImp realIterator;

  public:
    /**
     * \brief type of underlying implementation
     *
     * \warning Implementation details may change without prior notification.
     **/
    typedef IteratorImp Implementation;

    /**
     * \brief access to the underlying implementation
     *
     * \warning Implementation details may change without prior notification.
     **/
    Implementation &impl () { return realIterator; }
    /**
     * \brief access to the underlying implementation
     *
     * \warning Implementation details may change without prior notification.
     **/
    const Implementation &impl () const { return realIterator; }

    typedef typename Grid::template Codim< codim >::Entity Entity;

    /** \brief Type of the reference used when derefencing the Ptr */
    typedef typename std::conditional<
      std::is_lvalue_reference<
        decltype(realIterator.dereference())
        >::value,
      const Entity&,
      Entity
      >::type Reference;

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

    // The behavior when dereferencing the EntityIterator facade depends on
    // the way the grid implementation handles returning entities. The implementation
    // may either return a reference to an entity stored inside the EntityIterator
    // implementation or a temporary Entity object. This object has to be forwarded through
    // the facade to the user, which requires a little trickery, especially for operator->().
    //
    // In order to avoid confusing users reading the Doxygen documentation, we provide "clean"
    // function signatures to Doxygen and hide the actual implementations.


#ifdef DOXYGEN

    /** \brief Dereferencing operator. */
    const Entity& operator*() const;

    /** \brief Pointer operator. */
    const Entity& operator->() const;

#else // DOXYGEN

    /** \brief Dereferencing operator. */
    typename std::conditional<
      std::is_lvalue_reference<
        decltype(realIterator.dereference())
        >::value,
      const Entity&,
      Entity
      >::type
    operator*() const
    {
      return realIterator.dereference();
    }

    /** \brief Pointer operator. */
    decltype(handle_proxy_member_access(realIterator.dereference()))
    operator->() const
    {
      return handle_proxy_member_access(realIterator.dereference());
    }

#endif // DOXYGEN


    /** \brief Checks for equality. */
    bool operator==(const EntityIterator& rhs) const
    {
      return this->realIterator.equals(rhs.realIterator);
    }

    /** \brief Checks for inequality. */
    bool operator!=(const EntityIterator& rhs) const
    {
      return !this->realIterator.equals(rhs.realIterator);
    }


    /** \name Implementor's interface
     *  \{
     */

    /** \brief default construct (undefined) iterator */
    EntityIterator ( )
    {}

    /** \brief copy constructor from implementaton */
    EntityIterator ( const IteratorImp &imp )
      : realIterator( imp )
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

#endif // #ifndef DUNE_GRID_COMMON_ENTITYITERATOR_HH
