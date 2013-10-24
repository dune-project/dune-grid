// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_ITERABLEENTITYSET_HH
#define DUNE_GRID_UTILITY_ITERABLEENTITYSET_HH

#include <iterator>

#if HAVE_TBB
#include <tbb/tbb_stddef.h>
#endif

#include <dune/common/typetraits.hh>

#include <dune/grid/utility/entityrange.hh>
#include <dune/grid/utility/entityset.hh>
#include <dune/grid/utility/iteratoradapters.hh>

namespace Dune {

  //! EntitySet category for entity sets supporting \c size().
  struct IterableEntitySetTag : SizableEntitySetTag {};

  //! Interface for IterableEntitySets
  /**
   * This basically merges the interface of EntitySets and entityRanges.  \c
   * contains() should return true for every entity visited by iteration, and
   * no other entities.  \c size() should be equal to \c
   * std::distance(begin(),end()).  \c Size should be the same type as \c
   * std::iterator_traits<Iterator>::difference_type.
   */
  template<class Iterator_>
  class IterableEntitySetInterface :
    public SizableEntitySetInterface<
      typename remove_const<typename Iterator_::value_type>::type>,
    public EntityRangeInterface<Iterator_>
  {
  public:
    //! category
    typedef IterableEntitySetTag EntitySetCategory;

    using EntityRangeInterface<Iterator_>::Entity;
  };

  //! generic iterable entityset given a plain entity set and an iterator range
  /**
   * This implements iteration by iterating over the iterator range and
   * stopping only when \c \c contains() is true.
   */
  template<class EntitySet, class HostIterator>
  class IterableEntitySet
  {
  public:
    //! category
    typedef IterableEntitySetTag EntitySetCategory;

    //! type of entity
    typedef typename remove_const<
      typename std::iterator_traits<HostIterator>::value_type>::type Entity;
    //! Iterator
    typedef EntitySetIterator<EntitySet, HostIterator> Iterator;
    //! type of Size
    typedef typename std::iterator_traits<HostIterator>::difference_type Size;

    //! Construct
    /**
     * This stores copies of the passed-in values.
     */
    IterableEntitySet(const EntitySet &entitySet, const HostIterator &begin,
                      const HostIterator &end) :
      entitySet_(entitySet), begin_(begin), end_(end)
    { }

    //! Returns true if e is contained in the EntitySet
    bool contains(const Entity& e) const
    {
      return entitySet_.contains(e);
    }

    //! Number of Elements visited by an iterator
    Size size() const
    {
      return sizeImpl(EntitySet::EntitySetCategory());
    }

    //! Create a begin iterator
    Iterator begin() const
    {
      return Iterator(entitySet_, begin_, end_);
    }

    //! Create a end iterator
    Iterator end() const
    {
      return Iterator(entitySet_, end_);
    }

#if HAVE_TBB
    // TBB range support
    IterableEntitySet(IterableEntitySet &other, tbb::split split) :
      entitySet_(other.entitySet_, split), begin_(other.begin_),
      end_(other.end_)
    { }
    bool empty() const
    {
      return entitySet_.empty();
    }
    bool is_divisible() const
    {
      return entitySet_.is_divisible();
    }
#endif // HAVE_TBB

  private:
    Size sizeImpl(EntitySetTag) const
    {
      return std::distance(begin(),end());
    }
    Size sizeImpl(SizableEntitySetTag) const
    {
      return entitySet_.size();
    }

    EntitySet entitySet_;
    HostIterator begin_;
    HostIterator end_;
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_ITERABLEENTITYSET_HH
