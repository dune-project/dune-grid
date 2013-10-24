// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_ITERABLEENTITYSET_HH
#define DUNE_GRID_UTILITY_ITERABLEENTITYSET_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/utility/entityrange.hh>
#include <dune/grid/utility/entityset.hh>

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

} // namespace Dune

#endif // DUNE_GRID_UTILITY_ITERABLEENTITYSET_HH
