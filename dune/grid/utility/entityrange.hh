// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_ENTITYRANGE_HH
#define DUNE_GRID_UTILITY_ENTITYRANGE_HH

#include <iterator>

#include <dune/common/typetraits.hh>

namespace Dune {

  //! EntitySet category for simple entity ranges
  struct EntityRangeTag {};

  //! Interface for entity ranges
  /**
   * Validity of the iterators returned by \c begin() and \c end() is only
   * guaranteed as long as the EntityRange object itself is valid.
   */
  template<class Iterator_>
  struct EntityRangeInterface {
    //! category
    typedef EntityRangeTag EntitySetCategory;

    //! type of iterator
    typedef Iterator_ Iterator;
    //! type of entity
    typedef typename remove_const<
      typename std::iterator_traits<Iterator>::value_type
      >::type Entity;
    //! type used to count entites
    typedef typename std::iterator_traits<Iterator>::difference_type Size;

    //! begin iterator
    Iterator begin() const;
    //! end iterator
    Iterator end() const;

    //! number of entities in range
    /**
     * If possible, this should determine the size with better than O(n)
     * complexity.  If this is not possible, \c size() should beequivalent to
     * \c std::distance(begin(),end()).
     */
    Size size() const;
  };

  //! EntityRange defined by two iterators
  /**
   * \implements EntityRangeInterface
   */
  template<class Iterator_>
  class IteratorEntityRange {
  public:
    //! category
    typedef EntityRangeTag EntitySetCategory;

    //! type of iterator
    typedef Iterator_ Iterator;
    //! type of entity
    typedef typename remove_const<
      typename std::iterator_traits<Iterator>::value_type
      >::type Entity;
    //! type used to count entites
    typedef typename std::iterator_traits<Iterator>::difference_type Size;

    //! Constructor
    /**
     * The EntityRange stores copies of the iterators passed in here.
     */
    IteratorEntityRange(const Iterator &begin, const Iterator &end) :
      begin_(begin), end_(end)
    { }

    //! begin iterator
    const Iterator &begin() const
    {
      return begin_;
    }
    //! end iterator
    const Iterator &end() const
    {
      return end_;
    }

    //! number of entities in range
    /**
     * This does \c std::distance(begin(),end()) internally
     */
    Size size() const
    {
      return std::distance(begin_, end_);
    }

  private:
    Iterator begin_;
    Iterator end_;
  };

  template<int codim, class GV>
  IteratorEntityRange<typename GV::template Codim<codim>::Iterator>
  entityRange(const GV &gv)
  {
    return IteratorEntityRange<typename GV::template Codim<codim>::Iterator>
      (gv.template begin<codim>(), gv.template end<codim>());
  }

} // namespace Dune

#endif // DUNE_GRID_UTILITY_ENTITYRANGE_HH
