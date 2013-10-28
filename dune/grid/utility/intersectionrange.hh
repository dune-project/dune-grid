// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_INTERSECTIONRANGE_HH
#define DUNE_GRID_UTILITY_INTERSECTIONRANGE_HH

#include <iterator>

#include <dune/common/typetraits.hh>

namespace Dune {

  //! IntersectionRange defined by two iterators
  template<class Iterator_>
  class IteratorIntersectionRange {
  public:
    //! type of iterator
    typedef Iterator_ Iterator;
    //! type of intersection
    typedef typename remove_const<
      typename std::iterator_traits<Iterator>::value_type
      >::type Intersection;

    //! Constructor
    /**
     * The IntersectionRange stores copies of the iterators passed in here.
     */
    IteratorIntersectionRange(const Iterator &begin, const Iterator &end) :
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

  private:
    Iterator begin_;
    Iterator end_;
  };

  template<class GV>
  IteratorIntersectionRange<typename GV::IntersectionIterator>
  intersections(const GV &gv, const typename GV::template Codim<0>::Entity &e)
  {
    return IteratorIntersectionRange<typename GV::IntersectionIterator>
      (gv.ibegin(e), gv.iend(e));
  }

} // namespace Dune

#endif // DUNE_GRID_UTILITY_INTERSECTIONRANGE_HH
