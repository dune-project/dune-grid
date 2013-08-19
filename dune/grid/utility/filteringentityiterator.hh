// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_FILTERINGENTITYITERATOR_HH
#define DUNE_GRID_UTILITY_FILTERINGENTITYITERATOR_HH

#include <cassert>

#include <dune/grid/common/entitypointer.hh>

namespace Dune {

  template<class Filter, class HostIterator>
  class FilteringEntityIteratorImpl
  {
    Filter filter_;
    HostIterator host_;
    HostIterator hostEnd_;

    void findValid()
    {
      while(host_ != hostEnd_ && !filter_.contains(*host_))
        ++host_;
    }

  public:
    enum { codimension = HostIterator::codimension };
    typedef typename HostIterator::Entity Entity;

    FilteringEntityIteratorImpl(const Filter &filter,
                                const HostIterator &host,
                                const HostIterator &hostEnd) :
      filter_(filter), host_(host), hostEnd_(hostEnd)
    {
      findValid();
    }

    FilteringEntityIteratorImpl(const Filter &filter,
                                const HostIterator &hostEnd) :
      filter_(filter), host_(hostEnd), hostEnd_(hostEnd)
    { }

    Entity &dereference() const
    {
      return *host_;
    }

    bool equals(const FilteringEntityIteratorImpl &rhs) const
    {
      return host_ == rhs.host_;
    }

    void increment() {
      assert(host_ != hostEnd_);
      ++host_;
      findValid();
    };
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_FILTERINGENTITYITERATOR_HH
