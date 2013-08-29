// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_FILTERINGENTITYSET_HH
#define DUNE_GRID_UTILITY_FILTERINGENTITYSET_HH

#include <tbb/tbb_stddef.h>

#include <dune/grid/common/entityiterator.hh>

#include "entityfilter.hh"
#include "filteringentityiterator.hh"

namespace Dune {

  template<class Filter, class GV>
  class FilteringEntitySet
  {
  public:

    typedef GV GridView;
    enum {
      codim = Filter::Entity::codimension
    };

  private:
    typedef FilteringEntityIteratorImpl<
      Filter,
      typename GridView::template Codim<codim>::Iterator> IteratorImpl;

  public:
    //! Type of Elements contained in this EntitySet
    typedef typename Filter::Entity Element;

    //! Type of local coordinates with respect to the Element
    typedef typename Element::Geometry::LocalCoordinate LocalCoordinate;
    typedef typename Element::Geometry::GlobalCoordinate GlobalCoordinate;

    typedef Element value_type;

    //! A forward iterator
    typedef EntityIterator<codim, typename GridView::Grid, IteratorImpl>
      const_iterator;

    //! Same as const_iterator
    typedef const_iterator iterator;

    //! Construct FilteringEntitySet from scratch
    FilteringEntitySet(const Filter &filter, const GridView &gridView) :
      filter_(filter), gridView_(gridView)
    {}

    //! Returns true if e is contained in the EntitySet
    bool contains(const Element& e) const
    {
      return filter_.contains(e);
    }

    //! Number of Elements visited by an iterator
    size_t size() const
    {
      return filter_.size();
    }

    //! Create a begin iterator
    const_iterator begin() const
    {
      return IteratorImpl(filter_, gridView_.template begin<codim>(),
                          gridView_.template end<codim>());
    }

    //! Create a end iterator
    const_iterator end() const
    {
      return IteratorImpl(filter_, gridView_.template end<codim>());
    }

    const GridView& gridView() const
    {
      return gridView_;
    }

    // TBB range support
    FilteringEntitySet(FilteringEntitySet &other, tbb::split split) :
      filter_(other.filter_, split), gridView_(other.gridView_)
    { }
    bool empty() const
    {
      return filter_.empty();
    }
    bool is_divisible() const
    {
      return filter_.is_divisible();
    }

  private:
    Filter filter_;
    const GridView &gridView_;
  };

  template<class GV, int codim>
  class StridedEntitySet :
    public FilteringEntitySet<
      StridedEntityFilter<
        typename GV::IndexSet,
        typename GV::template Codim<codim>::Entity>,
      GV>
  {
    typedef StridedEntityFilter<
      typename GV::IndexSet,
      typename GV::template Codim<codim>::Entity> Filter;
    typedef FilteringEntitySet<Filter, GV> Base;

  public:
    StridedEntitySet(const GV &gv) :
      Base(Filter(gv.indexSet()), gv)
    { }
    StridedEntitySet(StridedEntitySet &other, tbb::split split) :
      Base(other, split)
    { }
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_FILTERINGENTITYSET_HH
