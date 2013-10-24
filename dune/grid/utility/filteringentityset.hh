// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_FILTERINGENTITYSET_HH
#define DUNE_GRID_UTILITY_FILTERINGENTITYSET_HH

#include <cstddef>

#if HAVE_TBB
#include <tbb/tbb_stddef.h>
#endif

#include <dune/grid/common/entityiterator.hh>

#include "entityset.hh"
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

#if HAVE_TBB
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
#endif // HAVE_TBB

  private:
    Filter filter_;
    const GridView &gridView_;
  };

  template<class GV, int codim>
  class StridedEntitySet :
    public FilteringEntitySet<
      StridedIndexEntitySet<
        typename GV::IndexSet,
        typename GV::template Codim<codim>::Entity>,
      GV>
  {
    typedef StridedIndexEntitySet<
      typename GV::IndexSet,
      typename GV::template Codim<codim>::Entity> Filter;
    typedef FilteringEntitySet<Filter, GV> Base;

  public:
    StridedEntitySet(const GV &gv,
                     typename GV::IndexSet::IndexType maxStride = 0) :
      Base(Filter(gv.indexSet(), maxStride), gv)
    { }
#if HAVE_TBB
    StridedEntitySet(StridedEntitySet &other, tbb::split split) :
      Base(other, split)
    { }
#endif // HAVE_TBB
  };

  template<class GridView, int codim>
  class StridedPartitioning
  {
  public:
    typedef StridedEntitySet<GridView, codim> EntitySet;

    StridedPartitioning(const GridView &gv, std::size_t max_stride = 0) :
      gv_(gv), max_stride_(max_stride)
    { }

    std::size_t partitions() const
    {
      return max_stride_;
    }
    // return number of elements visited by partition
    std::size_t size(std::size_t partition) const
    {
      std::size_t totalsize = gv_.size(codim);
      return totalsize / max_stride_ + (partition < totalsize % max_stride_);
    }

    // entityset for whole partitioning
    EntitySet entitySet() const
    {
      return EntitySet(gv_, max_stride_);
    }

  private:
    const GridView &gv_;
    std::size_t max_stride_;
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_FILTERINGENTITYSET_HH
