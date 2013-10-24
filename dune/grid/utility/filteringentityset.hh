// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_FILTERINGENTITYSET_HH
#define DUNE_GRID_UTILITY_FILTERINGENTITYSET_HH

#include <cstddef>

#if HAVE_TBB
#include <tbb/tbb_stddef.h>
#endif

#include <dune/grid/utility/iteratoradapters.hh>
#include <dune/grid/utility/iterableentityset.hh>

#include "entityset.hh"

namespace Dune {

  template<class GV, int codim>
  class StridedEntitySet :
    public IterableEntitySet<
      StridedIndexEntitySet<
        typename GV::IndexSet,
        typename GV::template Codim<codim>::Entity>,
      typename GV::template Codim<codim>::Iterator>
  {
    typedef StridedIndexEntitySet<
      typename GV::IndexSet,
      typename GV::template Codim<codim>::Entity> Filter;
    typedef IterableEntitySet<
      Filter, typename GV::template Codim<codim>::Iterator> Base;

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
