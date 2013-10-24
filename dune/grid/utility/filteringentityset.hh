// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_FILTERINGENTITYSET_HH
#define DUNE_GRID_UTILITY_FILTERINGENTITYSET_HH

#include <cstddef>

#include <dune/grid/utility/entityset.hh>
#include <dune/grid/utility/iterableentityset.hh>

namespace Dune {

  template<class GridView, int codim>
  class StridedPartitioning
  {
    typedef StridedIndexEntitySet<
      typename GridView::IndexSet,
      typename GridView::template Codim<codim>::Entity> Filter;

  public:
    typedef IterableEntitySet<
      Filter,
      typename GridView::template Codimy<codim>::Iterator> EntitySet;

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
      return EntitySet(Filter(gv_.indexSet(), max_stride_),
                       gv_.template begin<codim>(), gv_.template end<codim>());
    }

  private:
    const GridView &gv_;
    std::size_t max_stride_;
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_FILTERINGENTITYSET_HH
