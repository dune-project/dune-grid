// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_PARTITIONING_STRIDED_HH
#define DUNE_GRID_UTILITY_PARTITIONING_STRIDED_HH

#include <dune/grid/utility/entityset.hh>
#include <dune/grid/utility/iterableentityset.hh>

namespace Dune {

  //! Partitioning selecting strides based on the index
  template<class GridView, int codim>
  class StridedPartitioning
  {
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<codim>::Entity Entity;
    typedef typename GridView::template Codim<codim>::Iterator Iterator;
    typedef StridedIndexEntitySet<IndexSet, Entity> Filter;

  public:
    //! type of partitions
    typedef IterableEntitySet<Filter, Iterator> Partition;
    //! type used to count partitions
    typedef typename Filter::Size Size;

    //! construct
    /**
     * \param max_stride Maximum stride to allow when splitting.
     */
    StridedPartitioning(const GridView &gv, Size max_stride = 0) :
      gv_(gv), max_stride_(max_stride)
    { }

    //! return maximum number of partitions
    Size partitions() const
    {
      return max_stride_;
    }

    //! return a particular partition
    Partition partition(Size pId) const
    {
      return Partition(Filter(gv_.indexSet(), max_stride_, pId, max_stride_),
                       gv_.template begin<codim>(), gv_.template end<codim>());
    }

    //! whole partitioning as a partition object
    /**
     * The partition can be split by tbb
     */
    Partition everything() const
    {
      return Partition(Filter(gv_.indexSet(), max_stride_),
                       gv_.template begin<codim>(), gv_.template end<codim>());
    }

  private:
    const GridView &gv_;
    Size max_stride_;
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_PARTITIONING_STRIDED_HH
