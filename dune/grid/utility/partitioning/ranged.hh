// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_PARTITIONING_RANGED_HH
#define DUNE_GRID_UTILITY_PARTITIONING_RANGED_HH

#include <iterator>
#include <vector>

#if HAVE_TBB
#include <tbb/tbb_stddef.h>
#endif

#include <dune/grid/utility/entityrange.hh>

namespace Dune {

  //! Partioning base on remambering iterator ranges
  template<class GridView, int codim>
  class RangedPartitioning
  {
    typedef typename GridView::template Codim<codim>::Iterator Iterator;
  public:
    //! type of partitions
    class Partition;
    //! type used to count partitions
    typedef typename std::iterator_traits<
      typename std::vector<Iterator>::iterator
      >::difference_type Size;

    //! construct
    RangedPartitioning(const GridView &gv, Size partitions)
    {
      const auto end = gv.template end<0>();
      // GridView's size() would return the wrong count for non AllPartition Views/Iterators
      const auto totalsize = std::distance(gv.template begin<0>(), end);
      auto it = gv.template begin<0>();
      entry_points_.reserve(partitions+1);
      stride_ = totalsize / partitions;
      overflow_ = totalsize % partitions;
      for(Size p = 0; p < partitions; ++p)
      {
        entry_points_.push_back(it);
        std::advance(it, stride_ + (p < overflow_));
      }
      entry_points_.push_back(end);
    }

    //! return maximum number of partitions
    Size partitions() const
    {
      return entry_points_.size()-1;
    }

    //! whole partitioning as a partition object
    /**
     * The partition can be split by tbb
     */
    Partition everything() const
    {
      return Partition(*this, 0, partitions());
    }
    //! return a particular partition
    Partition partition(Size pId) const
    {
      return Partition(*this, pId, pId+1);
    }
    //! return a range of partitions
    /**
     * The returned partition object can be split by tbb
     */
    Partition partitions(Size first, Size last) const
    {
      return Partition(*this, first, last);
    }

  private:
    std::vector<Iterator> entry_points_;
    Size stride_;
    Size overflow_;
  };

  //! Partition type for ranged partitioning
  /**
   * \implements EntityRangeInterface
   *
   * This implements the interface for entity ranges.  It differs from
   * IteratorEntityRange in that it is TBB-splittable and thus needs to hold a
   * reference to the partitioning.
   */
  template<class GridView, int codim>
  class RangedPartitioning<GridView, codim>::Partition
  {
  public:
    //! category
    typedef EntityRangeTag EntitySetCategory;

    //! type of iterator
    typedef typename GridView::template Codim<codim>::Iterator Iterator;
    //! type of entity
    typedef typename GridView::template Codim<codim>::Entity Entity;
    //! type used to count entites
    typedef typename std::iterator_traits<Iterator>::difference_type Size;

    //! Construct partition from scratch
    Partition(const RangedPartitioning &partitioning,
              RangedPartitioning::Size firstPartition,
              RangedPartitioning::Size lastPartition) :
      partitioning_(partitioning), firstPartition_(firstPartition),
      lastPartition_(lastPartition)
    {}

    //! Create a begin iterator
    const Iterator &begin() const
    {
      return partitioning_.entry_points_[firstPartition_];
    }

    //! Create an end iterator
    const Iterator &end() const
    {
      return partitioning_.entry_points_[lastPartition_];
    }

    //! Number of Elements visited by an iterator
    Size size() const
    {
      return (lastPartition_-firstPartition_)*partitioning_.stride_
        + (firstPartition_ < partitioning_.overflow_
           ? partitioning_.overflow_ - firstPartition_
           : 0);
    }

#if HAVE_TBB
    //! Splitting Constructor
    /**
     * Construct second half of set, update \c other to represent first half.
     *
     * This is done by looking up the entry point in the partitioning that
     * lies roughly in the middle of the current \c begin() and \c end().
     */
    Partition(Partition &other, tbb::split) :
      partitioning_(other.partitioning_)
    {
      firstPartition_ = (other.firstPartition_ + other.lastPartition_ + 1) / 2;
      lastPartition_ = other.lastPartition_;
      other.lastPartition_ = firstPartition_;
    }
    //! check whether set is empty (for TBB)
    /**
     * This is equivalent to \c size()==0.
     */
    bool empty() const
    {
      return size() == 0;
    }
    //! check whether set can be split (for TBB)
    bool is_divisible() const
    {
      return lastPartition_ - firstPartition_ > 1;
    }
#endif // HAVE_TBB

  private:
    const RangedPartitioning &partitioning_;
    RangedPartitioning::Size firstPartition_;
    RangedPartitioning::Size lastPartition_;
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_PARTITIONING_RANGED_HH
