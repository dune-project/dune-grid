// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_RANGEDENTITYSET_HH
#define DUNE_GRID_UTILITY_RANGEDENTITYSET_HH

#include <cstddef>
#include <iterator>
#include <vector>

#if HAVE_TBB
#include <tbb/tbb_stddef.h>
#endif

namespace Dune {

  template<class GridView, int codim>
  class RangedPartitioning
  {
  public:
    class EntitySet;

    //! A forward iterator
    typedef typename GridView::template Codim<codim>::Iterator const_iterator;

    RangedPartitioning(const GridView &gv, std::size_t partitions) :
      gv_(gv)
    {
      auto totalsize = gv.size(codim);
      auto end = gv.template end<0>();
      auto it = gv.template begin<0>();
      entry_points_.reserve(partitions+1);
      stride_ = totalsize / partitions;
      overflow_ = totalsize % partitions;
      for(std::size_t p = 0; p < partitions; ++p)
      {
        entry_points_.push_back(it);
        std::advance(it, size(p));
      }
      entry_points_.push_back(end);
    }

    std::size_t partitions() const
    {
      return entry_points_.size()-1;
    }
    // return number of elements visited by partition
    std::size_t size(std::size_t partition) const
    {
      return stride_ + (partition<overflow_);
    }
    // return number of elements visited by a range of partitions
    std::size_t size(std::size_t first, std::size_t last) const
    {
      return (last-first)*stride_
        + (first < overflow_ ? overflow_ - first : 0);
    }

    // entityset for whole partitioning
    EntitySet entitySet() const
    {
      return EntitySet(*this, 0, partitions());
    }
    // entityset for a particular partition
    EntitySet entitySet(std::size_t partition) const
    {
      return EntitySet(*this, partition, partition+1);
    }
    // entityset for a range of partitions
    EntitySet entitySet(std::size_t first, std::size_t last) const
    {
      return EntitySet(*this, first, last);
    }

  private:
    const GridView &gv_;
    std::vector<const_iterator> entry_points_;
    std::size_t stride_;
    std::size_t overflow_;
  };

  template<class GV, int cd>
  class RangedPartitioning<GV, cd>::EntitySet
  {
  public:

    typedef GV GridView;
    enum {
      codim = cd
    };

    //! Type of Elements contained in this EntitySet
    typedef typename GridView::template Codim<codim>::Entity Element;

    //! Type of local coordinates with respect to the Element
    typedef typename Element::Geometry::LocalCoordinate LocalCoordinate;
    typedef typename Element::Geometry::GlobalCoordinate GlobalCoordinate;

    typedef Element value_type;

    //! A forward iterator
    typedef typename GridView::template Codim<codim>::Iterator const_iterator;

    //! Same as const_iterator
    typedef const_iterator iterator;

    //! Construct RangedEntitySet from scratch
    EntitySet(const RangedPartitioning &partitioning,
              std::size_t firstPartition, std::size_t lastPartition) :
      partitioning_(partitioning), firstPartition_(firstPartition),
      lastPartition_(lastPartition)
    {}

    // //! Returns true if e is contained in the EntitySet
    // // not implemented
    // template<class Element>
    // bool contains(const Element& e) const
    // {
    // }

    //! Number of Elements visited by an iterator
    size_t size() const
    {
      return partitioning_.size(firstPartition_, lastPartition_);
    }

    //! Create a begin iterator
    const const_iterator &begin() const
    {
      return partitioning_.entry_points_[firstPartition_];
    }

    //! Create a end iterator
    const const_iterator &end() const
    {
      return partitioning_.entry_points_[lastPartition_];
    }

    const GridView& gridView() const
    {
      return partitioning_.gv_;
    }

#if HAVE_TBB
    // TBB range support
    EntitySet(EntitySet &other, tbb::split) :
      partitioning_(other.partitioning_)
    {
      firstPartition_ = (other.firstPartition_ + other.lastPartition_ + 1) / 2;
      lastPartition_ = other.lastPartition_;
      other.lastPartition_ = firstPartition_;
    }
    bool empty() const
    {
      return size() == 0;
    }
    bool is_divisible() const
    {
      return lastPartition_ - firstPartition_ > 1;
    }
#endif // HAVE_TBB

  private:
    RangedPartitioning partitioning_;
    std::size_t firstPartition_;
    std::size_t lastPartition_;
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_RANGEDENTITYSET_HH
