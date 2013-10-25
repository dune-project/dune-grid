// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_PARTITIONING_MAPPED_HH
#define DUNE_GRID_UTILITY_PARTITIONING_MAPPED_HH

#include <algorithm>
#include <cstddef>
#include <vector>

#include <dune/geometry/type.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/utility/entityset.hh>
#include <dune/grid/utility/iterableentityset.hh>

namespace Dune {

  //! Partitioning based on a bitmap of the partitions
  /**
   * This stores the partition number in a vector indexed by a MCMGMapper.
   * All entities initially belong to partition 0, but the partition of an
   * entity may be changed using \c setPartitionId().
   */
  template<class GridView, int codim, class Size_ = std::size_t>
  class MappedPartitioning
  {
    template<int dimgrid>
    struct Layout {
      static bool contains(GeometryType gt) {
        return gt.dim() == dimgrid - codim;
      }
    };
    typedef MultipleCodimMultipleGeomTypeMapper<GridView, Layout> Mapper;
    typedef std::vector<Size_> Vector;
    typedef typename GridView::template Codim<0>::Entity Entity;
    typedef typename GridView::template Codim<0>::Iterator Iterator;
    typedef PartitionMapEntitySet<Entity, Mapper, Vector, Size_> Filter;

  public:
    //! type of partitions
    typedef IterableEntitySet<Filter, Iterator> Partition;
    //! type used to count partitions
    typedef Size_ Size;

    //! construct
    /**
     * Initially, all entites will belong to partition 0.
     */
    MappedPartitioning(const GridView &gv) :
      gv_(gv), mapper_(gv_), data_(mapper_.size(), 0)
    { }

    //! return maximum number of partitions
    /**
     * \note This iterates through the whole map to determine the maximum
     *       partition number used.
     */
    Size partitions() const
    {
      if(data_.size() == 0)
        return 0;
      else
        return *std::max_element(data_.begin(), data_.end())+1;
    }

    //! return a particular partition
    Partition partition(Size pId) const
    {
      return Partition(Filter(mapper_, data_, pId),
                       gv_.template begin<0>(), gv_.template end<0>());
    }

    //! get partition number of an entity
    Size getPartitionId(const Entity &e) const
    {
      return data_[mapper_.map(e)];
    }

    //! set partition number of an entity
    void setPartitionId(const Entity &e, Size pId)
    {
      data_[mapper_.map(e)] = pId;
    }

  private:
    GridView gv_;
    Mapper mapper_;
    std::vector<Size> data_;
  };


} // namespace Dune

#endif // DUNE_GRID_UTILITY_PARTITIONING_MAPPED_HH
