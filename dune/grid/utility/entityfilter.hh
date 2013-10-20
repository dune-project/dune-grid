// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_ENTITYFILTER_HH
#define DUNE_GRID_UTILITY_ENTITYFILTER_HH

#include <algorithm>
#include <cstddef>

#if HAVE_TBB
#include <tbb/tbb_stddef.h>
#endif

#include <dune/geometry/type.hh>

namespace Dune {

  template<class E>
  class EntityFilterInterface
  {
  public:
    typedef E Entity;

    bool contains(const Entity &e) const;
#if HAVE_TBB
    // TBB range support
    EntityFilterInterface(EntityFilterInterface &, tbb::split);
    bool empty() const;
    bool is_divisible() const;
#endif // HAVE_TBB
  };

  template<class IndexSet, class E>
  class StridedEntityFilter
  {
    const IndexSet *is_;
    typename IndexSet::IndexType stride_;
    typename IndexSet::IndexType offset_;
    typename IndexSet::IndexType maxStride_;

    typename IndexSet::IndexType sizePerGT(const Dune::GeometryType &gt) const
    {
      typename IndexSet::IndexType size = is_->size(Entity::codimension);
      return size / stride_ + (offset_ < size % stride_);
    }

  public:
    typedef E Entity;

    StridedEntityFilter(const IndexSet &is,
                        typename IndexSet::IndexType stride,
                        typename IndexSet::IndexType offset,
                        typename IndexSet::IndexType maxStride = 0) :
      is_(&is), stride_(stride), offset_(offset), maxStride_(maxStride)
    { }

    StridedEntityFilter(const IndexSet &is,
                        typename IndexSet::IndexType maxStride = 0) :
      is_(&is), stride_(1), offset_(0), maxStride_(maxStride)
    { }

    bool contains(const Entity &e) const
    {
      return is_->index(e) % stride_ == offset_;
    }

    typename IndexSet::IndexType size() const
    {
      typename IndexSet::IndexType total = 0;
      for(auto gt : is_->geomTypes(Entity::codimension))
        total += sizePerGT(gt);
      return total;
    }
#if HAVE_TBB
    // TBB range support
    // construct second half of set, update other to represent first half
    StridedEntityFilter(StridedEntityFilter &other, tbb::split) :
      is_(other.is_), stride_(other.stride_*2),
      offset_(other.offset_+other.stride_), maxStride_(other.maxStride_)
    {
      other.stride_ *= 2;
    }
    bool empty() const
    {
      return size() == 0;
    }
    bool is_divisible() const
    {
      if(maxStride_ && 2 * stride_ > maxStride_)
        return false;
      for(auto gt : is_->geomTypes(Entity::codimension))
        if(sizePerGT(gt) > 1)
          return true;
      return false;
    }
#endif // HAVE_TBB
  };

  template<class E, class Mapper, class Vector, class ID>
  class GeneralEntityFilter
  {
    const Mapper &mapper_;
    const Vector &data_;
    ID id_;

  public:
    typedef E Entity;

    GeneralEntityFilter(const Mapper &mapper, const Vector &data, ID id) :
      mapper_(mapper), data_(data), id_(id)
    { }

    bool contains(const Entity &e) const
    {
      return data_[mapper_.map(e)] == id_;
    }

    std::size_t size() const
    {
      return std::count(data_.begin(), data_.end(), id_);
    }
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_ENTITYFILTER_HH
