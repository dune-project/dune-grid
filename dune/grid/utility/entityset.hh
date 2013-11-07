// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_ENTITYSET_HH
#define DUNE_GRID_UTILITY_ENTITYSET_HH

#include <algorithm>
#include <cstddef>

#if HAVE_TBB
#include <tbb/tbb_stddef.h>
#endif

#include <dune/common/documentation.hh>

#include <dune/geometry/type.hh>

namespace Dune {

  //! EntitySet category for simple entity sets
  struct EntitySetTag {};

  //! Interface for EntitySets
  /**
   * An EntitySet defines a set of entities.  The member function contains can
   * be used to check whether an entity is contained in the set.
   */
  template<class E>
  class EntitySetInterface
  {
  public:
    //! category
    typedef EntitySetTag EntitySetCategory;
    //! type of entity handled by this set.
    typedef E Entity;

    //! check whether an entity is contained in the set
    bool contains(const Entity &e) const;
  };

  //! EntitySet category for entity sets supporting \c size().
  struct SizableEntitySetTag : EntitySetTag {};

  //! Interface for SizableIterableEntitySets
  template<class E>
  class SizableEntitySetInterface :
    public EntitySetInterface<E>
  {
  public:
    //! category
    typedef SizableEntitySetTag EntitySetCategory;

    //! type used for the size
    typedef ImplementationDefined Size;

    //! return number of entites in the set
    Size size() const;
  };

  // sizable entity set adapted from an index set
  template<class IndexSet, class Entity_>
  class IndexSetEntitySet
  {
  public:
    //! category
    typedef SizableEntitySetTag EntitySetCategory;

    //! type of entity handled by this set.
    typedef Entity_ Entity;

    //! type used for the size
    typedef typename IndexSet::IndexType Size;

    //! construct
    /**
     * The reference to the index set is stored internally.
     */
    IndexSetEntitySet(const IndexSet &is) :
      is_(is)
    { }

    //! check whether an entity is contained in the set
    bool contains(const Entity &e) const
    {
      return is_.contains(e);
    }

    //! return number of entites in the set
    Size size() const
    {
      return is_.size(Entity::codimension);
    }

  private:
    const IndexSet &is_;
  };

  //! EntitySet selecting entity based on the index
  /**
   * Given an index set \c is, an offset \c offset and a stride \c stride,
   * this Entity set selects every entity \c e for which \c
   * is.index(e)%stride==offset is true.
   *
   * \implements SizableEntitySetInterface
   */
  template<class IndexSet, class E>
  class StridedIndexEntitySet
  {
  public:
    //! type used for counting entities
    typedef typename IndexSet::IndexType Size;

    //! category
    typedef SizableEntitySetTag EntitySetCategory;
    //! type of entity
    typedef E Entity;

    //! Construct
    /**
     * \param is        Index set.
     * \param stride    Stride.
     * \param offset    Offset.
     * \param maxStride This is used when splitting the EntitySet.  If
     *                  splitting would make \c stride larger than \c
     *                  maxStride, \c is_divisible() will be false.
     */
    StridedIndexEntitySet(const IndexSet &is, Size stride, Size offset,
                          Size maxStride = 0) :
      is_(&is), stride_(stride), offset_(offset), maxStride_(maxStride)
    { }

    //! Construct
    /**
     * \param is        Index set.
     * \param maxStride This is used when splitting the EntitySet.  If
     *                  splitting would make \c stride larger than \c
     *                  maxStride, \c is_divisible() will be false.
     *
     * This sets the initial stride to 1 and the initial offset to 0.
     */
    StridedIndexEntitySet(const IndexSet &is, Size maxStride = 0) :
      is_(&is), stride_(1), offset_(0), maxStride_(maxStride)
    { }

    //! check whether an entity is contained in the set
    bool contains(const Entity &e) const
    {
      return is_->index(e) % stride_ == offset_;
    }

    //! return number of entites in the set
    Size size() const
    {
      Size total = 0;
      for(auto gt : is_->geomTypes(Entity::codimension))
        total += sizePerGT(gt);
      return total;
    }
#if HAVE_TBB
    //! Splitting Constructor
    /**
     * Construct second half of set, update \c other to represent first half.
     *
     * This is done by doubling the stride for both sets, and using the offset
     * of \c other plus the old stride value as the offset of the constructed
     * set.
     */
    StridedIndexEntitySet(StridedIndexEntitySet &other, tbb::split) :
      is_(other.is_), stride_(other.stride_*2),
      offset_(other.offset_+other.stride_), maxStride_(other.maxStride_)
    {
      other.stride_ *= 2;
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
      if(maxStride_ && 2 * stride_ > maxStride_)
        return false;
      for(auto gt : is_->geomTypes(Entity::codimension))
        if(sizePerGT(gt) > 1)
          return true;
      return false;
    }
#endif // HAVE_TBB
  private:
    const IndexSet *is_;
    Size stride_;
    Size offset_;
    Size maxStride_;

    Size sizePerGT(const Dune::GeometryType &gt) const
    {
      Size size = is_->size(Entity::codimension);
      return size / stride_ + (offset_ < size % stride_);
    }
  };

  //! EntitySet of all entities for which a Vector contains a certain value
  template<class E, class Mapper, class Vector, class ID>
  class PartitionMapEntitySet
  {
    const Mapper &mapper_;
    const Vector &data_;
    ID id_;

  public:
    //! category
    typedef EntitySetTag EntitySetCategory;
    //! type of entity handled by this set.
    typedef E Entity;

    //! Construct
    /**
     * \param mapper Some object fulfilling Dune's Mapper interface.
     * \param data   Some subscriptable object that supports
     *               data[mapper.map(e)].
     * \param id     Value to look for in the vector.
     */
    PartitionMapEntitySet(const Mapper &mapper, const Vector &data, ID id) :
      mapper_(mapper), data_(data), id_(id)
    { }

    //! check whether an entity is contained in the set
    bool contains(const Entity &e) const
    {
      return data_[mapper_.map(e)] == id_;
    }
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_ENTITYSET_HH
