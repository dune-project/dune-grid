// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_TRIPARTIT_HH
#define DUNE_GRID_UTILITY_TRIPARTIT_HH

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/utility/entityfilter.hh>
#include <dune/grid/utility/filteringentityset.hh>

namespace Dune {

  template<class GV, class Vector, class ID,
           class Mapper = MultipleCodimMultipleGeomTypeMapper
                            <GV, MCMGElementLayout> >
  class GeneralFilteredEntitySet :
    public FilteringEntitySet<
      GeneralEntityFilter<typename GV::template Codim<0>::Entity,
                          Mapper, Vector, ID>,
      GV>
  {
  public:
    typedef typename GV::template Codim<0>::Entity Entity;

  private:
    typedef GeneralEntityFilter<Entity, Mapper, Vector, ID> Filter;
    typedef FilteringEntitySet<Filter, GV> Base;

  public:
    GeneralFilteredEntitySet(const GV &gv, const Vector &data, ID id,
                             const Mapper &mapper) :
      Base(Filter(mapper, data, id), gv)
    { }
  };

  template<class GridView>
  class GeneralFilteredPartitioning
  {
    typedef MultipleCodimMultipleGeomTypeMapper<GridView, MCMGElementLayout>
      Mapper;
    typedef std::vector<std::size_t> Vector;

  public:
    typedef GeneralFilteredEntitySet<GridView, Vector, std::size_t,
                                     Mapper> EntitySet;
    typedef typename EntitySet::Entity Element;

    GeneralFilteredPartitioning(const GridView &gv) :
      gv_(gv), mapper_(gv_), data_(mapper_.size(), 0)
    { }

    std::size_t getPartition(const Element &e) const
    {
      return data_[mapper_.map(e)];
    }

    void setPartition(const Element &e, std::size_t p)
    {
      data_[mapper_.map(e)] = p;
    }

    std::size_t partitions() const
    {
      if(data_.size() == 0)
        return 0;
      else
        return *std::max_element(data_.begin(), data_.end())+1;
    }
    // return number of elements visited by partition
    std::size_t size(std::size_t partition) const
    {
      return std::count(data_.begin(), data_.end(), partition);
    }

    // entityset for whole partitioning
    EntitySet entitySet(std::size_t partition) const
    {
      return EntitySet(gv_, data_, partition, mapper_);
    }

  private:
    GridView gv_;
    Mapper mapper_;
    std::vector<std::size_t> data_;
  };

  class Coloring
  {
    Coloring(std::size_t ncolors) :
      color_partitions(ncolors)
    { }

    void set_colors(std::size_t ncolors)
    {
      color_partitions.clear();
      color_partitions.resize(ncolors);
    }

    const std::vector<std::size_t> &getPartitions(std::size_t color) const
    {
      return color_partitions[color];
    }

    void add_partition(std::size_t color, std::size_t partition)
    {
      color_partitions[color].push_back(partition);
    }

  private:
    std::vector<std::vector<std::size_t> > color_partitions;
  };

  struct OverlapError : Exception {};

  template<class GV, class Partitioning, class Mapper>
  class EquidistantPartitioner
  {
    static const std::size_t nooverlap =
      std::numeric_limits<std::size_t>::max();
    static const std::size_t someoverlap = nooverlap - 1;

    const GV &gv_;
    Partitioning &partitioning_;
    const Mapper &mapper_;
    std::vector<std::size_t> overlapMap_;
    std::size_t direction_;
    std::size_t overlap_;

    template<class EntitySet>
    class PartitioningContext {
      typedef typename EntitySet::Entity Entity;
      typedef typename Entity::ctype ctype;

      const GV &gv_;
      Partitioning &partitioning_;
      const Mapper &mapper_;
      std::vector<std::size_t> &overlapMap_;
      std::size_t direction_;
      std::size_t overlap_;
      const EntitySet &entitySet_;
      const std::vector<std::size_t> &newPartitionIds_;

      ctype minc_;
      ctype maxc_;

      void initMinMax()
      {
        using std::min;
        using std::max;
        minc_ = std::numeric_limits<ctype>::infinity();
        maxc_ = -minc_;
        for(const auto &e : entitySet_)
        {
          ctype p = e.geometry().center()[direction_];
          minc_ = min(minc_, p);
          maxc_ = max(maxc_, p);
        }
      }

      std::size_t pIndex(const Entity &e) const
      {
        using std::min;
        using std::max;
        ctype tmp = e.geometry().center()[direction_];
        tmp -= minc_;
        tmp /= maxc_ - minc_;
        tmp *= newPartitionIds_.size();
        tmp = std::floor(tmp);
        tmp = max(ctype(0), tmp);
        return min(newPartitionIds_.size()-1, std::size_t(tmp));
      };

      void markneighbors(const Entity &e, std::size_t myPIndex) {
        auto iend = gv_.iend(e);
        for(auto iit = gv_.ibegin(e); iit != iend; ++iit)
        {
          const auto &is = *iit;
          if(!is.neighbor())
            continue;
          auto outsidep = is.outside();
          const auto &outside = *outsidep;
          if(!entitySet_.contains(outside))
            // These were marked previously in a coarser partitioning
            continue;
          // Now it is legal to obtain the local index for this partition
          std::size_t outPIndex = pIndex(outside);
          if(outPIndex == myPIndex)
            // Only mark our overlap in other partitions
            continue;
          if(outPIndex + 1 < myPIndex || myPIndex + 1 < outPIndex)
            // We are actually at least two partitions from the one whose
            // overlap we are marking -- abandon
            DUNE_THROW(OverlapError, "Subpartition's overlap would cover " <<
                       "non-adjacent sibling.");
          auto &mark = overlapMap_[mapper_.map(outside)];
          if(mark == nooverlap)
            mark = myPIndex;
          else if(mark == myPIndex)
            /* nothing to do */;
          else  // conflict
            DUNE_THROW(OverlapError, "Subpartition's overlap would intersect "
                       << "other overlap");
        }
      };

      void checkOverlap()
      {
        if(overlap_ == 0)
          return;
        try {
          for(const auto &e : entitySet_)
          {
            auto myPIndex = pIndex(e);
            // check that we are either in the first or the last subpartition,
            // or outside the coarse overlap.  This avoids cases where the
            // first or last subpartition ends up empty, but the double
            // overlap between coarse border and the next subpartition goes
            // undetected.
            if(myPIndex > 0 && myPIndex + 1 < newPartitionIds_.size() &&
               overlapMap_[mapper_.map(e)] == someoverlap)
            {
              DUNE_THROW(OverlapError, "Inner subpartition intersects coarse "
                         << "overlap");
            }
            markneighbors(e, myPIndex);
          }

          for(std::size_t o = 1; o < overlap_; ++o)
            for(const auto &e : entitySet_) {
              std::size_t mark =
                overlapMap_[mapper_.map(e)];
              if(mark < newPartitionIds_.size())
                markneighbors(e, mark);
            }
        }
        catch(const OverlapError&)
        {
          // reset any marks to nooverlap
          for(const auto &e : entitySet_) {
            auto &mark = overlapMap_[mapper_.map(e)];
            if(mark != someoverlap)
              mark = nooverlap;
          }
          // rethrow
          throw;
        }
      }

      void finalize()
      {
        for(const auto &e : entitySet_) {
          auto &mark = overlapMap_[mapper_.map(e)];
          // translate overlap marks
          if(mark != nooverlap)
            mark = someoverlap;
          // set new partition number
          // be very careful here: Setting a new partition number may remove
          // the entity from the partition we are currently iterating over, so
          // the entitySet must not let itself confuse by this.
          partitioning_.setPartition(e, newPartitionIds_[pIndex(e)]);
        }
      }

    public:
      PartitioningContext(EquidistantPartitioner &partitioner,
                          const EntitySet &entitySet,
                          const std::vector<std::size_t> &newPartitionIds) :
        gv_(partitioner.gv_), partitioning_(partitioner.partitioning_),
        mapper_(partitioner.mapper_), overlapMap_(partitioner.overlapMap_),
        direction_(partitioner.direction_), overlap_(partitioner.overlap_),
        entitySet_(entitySet), newPartitionIds_(newPartitionIds)
      {
        initMinMax();
        checkOverlap();
        finalize();
      }

    };
  public:
    EquidistantPartitioner(const GV &gv, Partitioning &partitioning,
                           const Mapper &mapper, std::size_t direction,
                           std::size_t overlap = 1) :
      gv_(gv), partitioning_(partitioning), mapper_(mapper),
      overlapMap_(mapper_.size(), std::size_t(nooverlap)),
      direction_(direction), overlap_(overlap)
    { }

    template<class EntitySet>
    void partition(const EntitySet &entitySet,
                   const std::vector<std::size_t> &newPartitionIds)
    {
      // Constructor of temporary does all the work
      PartitioningContext<EntitySet>(*this, entitySet, newPartitionIds);
    }
  };

  template<class GV, class Partitioning>
  class RecursiveEquidistantPartitioner {
    typedef MultipleCodimMultipleGeomTypeMapper<GV, MCMGElementLayout> Mapper;
    typedef EquidistantPartitioner<GV, Partitioning, Mapper> DirPartitioner;

    const GV &gv_;
    Partitioning &partitioning_;
    Mapper mapper_;
    std::vector<DirPartitioner> dirPartitioners_;
    std::size_t partitions_;
    std::vector<std::size_t> lastDir_;
    std::vector<std::size_t> failCount_;
    std::vector<std::size_t> color_;

  public:
    RecursiveEquidistantPartitioner(const GV &gv, Partitioning &partitioning,
                                    std::size_t overlap = 1) :
      gv_(gv), partitioning_(partitioning), mapper_(gv_), partitions_(1),
      lastDir_(1, GV::dimension-1), failCount_(1, 0), color_(1, 0)
    {
      dirPartitioners_.reserve(GV::dimensionworld);
      for(std::size_t d = 0; d < GV::dimensionworld; ++d)
        dirPartitioners_.emplace_back(gv_, partitioning_, mapper_, d, overlap);
    }

    const Partitioning &partitioning() const
    {
      return partitioning_;
    }

    std::size_t color(std::size_t partition) const
    {
      return color_[partition];
    }

    std::vector<std::size_t> tryRefine(std::size_t pId, std::size_t dir,
                                       std::size_t divisor)
    {
      std::vector<std::size_t> newPartitionIds;
      newPartitionIds.reserve(divisor);
      newPartitionIds.push_back(pId);
      for(std::size_t p = partitions_; newPartitionIds.size() < divisor; ++p)
        newPartitionIds.push_back(p);
      dirPartitioners_[dir].partition(partitioning_.entitySet(pId),
                                      newPartitionIds);
      partitions_ += divisor-1;
      lastDir_.resize(partitions_, dir);
      lastDir_[pId] = dir;
      failCount_.resize(partitions_, 0);
      failCount_[pId] = 0;
      color_.resize(partitions_);
      for(std::size_t i = 1; i < divisor; ++i)
        color_[newPartitionIds[i]] = color_[pId] ^ ((i & 1) << dir);

      return std::move(newPartitionIds);
    }

    bool globalRefine() {
      bool result = false;
      std::size_t oldPartitions = partitions_;
      for(std::size_t p = 0; p < oldPartitions; ++p)
        if(failCount_[p] < GV::dimensionworld)
          try {
            tryRefine(p,
                      (lastDir_[p] + 1 + failCount_[p]) % GV::dimensionworld,
                      3);
            result = true;
          }
          catch(const OverlapError &) {
            failCount_[p]++;
          }
      return result;
    }
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_TRIPARTIT_HH
