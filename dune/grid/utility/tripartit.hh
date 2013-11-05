// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_TRIPARTIT_HH
#define DUNE_GRID_UTILITY_TRIPARTIT_HH

#include <cstddef>
#include <limits>
#include <list>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/utility/iterableentityset.hh>
#include <dune/grid/utility/partitioner/equidistant.hh>
#include <dune/grid/utility/partitioning/seedlist.hh>

namespace Dune {

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

  class RepartitionPolicyInterface {
  public:
    void setNumSubpartitions(std::size_t npartitions);

    template<class Entity>
    void addToSubPartition(std::size_t subpartition, const Entity &e);
  };

  template<class SeedPartitioning, class MapPartitioning>
  class HybridRepartitionPolicy {
    typedef typename SeedPartitioning::Partition::Entity::EntitySeed Seed;

    SeedPartitioning &seedPartitioning_;
    MapPartitioning &mapPartitioning_;
    std::size_t oldPartition_;
    std::vector<typename SeedPartitioning::Size> &newPartitions_;
    std::vector<std::list<Seed> > data_;

  public:
    HybridRepartitionPolicy(SeedPartitioning &seedPartitioning,
                            MapPartitioning &mapPartitioning,
                            std::size_t oldPartition,
                            std::vector<typename SeedPartitioning::Size> &newPartitions) :
      seedPartitioning_(seedPartitioning), mapPartitioning_(mapPartitioning),
      oldPartition_(oldPartition), newPartitions_(newPartitions)
    { }

    void setNumSubPartitions(std::size_t npartitions)
    {
      data_.resize(npartitions);
    }

    template<class Entity>
    void addToSubPartition(std::size_t subpartition, const Entity &e)
    {
      data_[subpartition].push_back(e.seed());
    }

    void commit() {
      newPartitions_ = seedPartitioning_.splitPartition(oldPartition_, data_);
      for(std::size_t partition : newPartitions_)
        for(const auto &e : seedPartitioning_.partition(partition))
          mapPartitioning_.setPartitionId(e, partition);
    }
  };

  struct OverlapError : Exception {};

  template<class GV, class Mapper>
  class OverlappedEquidistantPartitioner
  {
    static const std::size_t nooverlap =
      std::numeric_limits<std::size_t>::max();
    static const std::size_t someoverlap = nooverlap - 1;

    const GV &gv_;
    const Mapper &mapper_;
    std::vector<std::size_t> overlapMap_;
    std::size_t direction_;
    std::size_t overlap_;

    template<class EntitySet, class RepartitionPolicy>
    class PartitioningContext {
      typedef typename EntitySet::Entity Entity;
      typedef typename Entity::ctype ctype;

      const GV &gv_;
      const Mapper &mapper_;
      const EntitySet &entitySet_;
      RepartitionPolicy repartitionPolicy_;
      std::vector<std::size_t> &overlapMap_;
      std::size_t overlap_;
      EquidistantPartitioner<ctype> basePartitioner_;

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
          std::size_t outPIndex = basePartitioner_.partition(outside);
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

      void checkNonempty(const std::vector<std::size_t> &pSizes)
      {
        for(auto size : pSizes)
          if(size == 0)
            DUNE_THROW(OverlapError,
                       "Subpartioning results in empty partition");
      }

      void checkOverlap()
      {
        std::vector<std::size_t> pSizes(basePartitioner_.partitions(), 0);
        if(overlap_ == 0)
        {
          for(const auto &e : entitySet_)
            ++pSizes[basePartitioner_.partition(e)];
          checkNonempty(pSizes);
          return;
        }
        try {
          for(const auto &e : entitySet_)
          {
            auto myPIndex = basePartitioner_.partition(e);
            ++pSizes[myPIndex];
            // check that we are either in the first or the last subpartition,
            // or outside the coarse overlap.  This avoids cases where the
            // first or last subpartition ends up empty, but the double
            // overlap between coarse border and the next subpartition goes
            // undetected.
            if(myPIndex > 0 && myPIndex + 1 < basePartitioner_.partitions() &&
               overlapMap_[mapper_.map(e)] == someoverlap)
            {
              DUNE_THROW(OverlapError, "Inner subpartition intersects coarse "
                         << "overlap");
            }
            markneighbors(e, myPIndex);
          }
          checkNonempty(pSizes);

          for(std::size_t o = 1; o < overlap_; ++o)
            for(const auto &e : entitySet_) {
              std::size_t mark =
                overlapMap_[mapper_.map(e)];
              if(mark < basePartitioner_.partitions())
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
        repartitionPolicy_.setNumSubPartitions(basePartitioner_.partitions());
        for(const auto &e : entitySet_) {
          auto &mark = overlapMap_[mapper_.map(e)];
          // translate overlap marks
          if(mark != nooverlap)
            mark = someoverlap;
          // set new partition number
          // be very careful here: Setting a new partition number may remove
          // the entity from the partition we are currently iterating over, so
          // the entitySet must not let itself confuse by this.
          repartitionPolicy_.addToSubPartition(basePartitioner_.partition(e),
                                               e);
        }
        repartitionPolicy_.commit();
      }

    public:
      PartitioningContext(OverlappedEquidistantPartitioner &partitioner,
                          const EntitySet &entitySet,
                          const RepartitionPolicy &repartitionPolicy,
                          std::size_t divisor) :
        gv_(partitioner.gv_), mapper_(partitioner.mapper_),
        entitySet_(entitySet), repartitionPolicy_(repartitionPolicy),
        overlapMap_(partitioner.overlapMap_), overlap_(partitioner.overlap_),
        basePartitioner_(entitySet, partitioner.direction_, divisor)
      {
        checkOverlap();
        finalize();
      }

    };
  public:
    OverlappedEquidistantPartitioner(const GV &gv, const Mapper &mapper,
                                     std::size_t direction,
                                     std::size_t overlap = 1) :
      gv_(gv), mapper_(mapper),
      overlapMap_(mapper_.size(), std::size_t(nooverlap)),
      direction_(direction), overlap_(overlap)
    { }

    template<class EntitySet, class RepartitionPolicy>
    void partition(const EntitySet &entitySet,
                   const RepartitionPolicy &repartitionPolicy,
                   std::size_t divisor)
    {
      // Constructor of temporary does all the work
      PartitioningContext<EntitySet, RepartitionPolicy>(*this, entitySet,
                                                        repartitionPolicy,
                                                        divisor);
    }
  };

  template<class GV, class SeedPartitioning, class MapPartitioning>
  class RecursiveEquidistantPartitioner {
    typedef MultipleCodimMultipleGeomTypeMapper<GV, MCMGElementLayout> Mapper;
    typedef OverlappedEquidistantPartitioner<GV, Mapper> DirPartitioner;

    const GV &gv_;
    Mapper mapper_;
    SeedPartitioning &seedPartitioning_;
    MapPartitioning &mapPartitioning_;
    std::vector<DirPartitioner> dirPartitioners_;
    std::vector<std::size_t> lastDir_;
    std::vector<std::size_t> failCount_;
    std::vector<std::size_t> color_;

  public:
    RecursiveEquidistantPartitioner(const GV &gv,
                                    SeedPartitioning &seedPartitioning,
                                    MapPartitioning &mapPartitioning,
                                    std::size_t overlap = 1) :
      gv_(gv), mapper_(gv_), seedPartitioning_(seedPartitioning),
      mapPartitioning_(mapPartitioning), lastDir_(1, GV::dimension-1),
      failCount_(1, 0), color_(1, 0)
    {
      dirPartitioners_.reserve(GV::dimensionworld);
      for(std::size_t d = 0; d < GV::dimensionworld; ++d)
        dirPartitioners_.emplace_back(gv_, mapper_, d, overlap);
    }

    std::size_t color(std::size_t partition) const
    {
      return color_[partition];
    }

    std::size_t color(const typename GV::template Codim<0>::Entity &e) const
    {
      return color_[mapPartitioning_.getPartitionId(e)];
    }

    std::vector<typename SeedPartitioning::Size>
    tryRefine(std::size_t pId, std::size_t dir, std::size_t divisor)
    {
      typedef HybridEntitySet<typename MapPartitioning::Partition,
                              typename SeedPartitioning::Partition> EntitySet;
      typedef HybridRepartitionPolicy<SeedPartitioning,
                                      MapPartitioning> RepartitionPolicy;
      std::vector<typename SeedPartitioning::Size> newPartitionIds;
      dirPartitioners_[dir].partition
        (EntitySet(mapPartitioning_.partition(pId),
                   seedPartitioning_.partition(pId)),
         RepartitionPolicy(seedPartitioning_, mapPartitioning_, pId,
                           newPartitionIds),
         divisor);
      lastDir_.resize(seedPartitioning_.partitions(), dir);
      lastDir_[pId] = dir;
      failCount_.resize(seedPartitioning_.partitions(), 0);
      failCount_[pId] = 0;
      color_.resize(seedPartitioning_.partitions());
      for(std::size_t i = 1; i < divisor; ++i)
        color_[newPartitionIds[i]] = color_[pId] ^ ((i & 1) << dir);

      return std::move(newPartitionIds);
    }

    bool globalRefine() {
      bool result = false;
      // new partitions are added at the end; remember the old limit so we
      // don't try to refine partitions we just added.
      std::size_t oldPartitions = seedPartitioning_.partitions();
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
