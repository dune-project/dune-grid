// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_PARTITIONER_RECURSIVE_EQUIDISTANT_HH
#define DUNE_GRID_UTILITY_PARTITIONER_RECURSIVE_EQUIDISTANT_HH

#include <cstddef>
#include <deque>
#include <queue>
#include <utility>
#include <vector>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/utility/iterableentityset.hh>
#include <dune/grid/utility/partitioner/equidistant.hh>

namespace Dune {

  //! \brief Partitioner that can be used to recursively refine an equidistant
  //!        partitioning
  /**
   * This partitioner holds all the necessary data to generate a recursive
   * equidistant partitioning.  Users of the class have successively specify a
   * partition to refine, and how to refine it.
   */
  template<class GV, class SeedPartitioning_, class MapPartitioning_>
  class RecursiveEquidistantPartitioner {
    typedef MultipleCodimMultipleGeomTypeMapper<GV, MCMGElementLayout> Mapper;
    typedef OverlappedEquidistantPartitioner<typename GV::ctype>
      DirPartitioner;

  public:
    //! type of seed-list partitioning
    typedef SeedPartitioning_ SeedPartitioning;
    //! type of mapped partitioning
    typedef MapPartitioning_ MapPartitioning;

    //! contruct
    /**
     * \param gv               The grid view; used to obtain mappers and
     *                         intersections.
     * \param seedPartitioning Seed-list partitioning to fill.
     * \param mapPartitioning  Mapped partitioning to fill.
     * \param overlapSize      Size of overlap region.
     */
    RecursiveEquidistantPartitioner(const GV &gv,
                                    SeedPartitioning &seedPartitioning,
                                    MapPartitioning &mapPartitioning,
                                    std::size_t overlapSize = 1) :
      gv_(gv), mapper_(gv_), seedPartitioning_(seedPartitioning),
      mapPartitioning_(mapPartitioning), overlapSize_(overlapSize)
    {
      overlapMaps_.reserve(GV::dimensionworld);
      for(std::size_t d = 0; d < GV::dimensionworld; ++d)
        overlapMaps_.emplace_back(mapper_);
    }

    //! return the seed partitioning
    const SeedPartitioning seedPartitioning() const
    {
      return seedPartitioning_;
    }

    //! return the mapped partitioning
    const MapPartitioning mapPartitioning() const
    {
      return mapPartitioning_;
    }

    //! split a partition
    /**
     * \param pId           Id of partitions to split.
     * \param directions    Direction in wich to split.
     * \param subPartitions Number of sub-partitions to generate.
     *
     * \returns List of new sub-partition ids.  The first entry will be equal
     *          to \c pId, the remaining entries are allocated sequentially
     *          directly after the previously existing partition ids.
     *
     * \throw InvalidPartitionerError If the generated partitioning leads to
     *                                double overlap or would otherwise be
     *                                invalid.
     */
    std::vector<typename SeedPartitioning::Size>
    splitPartition(std::size_t pId, std::size_t direction,
                   std::size_t subPartitions)
    {
      // construct iterable entity set
      typedef HybridEntitySet<typename MapPartitioning::Partition,
                              typename SeedPartitioning::Partition> EntitySet;
      EntitySet es(mapPartitioning_.partition(pId),
                   seedPartitioning_.partition(pId));

      // partition in the desired direction
      DirPartitioner dirPartitioner(gv_, es, direction, subPartitions,
                                    overlapSize_);

      // update seed partitioning
      auto newPartitionIds =
        seedPartitioning_.splitPartition(pId, dirPartitioner);

      // update mapped partitioning
      for(auto p : newPartitionIds)
        for(const auto &e : seedPartitioning_.partition(p))
          mapPartitioning_.setPartitionId(e, p);

      return std::move(newPartitionIds);
    }

  private:
    const GV &gv_;
    Mapper mapper_;
    SeedPartitioning &seedPartitioning_;
    MapPartitioning &mapPartitioning_;
    std::vector<OverlapMap<Mapper> > overlapMaps_;
    std::size_t overlapSize_;
  };

  //! recursively refine until a certain number of partitions is reached
  template<class RecursivePartitioner>
  class RecursiveIsotropicRefiner
  {
    static const std::size_t dim = RecursivePartitioner::SeedPartitioning::
      Partition::Iterator::value_type::dimensionworld;

  public:
    //! type used to index partitions
    typedef typename RecursivePartitioner::SeedPartitioning::Size
      PartitionIndex;

    //! construct
    /**
     * \param partitioner   Partititioner to apply the refinement to.  This
     *                      can currently be a RecursiveEquidistantPartitioner
     *                      or a RecursiveOddPartitioner.
     * \param subPartitions Number of sub-partitions to split into each time a
     *                      partition is split.
     */
    RecursiveIsotropicRefiner(RecursivePartitioner &partitioner,
                              std::size_t subPartitions) :
      partitioner_(partitioner), subPartitions_(subPartitions),
      direction_(partitioner_.seedPartitioning().partitions(), dim-1),
      triedDirections_(partitioner_.seedPartitioning().partitions(), 0)
    {
      for(PartitionIndex p = 0;
          p < partitioner_.seedPartitioning().partitions();
          ++p)
      {
        eligiblePartitions_.emplace
          (partitioner_.seedPartitioning().partition(p).size(), p);
      }
    }

    //! try to refine the largest eligible partition
    /**
     * The initial try will be in the direction following the one in which the
     * candidates parent partition was partitioned.  Any further attempts will
     * try the next direction.  When all directions have been exhausted, the
     * partition looses eligible status.
     *
     * \return \c true if eligible partitions are left after the try, \c false
     *         otherwise.
     */
    bool tryRefine() {
      if(eligiblePartitions_.empty())
        return false;
      auto p = eligiblePartitions_.top().second;
      auto dir = (direction_[p] + triedDirections_[p] + 1) % dim;
      try {
        auto newPartitions =
          partitioner_.splitPartition(p, dir, subPartitions_);
        direction_.resize(direction_.size() + subPartitions_ - 1, dir);
        direction_[p] = dir;
        triedDirections_.resize(direction_.size(), 0);
        triedDirections_[p] = 0;
        eligiblePartitions_.pop();
        for(auto pp : newPartitions)
          eligiblePartitions_.emplace
            (partitioner_.seedPartitioning().partition(pp).size(), pp);
      }
      catch(const InvalidPartitionerError &) {
        // increment number of tries, remove from eligible partitions when all
        // directions were tried
        if(++triedDirections_[p] == dim)
          eligiblePartitions_.pop();
      }
      return !eligiblePartitions_.empty();
    }

    //! try to refine until there are \c targetPartitions partitions
    /**
     * This repeatedly calls tryRefine until there are either no eligible
     * partitions left or the number of partitions has reached \c
     * targetPartitions.
     *
     * \returns Whether there are eligible partitions left.
     */
    bool refine(PartitionIndex targetPartitions)
    {
      direction_.reserve(targetPartitions + subPartitions_ - 1);
      triedDirections_.reserve(targetPartitions + subPartitions_ - 1);
      while(partitioner_.seedPartitioning().partitions() < targetPartitions &&
            tryRefine())
      { }
      return !eligiblePartitions_.empty();
    }

  private:
    RecursivePartitioner &partitioner_;
    // how many partitions to try to partition into on each try
    std::size_t subPartitions_;
    // the direction of the partitioning this partition resulted from
    std::vector<std::size_t> direction_;
    // stored the count of directions we have tried
    std::vector<std::size_t> triedDirections_;
    // queue of the largest partitions, the size is first, the partition
    // number second (so we can reuse the default comparison)
    typedef std::pair<std::size_t, PartitionIndex> QValue;
    std::priority_queue<QValue, std::deque<QValue> > eligiblePartitions_;
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_PARTITIONER_RECURSIVE_EQUIDISTANT_HH
