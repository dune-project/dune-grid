// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_TRIPARTIT_HH
#define DUNE_GRID_UTILITY_TRIPARTIT_HH

#include <cstddef>
#include <vector>

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

  template<class GV, class SeedPartitioning, class MapPartitioning>
  class RecursiveEquidistantPartitioner {
    typedef MultipleCodimMultipleGeomTypeMapper<GV, MCMGElementLayout> Mapper;
    typedef OverlappedEquidistantPartitioner<typename GV::ctype>
      DirPartitioner;

    const GV &gv_;
    Mapper mapper_;
    SeedPartitioning &seedPartitioning_;
    MapPartitioning &mapPartitioning_;
    std::vector<OverlapMap<Mapper> > overlapMaps_;
    std::vector<std::size_t> lastDir_;
    std::vector<std::size_t> failCount_;
    std::vector<std::size_t> color_;
    std::size_t overlap_;

  public:
    RecursiveEquidistantPartitioner(const GV &gv,
                                    SeedPartitioning &seedPartitioning,
                                    MapPartitioning &mapPartitioning,
                                    std::size_t overlap = 1) :
      gv_(gv), mapper_(gv_), seedPartitioning_(seedPartitioning),
      mapPartitioning_(mapPartitioning), lastDir_(1, GV::dimension-1),
      failCount_(1, 0), color_(1, 0), overlap_(overlap)
    {
      overlapMaps_.reserve(GV::dimensionworld);
      for(std::size_t d = 0; d < GV::dimensionworld; ++d)
        overlapMaps_.emplace_back(mapper_);
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
      // construct iterable entity set
      typedef HybridEntitySet<typename MapPartitioning::Partition,
                              typename SeedPartitioning::Partition> EntitySet;
      EntitySet es(mapPartitioning_.partition(pId),
                   seedPartitioning_.partition(pId));

      // partition in the desired direction
      DirPartitioner dirPartitioner(gv_, es, dir, divisor, overlap_);

      // update seed partitioning
      auto newPartitionIds =
        seedPartitioning_.splitPartition(pId, dirPartitioner);

      // update mapped partitioning
      for(auto p : newPartitionIds)
        for(const auto &e : seedPartitioning_.partition(p))
          mapPartitioning_.setPartitionId(e, p);

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
          catch(const InvalidPartitionerError &) {
            failCount_[p]++;
          }
      return result;
    }
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_TRIPARTIT_HH
