// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_PARTITIONER_RECURSIVE_ODD_HH
#define DUNE_GRID_UTILITY_PARTITIONER_RECURSIVE_ODD_HH

#include <cstddef>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/grid/utility/partitioner/equidistant.hh>
#include <dune/grid/utility/partitioner/recursive-equidistant.hh>

namespace Dune {

  //! \brief \c RecursiveEquidistantPartitioner with coloring
  /**
   * This is basically the same as \c RecursiveEquidistantPartitioner, except
   * that partitions cann only be split into an odd number of sub-partitions.
   * This has the nice benefit that a coloring can be easily provided.
   */
  template<class GV, class SeedPartitioning, class MapPartitioning>
  class RecursiveOddPartitioner :
    public RecursiveEquidistantPartitioner<GV, SeedPartitioning,
                                           MapPartitioning>
  {
    typedef RecursiveEquidistantPartitioner<GV, SeedPartitioning,
                                            MapPartitioning> Base;
  public:
    //! contruct
    /**
     * \param gv               The grid view; used to obtain mappers and
     *                         intersections.
     * \param seedPartitioning Seed-list partitioning to fill.
     * \param mapPartitioning  Mapped partitioning to fill.
     * \param overlapSize      Size of overlap region.
     */
    RecursiveOddPartitioner(const GV &gv, SeedPartitioning &seedPartitioning,
                            MapPartitioning &mapPartitioning,
                            std::size_t overlapSize = 1) :
      Base(gv, seedPartitioning, mapPartitioning, overlapSize),
      colors_(1, 0)
    { }

    //! split a partition
    /**
     * \param pId           Id of partitions to split.
     * \param directions    Direction in wich to split.
     * \param subPartitions Number of sub-partitions to generate.  Must be odd.
     *
     * \returns List of new sub-partition ids.  The first entry will be equal
     *          to \c pId, the remaining entries are allocated sequentially
     *          directly after the previously existing partition ids.
     *
     * \throw InvalidPartitionerError If the generated partitioning leads to
     *                                double overlap or would otherwise be
     *                                invalid, or \c subPartitiond is even.
     */
    std::vector<typename SeedPartitioning::Size>
    splitPartition(std::size_t pId, std::size_t direction,
                   std::size_t subPartitions)
    {
      if(subPartitions % 2 != 1)
        DUNE_THROW(InvalidPartitionerError, "Attempt to split into " <<
                   subPartitions << " sub-partitions in "
                   "RecursiveOddPartitioner; number of sub-partitions must be "
                   "odd.");
      auto newPartitions = Base::splitPartition(pId, direction, subPartitions);
      auto color = colors_[pId];
      auto toggle = 1 << direction;
      colors_.resize(this->seedPartitioning().partitions());
      for(auto p : newPartitions) {
        colors_[p] = color;
        color ^= toggle;
      }
      return std::move(newPartitions);
    }

    //! get color of a certain partition
    /**
     * There are always 2^dimw colors.
     */
    std::size_t color(std::size_t pId) const
    {
      return colors_[pId];
    }

  private:
    // make sure type for color is unsigned, we're using bit-arithmetic
    std::vector<std::size_t> colors_;
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_PARTITIONER_RECURSIVE_ODD_HH
