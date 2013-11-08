// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_PARTITIONER_TENSOR_EQUIDISTANT_HH
#define DUNE_GRID_UTILITY_PARTITIONER_TENSOR_EQUIDISTANT_HH

#include <cstddef>
#include <vector>

#include <dune/grid/utility/iterableentityset.hh>
#include <dune/grid/utility/partitioner/equidistant.hh>

namespace Dune {

  //! Apply equidistant partitioning seperately in each direction
  template<class ctype>
  class TensorEquidistantPartitioner {
    typedef OverlappedEquidistantPartitioner<ctype> DirPartitioner;

    template<class GV, class IterableEntitySet, class Vec>
    void init(const GV &gv, const IterableEntitySet &es,
              const Vec &subPartitions, std::size_t overlapSize)
    {
      dirPartitioners_.reserve(GV::dimensionworld);
      for(std::size_t d = 0; d < GV::dimensionworld; ++d)
        dirPartitioners_.emplace_back(gv, es, d, subPartitions[d],
                                      overlapSize);
    }

  public:
    //! contruct
    /**
     * \param gv            The grid view; used to obtain mappers and
     *                      intersections.
     * \param es            The entity set defining the domain.
     * \param subPartitions Vector of subpartitions for each direction.
     *                      Should be something subscriptable.
     * \param overlapSize   Size of overlap region.
     */
    template<class GV, class IterableEntitySet, class Vec>
    TensorEquidistantPartitioner(const GV &gv, const IterableEntitySet &es,
                                 const Vec &subPartitions,
                                 std::size_t overlapSize = 1)
    {
      init(gv, es, subPartitions, overlapSize);
    }

    //! contruct
    /**
     * \param gv            The grid view; used to obtain mappers and
     *                      intersections; also defined the domain.
     * \param subPartitions Vector of subpartitions for each direction.
     *                      Should be something subscriptable.
     * \param overlapSize   Size of overlap region.
     */
    template<class GV, class Vec>
    TensorEquidistantPartitioner(const GV &gv, const Vec &subPartitions,
                                 std::size_t overlapSize = 1)
    {
      init(gv, entities<0>(gv), subPartitions, overlapSize);
    }

    //! Get partition number of an entity
    template<class Entity>
    std::size_t partition(const Entity &e) const
    {
      std::size_t base = 1;
      std::size_t index = 0;
      for(const auto &p : dirPartitioners_)
      {
        index += base * p.partition(e);
        base *= p.partitions();
      }
      return index;
    }

    //! get overall number of partitions
    std::size_t partitions() const
    {
      std::size_t prod = 1;
      for(const auto &p : dirPartitioners_)
        prod *= p.partitions();
      return prod;
    }

    //! get the color of a partition
    std::size_t color(std::size_t partition) const
    {
      std::size_t base = 1;
      std::size_t color = 0;
      for(const auto &p : dirPartitioners_)
      {
        auto partitions = p.partitions();
        auto dPartition = partition % partitions;
        partition /= partitions;
        color += p.color(dPartition) * p.colors();
        base *= p.colors();
      }

      return color;
    }

    //! get the number of colors in the partitioner
    /**
     * \returns 1 if there is just one partition, 2 otherwise.
     */
    std::size_t colors() const
    {
      std::size_t prod = 1;
      for(const auto &p : dirPartitioners_)
        prod *= p.colors();
      return prod;
    }

  private:
    std::vector<DirPartitioner> dirPartitioners_;
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_PARTITIONER_TENSOR_EQUIDISTANT_HH
