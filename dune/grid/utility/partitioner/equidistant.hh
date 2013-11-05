// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_PARTITIONER_EQUIDISTANT_HH
#define DUNE_GRID_UTILITY_PARTITIONER_EQUIDISTANT_HH

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace Dune {

  //! A simple equidistant partitioner
  /**
   * This simply slices the domain into roughly equally spaced slices along
   * one coordinate axis, based solely on the center coordinates of the
   * codim<0> entities.
   *
   * \note The partitioner may yield some empty partitions.
   */
  template<class ctype>
  class EquidistantPartitioner
  {
  public:
    //! construct
    /**
     * The entityrange is used to figure out minimum and maximum coordinate.
     */
    template<class EntityRange>
    EquidistantPartitioner(const EntityRange &er, std::size_t direction,
                           std::size_t partitions) :
      minc_(std::numeric_limits<ctype>::infinity()), maxc_(-minc_),
      direction_(direction), partitions_(partitions)
    {
      using std::min;
      using std::max;
      for(const auto &e : er)
      {
        auto c = e.geometry().center()[direction_];
        minc_ = min(minc_, c);
        maxc_ = max(maxc_, c);
      }
    }

    //! Get partition number of an entity
    template<class Entity>
    std::size_t partition(const Entity &e) const
    {
      using std::min;
      using std::max;
      using std::floor;
      auto tmp = e.geometry().center()[direction_];
      tmp -= minc_;
      tmp /= maxc_ - minc_;
      tmp *= partitions_;
      tmp = floor(tmp);
      tmp = max(ctype(0), tmp);
      return min(partitions_-1, std::size_t(tmp));
    }

    //! get overall number of partitions
    std::size_t partitions() const
    {
      return partitions_;
    }

  private:
    ctype minc_;
    ctype maxc_;
    std::size_t direction_;
    std::size_t partitions_;
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_PARTITIONER_EQUIDISTANT_HH
