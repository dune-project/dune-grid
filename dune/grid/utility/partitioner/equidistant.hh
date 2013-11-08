// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_PARTITIONER_EQUIDISTANT_HH
#define DUNE_GRID_UTILITY_PARTITIONER_EQUIDISTANT_HH

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/utility/intersectionrange.hh>

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

    //! get the color of a partition
    std::size_t color(std::size_t partition) const
    {
      return partition % 2;
    }

    //! get the number of colors in the partitioner
    /**
     * \returns 1 if there is just one partition, 2 otherwise.
     */
    std::size_t colors() const
    {
      return 1 + (partitions() > 1);
    }

  private:
    ctype minc_;
    ctype maxc_;
    std::size_t direction_;
    std::size_t partitions_;
  };

  //! Hold overlap information
  /**
   * This simply derives from std::vector and adds an overload to operator[]
   * to allow direct subscription by entities.
   *
   * There are two member constants \c nooverlap and \c someoverlap, which are
   * chosed such that they occupy the largest positive values.  The convention
   * is that the overlap map can store either \c nooverlap or \c someoverlap,
   * or some partition number.  \see OverlappedEquidistantPartitioner.
   */
  template<class Mapper>
  class OverlapMap :
    public std::vector<std::size_t>
  {
    typedef std::vector<std::size_t> Base;

  public:
    static const std::size_t nooverlap =
      std::numeric_limits<std::size_t>::max();
    static const std::size_t someoverlap = nooverlap - 1;

    OverlapMap(const Mapper &mapper) :
      Base(mapper.size(), nooverlap),
      mapper_(mapper)
    { }

    using Base::operator[];

    template<class Entity>
    typename std::enable_if<!std::is_integral<Entity>::value,
                            std::size_t &>::type
    operator[](const Entity &e)
    {
      return (*this)[mapper_.map(e)];
    }
    template<class Entity>
    typename std::enable_if<!std::is_integral<Entity>::value,
                            const std::size_t &>::type
    operator[](const Entity &e) const
    {
      return (*this)[mapper_.map(e)];
    }

  private:
    Mapper mapper_;
  };

  // define the constants here, otherwise using them e.g. as the
  // initialization value of a std::vector does not work, because that
  // requires binding a const reference to them
  template<class Mapper>
  const std::size_t OverlapMap<Mapper>::nooverlap;
  template<class Mapper>
  const std::size_t OverlapMap<Mapper>::someoverlap;

  //! \brief exception thrown from the contructor of partitioners when a valid
  //!        partitioning for the requested parameters cannot be obtained
  struct InvalidPartitionerError : Exception {};

  template<class ctype>
  class OverlappedEquidistantPartitioner :
    public EquidistantPartitioner<ctype>
  {
    typedef EquidistantPartitioner<ctype> Base;

    //markneighbors
    template<class GV, class EntitySet, class OverlapMap>
    void markneighbors(const GV &gv, const EntitySet &es,
                       OverlapMap &overlapMap,
                       const typename GV::template Codim<0>::Entity &e,
                       std::size_t myPIndex) const
    {
      for(const auto &is : intersections(gv, e))
      {
        if(!is.neighbor())
          continue;
        auto outsidep = is.outside();
        const auto &outside = *outsidep;
        if(!es.contains(outside))
          // These were marked previously in a coarser partitioning
          continue;
        // Now it is legal to obtain the local index for this partition
        std::size_t outPIndex = this->partition(outside);
        if(outPIndex == myPIndex)
          // Only mark our overlap in other partitions
          continue;
        if(outPIndex + 1 < myPIndex || myPIndex + 1 < outPIndex)
          // We are actually at least two partitions from the one whose
          // overlap we are marking -- abandon
          DUNE_THROW(InvalidPartitionerError, "Subpartition's overlap would "
                     "cover non-adjacent sibling.");
        auto &mark = overlapMap[outside];
        if(mark == OverlapMap::nooverlap)
          mark = myPIndex;
        else if(mark == myPIndex)
        { /* nothing to do */ }
        else  // conflict
          DUNE_THROW(InvalidPartitionerError, "Subpartition's overlap would "
                     "intersect other overlap");
      }
    };

    template<class GV, class IterableEntitySet>
    void checkOverlap(const GV &gv, const IterableEntitySet &es,
                      std::size_t overlapSize) const
    {
      typedef MultipleCodimMultipleGeomTypeMapper<GV, MCMGElementLayout>
        Mapper;

      OverlapMap<Mapper> overlapMap((Mapper(gv)));
      markOverlap(gv, es, overlapSize, overlapMap, true);
    }

    //! Mark overlap, and as a side effect ensure partitioning is valid
    /**
     * This marks overlap in the domain of \c es.  On entry, overlapMap should
     * contain either \c OverlapMap::nooverlap or \c OverlapMap::someoverlap
     * for each entity in \c es.  This function sets \c overlapMap[e] to some
     * partition number p if \c e is in the overlap of p.  This is done by
     * growing the overlap layer by layer.  If at any time \c overlapMap[e]
     * contains a value other than p or \c OverlapMap::nooverlap prior to
     * setting it to p, an \c InvalidPartitionerError is thrown, and all
     * partition numbers are reset to \c OverlapMap::nooverlap, restoring the
     * overlapMap to its initial state in the domain of es.  If no conflicts
     * are detected, the overlapMap is finalized by replacing all partition
     * numbers by \c OverlapMap::someoverlap.
     *
     * If \c checkOnly is set to true, the overlapMap is assumed to be a
     * temporary and it is neither finalized nor restored.
     */
    template<class GV, class IterableEntitySet, class OverlapMap>
    void markOverlap(const GV &gv, const IterableEntitySet &es,
                     std::size_t overlapSize, OverlapMap &overlapMap,
                     bool checkOnly = false) const
    {
      if(overlapSize == 0)
        // nothing to check
        return;
      try {
        std::vector<std::size_t> pSizes(this->partitions(), 0);
        for(const auto &e : es)
        {
          auto myPIndex = this->partition(e);
          ++pSizes[myPIndex];
          // check that we are either in the first or the last subpartition,
          // or outside the coarse overlap.  This avoids cases where the first
          // or last subpartition ends up empty, but the double overlap
          // between coarse border and the next subpartition goes undetected.
          if(myPIndex > 0 && myPIndex + 1 < this->partitions() &&
             overlapMap[e] == OverlapMap::someoverlap)
          {
            DUNE_THROW(InvalidPartitionerError, "Inner subpartition "
                       "intersects coarse overlap");
          }
          markneighbors(gv, es, overlapMap, e, myPIndex);
        }

        // check that partitions are nonempty, for empty partitions we may be
        // unable to detect double overlap.
        for(auto size : pSizes)
          if(size == 0)
            DUNE_THROW(InvalidPartitionerError,
                       "Subpartioning results in empty partition");

        // mark higher layers of overlap
        // note: this may currently mark a larger overlap than actually
        // desired.
        for(std::size_t o = 1; o < overlapSize; ++o)
          for(const auto &e : es) {
            std::size_t mark = overlapMap[e];
            if(mark < this->partitions())
              markneighbors(gv, es, overlapMap, e, mark);
          }

        // translate overlap marks
        if(!checkOnly)
          for(const auto &e : es) {
            auto &mark = overlapMap[e];
            if(mark != OverlapMap::nooverlap)
              mark = OverlapMap::someoverlap;
          }
      }
      catch(const InvalidPartitionerError &) {
        // reset any marks to nooverlap
        if(!checkOnly)
          for(const auto &e : es) {
            auto &mark = overlapMap[e];
            if(mark != OverlapMap::someoverlap)
              mark = OverlapMap::nooverlap;
          }
        throw;
      }
    }

  public:
    //! constructor
    /**
     * Try to partition into \c partitions partition in direction \c
     * direction, ensuring an overlap of \c overlapSize is observed.  This
     * uses \c overlapMap to mark the overlap; it may contain overlap marks
     * from a coarser partioning, which will be taken into account.  On exit,
     * \c overlapMap will be updated with overlap marks for the new
     * partitions.  If partitioning fails, overlapMap is restored to its old
     * state and an \c InvalidPartitionerError is thrown.
     *
     * \c gv is used to obtain intersection iterators and \c es defined the
     * domain in which to partition.  Overlap marks are never access (neither
     * read nor write) outside of this domain.
     */
    template<class GV, class IterableEntitySet, class OverlapMap>
    OverlappedEquidistantPartitioner(const GV &gv, const IterableEntitySet &es,
                                     std::size_t direction,
                                     std::size_t partitions,
                                     std::size_t overlapSize,
                                     OverlapMap &overlapMap) :
      Base(es, direction, partitions)
    {
      markOverlap(gv, es, overlapSize, overlapMap);
    }

    //! constructor
    /**
     * Just like the other constructor, but use a temporary overlapMap
     */
    template<class GV, class IterableEntitySet>
    OverlappedEquidistantPartitioner(const GV &gv, const IterableEntitySet &es,
                                     std::size_t direction,
                                     std::size_t partitions,
                                     std::size_t overlapSize) :
      Base(es, direction, partitions)
    {
      checkOverlap(gv, es, overlapSize);
    }

  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_PARTITIONER_EQUIDISTANT_HH
