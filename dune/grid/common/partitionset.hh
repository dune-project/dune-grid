// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_PARTITIONSET_HH
#define DUNE_GRID_COMMON_PARTITIONSET_HH

#include <dune/common/typetraits.hh>
#include <dune/grid/common/gridenums.hh>

namespace Dune {

  /**
   * \addtogroup gridpartitions Parallel Grid Partitions
   * \{
   */

  namespace {

    // Simple TMP to deduce partition iterator type from set of partitions.
    template<unsigned int partitions>
    struct derive_partition_iterator_type
    {
      // We did not match any specialization, bail out...
      static_assert(AlwaysFalse<std::integral_constant<unsigned int,partitions> >::value,
                    "There is no partition iterator for this combination of entity partitions");
    };


    // specializations of derive_partition_iterator_type for existing PartitionIteratorTypes

    template<>
    struct derive_partition_iterator_type<
      (1 << InteriorEntity)
        >
        : public std::integral_constant<PartitionIteratorType,Interior_Partition>
    {};

    template<>
    struct derive_partition_iterator_type<
      (1 << InteriorEntity) |
      (1 << BorderEntity)
      >
      : public std::integral_constant<PartitionIteratorType,InteriorBorder_Partition>
    {};

    template<>
    struct derive_partition_iterator_type<
      (1 << InteriorEntity) |
      (1 << BorderEntity) |
      (1 << OverlapEntity)
      >
      : public std::integral_constant<PartitionIteratorType,Overlap_Partition>
    {};

    template<>
    struct derive_partition_iterator_type<
      (1 << InteriorEntity) |
      (1 << BorderEntity) |
      (1 << OverlapEntity) |
      (1 << FrontEntity)
      >
      : public std::integral_constant<PartitionIteratorType,OverlapFront_Partition>
    {};

    template<>
    struct derive_partition_iterator_type<
      (1 << InteriorEntity) |
      (1 << BorderEntity) |
      (1 << OverlapEntity) |
      (1 << FrontEntity) |
      (1 << GhostEntity)
      >
      : public std::integral_constant<PartitionIteratorType,All_Partition>
    {};

    template<>
    struct derive_partition_iterator_type<
      (1 << GhostEntity)
        >
        : public std::integral_constant<PartitionIteratorType,Ghost_Partition>
    {};

  } // anonymous namespace


  //! A set of PartitionType values.
  /**
   * PartitionSet cotains a set of PartitionType values fixed at compile time. The contents
   * of the set is encoded in the template parameter partitions, but the exact semantics are
   * an implementation detail. PartitionSets can be combined by adding them up. They also support
   * removing partitions by subtracting them.
   *
   * \tparam partitions  Implementation-defined representation of the partition set.
   */
  template<unsigned int partitions>
  struct PartitionSet
  {
    //! \private The actual representation should not be used outside of the implementation.
    static const unsigned int value = partitions;


    //! Returns a new PartitionSet that also contains the partitions in set.
    template<unsigned int p>
    struct PartitionSet<partitions | p>
    operator+(const PartitionSet<p>& set)
    {
      return PartitionSet<partitions | p>();
    }

    //! Returns a new PartitionSet that does not contain the partitions in set.
    template<unsigned int p>
    struct PartitionSet<partitions & ~p>
    operator-(const PartitionSet<p>& set)
    {
      return PartitionSet<partitions & ~p>();
    }

    //! Writes the PartitionSet to an output stream.
    template<unsigned int partitions>
    friend std::ostream& operator<<(std::ostream& os, const PartitionSet<partitions>&)
    {
      unsigned int set = partitions;
      os << "partition set {";
      bool first = true;
      for (unsigned int p = 0; set; set &= ~(1 << p++))
        {
          if (!(set & (1 << p)))
            continue;
          if (!first)
            os << ",";
          first = false;
          os << static_cast<PartitionType>(p);
        }
      os << "}";
      return os;
    }

#if HAVE_CONSTEXPR || DOXYGEN

    //! Returns the PartitionIteratorType that can be used to iterate over the partitions in the set.
    /**
     *
     * \throws  raises a static assertion if the partitions do not correspond to a valid PartitionIteratorType.
     * \since   GCC 4.6 (constexpr)
     */
    static constexpr PartitionIteratorType partitionIterator()
    {
      return derive_partition_iterator_type<partitions>::value;
    }

#endif // HAVE_CONSTEXPR

  };

  //! Creates a PartitionSet for the given PartitionType.
  /**
   * \related PartitionSet
   */
  template<PartitionType p>
  PartitionSet<(1 << p)> partitionSet()
  {
    return PartitionSet<(1 << p)>();
  }

  //! Predefined PartitionSets for commonly used combinations of parallel grid PartitionTypes.
  namespace Partitions {


#ifdef DOXYGEN

    //! PartitionSet for the interior partition.
    PartitionSet<...> interior;

    //! PartitionSet for the border partition.
    PartitionSet<...> border;

    //! PartitionSet for the overlap partition.
    PartitionSet<...> overlap;

    //! PartitionSet for the front partition.
    PartitionSet<...> front;

    //! PartitionSet for the ghost partition.
    PartitionSet<...> ghost;

    //! PartitionSet for the interior and border partitions.
    PartitionSet<...> interiorBorder;

    //! PartitionSet for the interior, border and overlap partitions.
    PartitionSet<...> interiorBorderOverlap;

    //! PartitionSet for the interior, border, overlap and front partitions.
    PartitionSet<...> interiorBorderOverlapFront;

    //! PartitionSet for all partitions.
    PartitionSet<...> all;

#else // DOXYGEN

    namespace {

      // place global objects in anonymous namespace to ensure that visibility is
      // restricted to the current translation unit, making it easier for the compiler
      // to eliminate the actual objects and to avoid linking problems

      const decltype(partitionSet<InteriorEntity>()) interior = {};
      const decltype(partitionSet<BorderEntity>()) border = {};
      const decltype(partitionSet<OverlapEntity>()) overlap = {};
      const decltype(partitionSet<FrontEntity>()) front = {};
      const decltype(partitionSet<GhostEntity>()) ghost = {};

      const decltype(interior + border) interiorBorder = {};
      const decltype(interior + border + overlap) interiorBorderOverlap = {};
      const decltype(interior + border + overlap + front) interiorBorderOverlapFront = {};
      const decltype(interior + border + overlap + front + ghost) all = {};

    }

#endif // DOXYGEN

  } // namespace Partitions

  /**
   * \} group Parallel Grid Partitions
   */

} // namespace Dune

#endif // DUNE_GRID_COMMON_PARTITIONSET_HH
