// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_PARTITIONSET_HH
#define DUNE_GRID_COMMON_PARTITIONSET_HH

#include <dune/common/typetraits.hh>
#include <dune/grid/common/gridenums.hh>

namespace Dune {

  /**
   * \addtogroup GIRelatedTypes
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
    operator+(const PartitionSet<p>& set) const
    {
      return PartitionSet<partitions | p>();
    }

    //! Returns a new PartitionSet that does not contain the partitions in set.
    template<unsigned int p>
    struct PartitionSet<partitions & ~p>
    operator-(const PartitionSet<p>& set) const
    {
      return PartitionSet<partitions & ~p>();
    }

    //! Writes the PartitionSet to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const PartitionSet&)
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

    //! Returns the PartitionIteratorType that can be used to iterate over the partitions in the set.
    /**
     *
     * \throws  raises a static assertion if the partitions do not correspond to a valid PartitionIteratorType.
     */
    static constexpr PartitionIteratorType partitionIterator()
    {
      return derive_partition_iterator_type<partitions>::value;
    }

    //! Tests whether the given PartitionType is contained in this set.
    static constexpr bool contains(PartitionType pt)
    {
      return partitions & (1 << pt);
    }

    //! Tests whether the given PartitionSet is contained in this set.
    template<unsigned int contained_partitions>
    static constexpr bool contains(PartitionSet<contained_partitions>)
    {
      return (partitions & contained_partitions) == contained_partitions;
    }

    //! Tests whether two PartitionsSet are equal.
    template<unsigned int p2>
    constexpr bool operator==(PartitionSet<p2>) const
    {
      return partitions == p2;
    }

    //! Tests whether two PartitionsSet are not equal.
    template<unsigned int p2>
    constexpr bool operator!=(PartitionSet<p2>) const
    {
      return partitions != p2;
    }

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

    //! Type of PartitionSet for the interior partition.
    typedef PartitionSet<...> Interior;

    //! Type of PartitionSet for the border partition.
    typedef PartitionSet<...> Border;

    //! Type of PartitionSet for the overlap partition.
    typedef PartitionSet<...> Overlap;

    //! Type of PartitionSet for the front partition.
    typedef PartitionSet<...> Front;

    //! Type of PartitionSet for the ghost partition.
    typedef PartitionSet<...> Ghost;

    //! Type of PartitionSet for the interior and border partitions.
    typedef PartitionSet<...> InteriorBorder;

    //! Type of PartitionSet for the interior, border and overlap partitions.
    typedef PartitionSet<...> InteriorBorderOverlap;

    //! Type of PartitionSet for the interior, border, overlap and front partitions.
    typedef PartitionSet<...> InteriorBorderOverlapFront;

    //! Type of PartitionSet for all partitions.
    typedef PartitionSet<...> All;


    //! PartitionSet for the interior partition.
    Interior interior;

    //! PartitionSet for the border partition.
    Border border;

    //! PartitionSet for the overlap partition.
    Overlap overlap;

    //! PartitionSet for the front partition.
    Front front;

    //! PartitionSet for the ghost partition.
    Ghost ghost;

    //! PartitionSet for the interior and border partitions.
    InteriorBorder interiorBorder;

    //! PartitionSet for the interior, border and overlap partitions.
    InteriorBorderOverlap interiorBorderOverlap;

    //! PartitionSet for the interior, border, overlap and front partitions.
    InteriorBorderOverlapFront interiorBorderOverlapFront;

    //! PartitionSet for all partitions.
     All all;

#else // DOXYGEN

    // First declare the types and objects for individual partitions

    typedef decltype(partitionSet<InteriorEntity>()) Interior;
    typedef decltype(partitionSet<BorderEntity>()) Border;
    typedef decltype(partitionSet<OverlapEntity>()) Overlap;
    typedef decltype(partitionSet<FrontEntity>()) Front;
    typedef decltype(partitionSet<GhostEntity>()) Ghost;

    namespace {

      // place global objects in anonymous namespace to ensure that visibility is
      // restricted to the current translation unit, making it easier for the compiler
      // to eliminate the actual objects and to avoid linking problems

      const Interior interior = {};
      const Border border = {};
      const Overlap overlap = {};
      const Front front = {};
      const Ghost ghost = {};

    }

    // Now we can declare the partition sets that are a result of combining partitions

    typedef decltype(interior + border) InteriorBorder;
    typedef decltype(interior + border + overlap) InteriorBorderOverlap;
    typedef decltype(interior + border + overlap + front) InteriorBorderOverlapFront;
    typedef decltype(interior + border + overlap + front + ghost) All;

    namespace {

      // again, place the global objects in an anonymous namespace

      const InteriorBorder interiorBorder = {};
      const InteriorBorderOverlap interiorBorderOverlap = {};
      const InteriorBorderOverlapFront interiorBorderOverlapFront = {};
      const All all = {};

    }

#endif // DOXYGEN

  } // namespace Partitions

  /**
   * \} group Parallel Grid Partitions
   */

} // namespace Dune

#endif // DUNE_GRID_COMMON_PARTITIONSET_HH
