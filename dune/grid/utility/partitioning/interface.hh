// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_PARTITIONING_INTERFACE_HH
#define DUNE_GRID_UTILITY_PARTITIONING_INTERFACE_HH

#include <dune/common/documentation.hh>

namespace Dune {

  //! Basic interface for partitionings
  class PartitioningInterface
  {
  public:
    //! type of partitions
    typedef ImplementationDefined Partition;
    //! type used to count partitions
    typedef ImplementationDefined Size;

    //! return number of partitions
    Size partitions() const;

    //! return a particular partition
    Partition partition(Size pId) const;
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_PARTITIONING_INTERFACE_HH
