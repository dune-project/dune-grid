// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_GRIDVIEWPARTITIONING_INTERFACE_HH
#define DUNE_GRID_UTILITY_GRIDVIEWPARTITIONING_INTERFACE_HH

#include <cstddef>

#include <dune/common/documentation.hh>

namespace Dune {

  template<class HostGridView>
  class GridViewPartitioningInterface
  {
  public:
    typedef ImplementationDefined GridView;
    std::size_t partitions() const;
    GridView gridview(std::size_t partition) const;
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_GRIDVIEWPARTITIONING_INTERFACE_HH
