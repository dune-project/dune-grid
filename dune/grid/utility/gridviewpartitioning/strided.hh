// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_GRIDVIEWPARTITIONING_STRIDED_HH
#define DUNE_GRID_UTILITY_GRIDVIEWPARTITIONING_STRIDED_HH

#include <cstddef>

#include <dune/grid/common/gridview.hh>

#include <dune/grid/utility/entityfilter.hh>
#include <dune/grid/utility/filteringgridview.hh>

namespace Dune {

  template<class HostView>
  class StridedGridViewPartitioning
  {
    typedef StridedEntityFilter<
      typename HostView::IndexSet,
      typename HostView::template Codim<0>::Entity> Filter;

    HostView host_;
    std::size_t partitions_;

  public:
    StridedGridViewPartitioning(const HostView &host,
                                std::size_t partitions) :
      host_(host), partitions_(partitions)
    { }

    typedef FilteringGridView<Filter, HostView> GridView;
    std::size_t partitions() const
    {
      return partitions_;
    }
    GridView gridView(std::size_t partition) const
    {
      return GridView(Filter(host_.indexSet(), partitions_, partition),
                      host_);
    }
  };

} // namespace Dune

#endif // DUNE_GRID_UTILITY_GRIDVIEWPARTITIONING_STRIDED_HH
