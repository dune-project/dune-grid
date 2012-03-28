// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_PERSISTENTCONTAINER_HH
#define DUNE_GEOGRID_PERSISTENTCONTAINER_HH

#include <dune/grid/geometrygrid/declaration.hh>
#include <dune/grid/utility/persistentcontainerwrapper.hh>

namespace Dune
{

  // PersistentContainer for GeometryGrid
  // ------------------------------------

  template< class HostGrid, class CoordFunction, class Data, class Allocator >
  class PersistentContainer< GeometryGrid< HostGrid, CoordFunction, Allocator >, Data, Allocator >
    : public PersistentContainerWrapper< GeometryGrid< HostGrid, CoordFunction, Allocator >, Data, Allocator >
  {
    typedef PersistentContainerWrapper< GeometryGrid< HostGrid, CoordFunction, Allocator >, Data, Allocator > Base;

  public:
    typedef typename Base::Grid Grid;

    PersistentContainer ( const Grid &grid, const int codim, const Allocator &allocator = Allocator() )
      : Base( grid, codim, allocator )
    {}
  };

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_PERSISTENTCONTAINER_HH
