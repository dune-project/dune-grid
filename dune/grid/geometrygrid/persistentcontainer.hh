// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_PERSISTENTCONTAINER_HH
#define DUNE_GEOGRID_PERSISTENTCONTAINER_HH

#include <dune/grid/geometrygrid/declaration.hh>
#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/grid/utility/persistentcontainerwrapper.hh>

namespace Dune
{

  // PersistentContainer for GeometryGrid
  // ------------------------------------

  template< class HostGrid, class CoordFunction, class Allocator, class T >
  class PersistentContainer< GeometryGrid< HostGrid, CoordFunction, Allocator >, T >
    : public PersistentContainerWrapper< GeometryGrid< HostGrid, CoordFunction, Allocator >, T >
  {
    typedef PersistentContainerWrapper< GeometryGrid< HostGrid, CoordFunction, Allocator >, T > Base;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    PersistentContainer ( const Grid &grid, int codim, const Value &value = Value() )
      : Base( grid, codim, value )
    {}
  };

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_PERSISTENTCONTAINER_HH
