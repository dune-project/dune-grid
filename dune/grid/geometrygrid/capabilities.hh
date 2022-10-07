// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_CAPABILITIES_HH
#define DUNE_GEOGRID_CAPABILITIES_HH

#include <cassert>
#include <type_traits>
#include <utility>

#include <dune/common/hybridutilities.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/geometrygrid/declaration.hh>

namespace Dune
{

  // Capabilities
  // ------------

  namespace Capabilities
  {

    // Capabilities from dune-grid
    // ---------------------------

    template< class HostGrid, class CoordFunction, class Allocator >
    struct hasSingleGeometryType< GeometryGrid< HostGrid, CoordFunction, Allocator > >
    {
      static const bool v = hasSingleGeometryType< HostGrid > :: v;
      static const unsigned int topologyId = hasSingleGeometryType< HostGrid > :: topologyId;
    };


    template< class HostGrid, class CoordFunction, class Allocator, int codim >
    struct hasEntity< GeometryGrid< HostGrid, CoordFunction, Allocator >, codim >
    {
      static const bool v = true;
    };


    template< class HostGrid, class CoordFunction, class Allocator, int codim >
    struct hasEntityIterator< GeometryGrid< HostGrid, CoordFunction, Allocator >, codim >
    {
      static const bool v = hasEntityIterator<HostGrid, codim>::v;
    };


    template< class HostGrid, class CoordFunction, class Allocator, int codim >
    struct canCommunicate< GeometryGrid< HostGrid, CoordFunction, Allocator >, codim >
    {
      static const bool v = canCommunicate< HostGrid, codim >::v && hasEntity< HostGrid, codim >::v;
    };


    template< class HostGrid, class CoordFunction, class Allocator >
    struct hasBackupRestoreFacilities< GeometryGrid< HostGrid, CoordFunction, Allocator > >
    {
      static const bool v = hasBackupRestoreFacilities< HostGrid >::v && std::is_default_constructible< CoordFunction >::value;
    };

    template< class HostGrid, class CoordFunction, class Allocator >
    struct isLevelwiseConforming< GeometryGrid< HostGrid, CoordFunction, Allocator > >
    {
      static const bool v = isLevelwiseConforming< HostGrid >::v;
    };

    template< class HostGrid, class CoordFunction, class Allocator >
    struct isLeafwiseConforming< GeometryGrid< HostGrid, CoordFunction, Allocator > >
    {
      static const bool v = isLeafwiseConforming< HostGrid >::v;
    };

    template< class HostGrid, class CoordFunction, class Allocator >
    struct threadSafe< GeometryGrid< HostGrid, CoordFunction, Allocator > >
    {
      static const bool v = false;
    };

    template< class HostGrid, class CoordFunction, class Allocator >
    struct viewThreadSafe< GeometryGrid< HostGrid, CoordFunction, Allocator > >
    {
      static const bool v = false;
    };




    // hasHostEntity
    // -------------

    template< class Grid, int codim >
    struct hasHostEntity;

    template< class Grid, int codim >
    struct hasHostEntity< const Grid, codim >
    {
      static const bool v = hasHostEntity< Grid, codim >::v;
    };

    template< class HostGrid, class CoordFunction, class Allocator, int codim >
    struct hasHostEntity< GeometryGrid< HostGrid, CoordFunction, Allocator >, codim >
    {
      static const bool v = hasEntity< HostGrid, codim >::v;
    };



    // CodimCache
    // ----------

    template< class Grid >
    class CodimCache
    {
      static const int dimension = Grid::dimension;

      bool hasHostEntity_[ Grid::dimension + 1 ];

      CodimCache ()
      {
        Hybrid::forEach( std::make_index_sequence< dimension+1 >{},
          [ & ]( auto i ){ hasHostEntity_[ i ] = Capabilities::hasHostEntity< Grid, i >::v; } );
      }

      static CodimCache &instance ()
      {
        static CodimCache singleton;
        return singleton;
      }

    public:
      static bool hasHostEntity ( int codim )
      {
        assert( (codim >= 0) && (codim <= dimension) );
        return instance().hasHostEntity_[ codim ];
      }
    };

  } // namespace Capabilities

} // namespace Dune

#endif // #ifndef DUNE_GEOGRID_CAPABILITIES_HH
