// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_CAPABILITIES_HH
#define DUNE_GEOGRID_CAPABILITIES_HH

#include <cassert>

#include <dune/common/forloop.hh>

#include <dune/grid/common/capabilities.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class HostGrid, class CoordFunction, class Allocator >
  class GeometryGrid;



  // Capabilities
  // ------------

  namespace Capabilities
  {

    // Capabilities from dune-grid
    // ---------------------------

    template< class HostGrid, class CoordFunction, class Allocator, int codim >
    struct hasEntity< GeometryGrid< HostGrid, CoordFunction, Allocator >, codim >
    {
      static const bool v = true;
    };


    template< class HostGrid, class CoordFunction, class Allocator >
    struct isParallel< GeometryGrid< HostGrid, CoordFunction, Allocator > >
    {
      static const bool v = isParallel< HostGrid >::v;
    };


    template< class HostGrid, class CoordFunction, class Allocator, int codim >
    struct canCommunicate< GeometryGrid< HostGrid, CoordFunction, Allocator >, codim >
    {
      static const bool v = canCommunicate< HostGrid, codim >::v;
    };


    template< class HostGrid, class CoordFunction, class Allocator >
    struct hasBackupRestoreFacilities< GeometryGrid< HostGrid, CoordFunction, Allocator > >
    {
      static const bool v = hasBackupRestoreFacilities< HostGrid >::v;
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



    // non-standard capabilities
    // -------------------------

    template< class Grid >
    struct hasHierarchicIndexSet;

    template< class HostGrid, class CoordFunction, class Allocator >
    struct hasHierarchicIndexSet< GeometryGrid< HostGrid, CoordFunction, Allocator > >
    {
      static const bool v = hasHierarchicIndexSet< HostGrid >::v;
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

      template< int codim >
      struct BuildCache;

      bool hasHostEntity_[ Grid::dimension + 1 ];

      CodimCache ()
      {
        Dune::ForLoop< BuildCache, 0, dimension >::apply( hasHostEntity_ );
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

    template< class Grid >
    template< int codim >
    struct CodimCache< Grid >::BuildCache
    {
      static void apply ( bool (&hasHostEntity)[ dimension + 1 ] )
      {
        hasHostEntity[ codim ] = Capabilities::hasHostEntity< Grid, codim >::v;
      }
    };

  }

}

#endif // #ifndef DUNE_GEOGRID_CAPABILITIES_HH
