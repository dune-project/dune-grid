// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_CAPABILITIES_HH
#define DUNE_GEOGRID_CAPABILITIES_HH

#include <cassert>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/genericgeometry/misc.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class HostGrid, class CoordFunction >
  class GeometryGrid;



  // Capabilities
  // ------------

  /** \brief capability structures
   *
   *  In DUNE, capabilities are flags indicating support for some feature.
   *  They are set by specializing a structure for a certain implementation.
   *
   *  An example from dune-grid is:
   *  \code
   *  template< class Grid >
   *  struct isParallel
   *  {
   *    static bool v = false;
   *  };
   *  \endcode
   *  Here, the capability isParallel is defined and defaults to false.
   */
  namespace Capabilities
  {

    // Capabilities from dune-grid
    // ---------------------------

    template< class HostGrid, class CoordFunction, int codim >
    struct hasEntity< GeometryGrid< HostGrid, CoordFunction >, codim >
    {
      static const bool v = true;
      //static const bool v = ((codim == 0) || (codim == HostGrid :: dimension));
    };


    template< class HostGrid, class CoordFunction >
    struct isParallel< GeometryGrid< HostGrid, CoordFunction > >
    {
      static const bool v = isParallel< HostGrid > :: v;
    };


    template< class HostGrid, class CoordFunction >
    struct hasHangingNodes< GeometryGrid< HostGrid, CoordFunction > >
    {
      static const bool v = hasHangingNodes< HostGrid > :: v;
    };

    template< class HostGrid, class CoordFunction >
    struct hasBackupRestoreFacilities< GeometryGrid< HostGrid, CoordFunction > >
    {
      static const bool v = false;
    };

    template< class HostGrid, class CoordFunction >
    struct isLevelwiseConforming< GeometryGrid< HostGrid, CoordFunction > >
    {
      static const bool v = isLevelwiseConforming< HostGrid > :: v;
    };

    template< class HostGrid, class CoordFunction >
    struct isLeafwiseConforming< GeometryGrid< HostGrid, CoordFunction > >
    {
      static const bool v = isLeafwiseConforming< HostGrid > :: v;
    };

    template< class HostGrid, class CoordFunction >
    struct IsUnstructured< GeometryGrid< HostGrid, CoordFunction > >
    {
      static const bool v = true;
    };

    template< class HostGrid, class CoordFunction >
    struct threadSafe< GeometryGrid< HostGrid, CoordFunction > >
    {
      static const bool v = false;
    };

    template< class HostGrid, class CoordFunction >
    struct viewThreadSafe< GeometryGrid< HostGrid, CoordFunction > >
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
      static const bool v = hasHostEntity< Grid, codim > :: v;
    };

    template< class HostGrid, class CoordFunction, int codim >
    struct hasHostEntity< GeometryGrid< HostGrid, CoordFunction >, codim >
    {
      static const bool v = hasEntity< HostGrid, codim > :: v;
    };



    // CodimCache
    // ----------

    template< class Grid >
    class CodimCache
    {
      static const int dimension = Grid :: dimension;

      template< int codim >
      struct BuildCache;

      bool hasHostEntity_[ Grid :: dimension + 1 ];

      CodimCache ()
      {
        GenericGeometry :: ForLoop< BuildCache, 0, dimension >
        :: apply( hasHostEntity_ );
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
    struct CodimCache< Grid > :: BuildCache
    {
      static void apply ( bool (&hasHostEntity)[ dimension + 1 ] )
      {
        hasHostEntity[ codim ] = Capabilities :: hasHostEntity< Grid, codim > :: v;
      }
    };

  }

}

#endif
