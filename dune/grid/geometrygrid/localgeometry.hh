// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_LOCALGEOMETRY_HH
#define DUNE_GEOGRID_LOCALGEOMETRY_HH

#include <dune/grid/genericgeometry/geometry.hh>

#include <dune/grid/geometrygrid/cornerstorage.hh>
#include <dune/grid/geometrygrid/numbering.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class Grid, int codim >
    struct LocalGeometryProvider;

  }



  // External Forward Declarations
  // -----------------------------

  template< class HostGrid, class CoordFunction, class Numbering, class Alloc >
  class GeometryGrid;



  namespace GenericGeometry
  {

    // LocalGeometryTraits for GeometryGrid
    // ------------------------------------

    template< class HostGrid, class CoordFunction, class Numbering, class Alloc >
    struct LocalGeometryTraits< GeometryGrid< HostGrid, CoordFunction, Numbering, Alloc > >
    {
      typedef GeometryGrid< HostGrid, CoordFunction, Numbering, Alloc > Grid;

      typedef DuneCoordTraits< typename HostGrid::ctype > CoordTraits;

      static const int dimGrid = HostGrid::dimension;
      static const int dimWorld = dimGrid;

      static const bool hybrid = true;
      // static const GeometryType::BasicType dunetype = GeometryType::simplex;

      static const GeometryType::BasicType linetype = GeometryType::simplex;

      template< class Topology >
      struct Mapping
      {
        typedef GeoGrid::CornerStorage< Topology, const Grid > CornerStorage;
        typedef CornerMapping< CoordTraits, Topology, dimWorld, CornerStorage > type;
      };

      struct Caching
      {
        static const EvaluationType evaluateJacobianTransposed = ComputeOnDemand;
        static const EvaluationType evaluateJacobianInverseTransposed = ComputeOnDemand;
        static const EvaluationType evaluateIntegrationElement = ComputeOnDemand;
        static const EvaluationType evaluateNormal = ComputeOnDemand;
      };

      typedef Alloc Allocator;
    };

  }



  namespace GeoGrid
  {

    template< class HostGrid, class Numbering, int codim >
    struct LocalCoordVector
    {
      typedef typename HostGrid::ctype ctype;

      static const int dimension = HostGrid::dimension;

      typedef FieldVector< ctype, dimension > Coordinate;

      typedef typename Numbering::template EntityNumbering< codim > EntityNumbering;
      typedef typename HostGrid::template Codim< codim >::LocalGeometry HostLocalGeometry;

      LocalCoordVector ( const HostLocalGeometry &hostLocalGeo, const EntityNumbering &numbering )
        : hostLocalGeo_( hostLocalGeo ),
          numbering_( numbering )
      {}

      template< unsigned int numCorners >
      void calculate ( Coordinate (&corners)[ numCorners ] ) const
      {
        const int subcodim = dimension - codim;
        assert( numCorners == hostLocalGeo_.numCorners() );
        for( unsigned int i = 0; i < numCorners; ++i )
        {
          const unsigned int j = numbering_.template map< Numbering::Forward >( subcodim, i );
          corners[ i ] = hostLocalGeo_.corner( j );
        }
      }

    private:
      const HostLocalGeometry &hostLocalGeo_;
      Numbering numbering_;
    };



    // LocalGeometryProvider
    // ---------------------

    template< class HostGrid, class CoordFunction, class Numbering, class Alloc, int codim >
    struct LocalGeometryProvider< GeometryGrid< HostGrid, CoordFunction, Numbering, Alloc >, codim >
    {
      typedef GeometryGrid< HostGrid, CoordFunction, Numbering, Alloc > Grid;

      static const int dimension = HostGrid::dimension;

      typedef typename HostGrid::template Codim< codim >::LocalGeometry HostLocalGeometry;
      typedef Dune::Geometry< dimension-codim, dimension, const Grid, GenericGeometry::LocalGeometry >
      LocalGeometry;

      typedef typename Numbering::template EntityNumbering< codim > EntityNumbering;

      typedef Alloc Allocator;

    private:
      typedef GeoGrid::LocalCoordVector< HostGrid, Numbering, codim > LocalCoordVector;
      typedef GenericGeometry::LocalGeometry< dimension-codim, dimension, const Grid > LocalGeometryImpl;

    public:
      explicit LocalGeometryProvider ( const Allocator &allocator = Allocator() )
        : allocator_( allocator ),
          geo_( LocalGeometryImpl() )
      {}

      const LocalGeometry &
      operator() ( const HostLocalGeometry &hostLocalGeo, const EntityNumbering &numbering ) const
      {
        LocalGeometryImpl &geo = Grid::getRealImplementation( geo_ );
        if( !geo )
        {
          const unsigned int topologyId = GenericGeometry::topologyId( hostLocalGeo.type() );
          LocalCoordVector coords( hostLocalGeo, numbering );
          geo = LocalGeometryImpl( topologyId, coords, allocator_ );
        }
        return geo_;
      }

    private:
      Allocator allocator_;
      mutable LocalGeometry geo_;
    };



    // LocalGeometryProvider for IdenticalNumbering
    // --------------------------------------------

    template< class HostGrid, class CoordFunction, class Alloc, int codim >
    struct LocalGeometryProvider< GeometryGrid< HostGrid, CoordFunction, IdenticalNumbering, Alloc >, codim >
    {
      typedef GeometryGrid< HostGrid, CoordFunction, IdenticalNumbering, Alloc > Grid;

      static const int dimension = HostGrid::dimension;

      typedef typename HostGrid::template Codim< codim >::LocalGeometry HostLocalGeometry;
      typedef HostLocalGeometry LocalGeometry;

      typedef typename IdenticalNumbering::template EntityNumbering< codim > EntityNumbering;

      typedef Alloc Allocator;

      explicit LocalGeometryProvider ( const Allocator &allocator = Allocator() )
      {}

      const LocalGeometry &
      operator() ( const HostLocalGeometry &hostLocalGeo, const EntityNumbering &numbering ) const
      {
        return hostLocalGeo;
      }
    };

  }

}

#endif // #ifndef DUNE_GEOGRID_LOCALGEOMETRY_HH
