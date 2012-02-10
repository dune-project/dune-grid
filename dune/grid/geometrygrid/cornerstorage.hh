// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_CORNERSTORAGE_HH
#define DUNE_GEOGRID_CORNERSTORAGE_HH

#include <dune/geometry/genericgeometry/geometry.hh>

#include <dune/grid/geometrygrid/hostcorners.hh>
#include <dune/grid/geometrygrid/coordfunction.hh>

namespace Dune
{

  namespace GeoGrid
  {

    // CoordFunctionCaller
    // -------------------

    template< class HostEntity, class CoordFunctionInterface >
    class CoordFunctionCaller;

    template< class HostEntity, class ct, unsigned int dimD, unsigned int dimR, class Impl >
    class CoordFunctionCaller< HostEntity, AnalyticalCoordFunctionInterface< ct, dimD, dimR, Impl > >
    {
      typedef AnalyticalCoordFunctionInterface< ct, dimD, dimR, Impl > CoordFunctionInterface;
      typedef CoordFunctionCaller< HostEntity, CoordFunctionInterface > This;

      static const int codimension = HostEntity::codimension;

    public:
      typedef typename CoordFunctionInterface::RangeVector RangeVector;

      CoordFunctionCaller ( const HostEntity &hostEntity,
                            const CoordFunctionInterface &coordFunction )
        : hostCorners_( hostEntity ),
          coordFunction_( coordFunction )
      {}

      void evaluate ( unsigned int i, RangeVector &y ) const
      {
        coordFunction_.evaluate( hostCorners_.corner( i ), y );
      }

      GeometryType type () const
      {
        return hostCorners_.type();
      }

      unsigned int numCorners () const
      {
        return hostCorners_.numCorners();
      }

    private:
      const HostCorners< HostEntity > hostCorners_;
      const CoordFunctionInterface &coordFunction_;
    };

    template< class HostEntity, class ct, unsigned int dimR, class Impl >
    class CoordFunctionCaller< HostEntity, DiscreteCoordFunctionInterface< ct, dimR, Impl > >
    {
      typedef DiscreteCoordFunctionInterface< ct, dimR, Impl > CoordFunctionInterface;
      typedef CoordFunctionCaller< HostEntity, CoordFunctionInterface > This;

      typedef typename CoordFunctionInterface::RangeVector RangeVector;

    public:
      CoordFunctionCaller ( const HostEntity &hostEntity,
                            const CoordFunctionInterface &coordFunction )
        : hostEntity_( hostEntity ),
          coordFunction_( coordFunction )
      {}

      void evaluate ( unsigned int i, RangeVector &y ) const
      {
        coordFunction_.evaluate( hostEntity_, i, y );
      }

      GeometryType type () const
      {
        return hostEntity_.type();
      }

      unsigned int numCorners () const
      {
        return hostEntity_.geometry().corners();
      }

    private:
      const HostEntity &hostEntity_;
      const CoordFunctionInterface &coordFunction_;
    };



    // CoordVector
    // -----------

    template< int mydim, class Grid, bool fake >
    class CoordVector;


    template< int mydim, class Grid >
    class CoordVector< mydim, Grid, false >
    {
      typedef typename remove_const< Grid >::type::Traits Traits;

      typedef typename Traits::ctype ctype;

      static const int dimension = Traits::dimension;
      static const int mydimension = mydim;
      static const int codimension = dimension - mydimension;
      static const int dimensionworld = Traits::dimensionworld;

      typedef FieldVector< ctype, dimensionworld > Coordinate;

      typedef typename Traits::HostGrid HostGrid;
      typedef typename Traits::CoordFunction CoordFunction;

      typedef typename HostGrid::template Codim< codimension >::Entity HostEntity;

      typedef GeoGrid :: CoordFunctionCaller< HostEntity, typename CoordFunction::Interface >
      CoordFunctionCaller;

    public:
      CoordVector ( const HostEntity &hostEntity,
                    const CoordFunction &coordFunction )
        : coordFunctionCaller_( hostEntity, coordFunction )
      {}

      template< unsigned int numCorners >
      void calculate ( Coordinate (&corners)[ numCorners ] ) const
      {
        assert( numCorners == coordFunctionCaller_.numCorners() );
        for( unsigned int i = 0; i < numCorners; ++i )
          coordFunctionCaller_.evaluate( i, corners[ i ] );
      }

    private:
      const CoordFunctionCaller coordFunctionCaller_;
    };


    template< int mydim, class Grid >
    class CoordVector< mydim, Grid, true >
    {
      typedef typename remove_const< Grid > :: type :: Traits Traits;

      typedef typename Traits::ctype ctype;

      static const int dimension = Traits::dimension;
      static const int mydimension = mydim;
      static const int codimension = dimension - mydimension;
      static const int dimensionworld = Traits::dimensionworld;

      typedef FieldVector< ctype, dimensionworld > Coordinate;

      typedef typename Traits::HostGrid HostGrid;
      typedef typename Traits::CoordFunction CoordFunction;

      typedef typename HostGrid::template Codim< 0 >::Entity HostElement;

      typedef GeoGrid::CoordFunctionCaller< HostElement, typename CoordFunction::Interface >
      CoordFunctionCaller;

    public:
      CoordVector ( const HostElement &hostElement,
                    const unsigned int subEntity,
                    const CoordFunction &coordFunction )
        : coordFunctionCaller_( hostElement, coordFunction ),
          subEntity_( subEntity )
      {}

      template< unsigned int numCorners >
      void calculate ( Coordinate (&corners)[ numCorners ] ) const
      {
        const GeometryType type = coordFunctionCaller_.type();
        const GenericReferenceElement< ctype, dimension > &refElement
          = GenericReferenceElements< ctype, dimension >::general( type );
        assert( numCorners == refElement.size( subEntity_, codimension, dimension ) );

        for( unsigned int i = 0; i < numCorners; ++i )
        {
          const unsigned int j = refElement.subEntity( subEntity_, codimension, i, dimension );
          coordFunctionCaller_.evaluate( j, corners[ i ] );
        }
      }

    private:
      const CoordFunctionCaller coordFunctionCaller_;
      const unsigned int subEntity_;
    };



    // IntersectionCoordVector
    // -----------------------

    template< class Grid >
    class IntersectionCoordVector
    {
      typedef typename remove_const< Grid >::type::Traits Traits;

      typedef typename Traits::ctype ctype;

      static const int dimension = Traits::dimension;
      static const int codimension = 1;
      static const int mydimension = dimension-codimension;
      static const int dimensionworld = Traits::dimensionworld;

      typedef FieldVector< ctype, dimensionworld > Coordinate;

      typedef typename Traits::HostGrid HostGrid;

      typedef typename Traits::template Codim< 0 >::GeometryImpl ElementGeometryImpl;
      typedef typename Traits::template Codim< codimension >::LocalGeometry HostLocalGeometry;

    public:
      IntersectionCoordVector ( const ElementGeometryImpl &elementGeometry,
                                const HostLocalGeometry &hostLocalGeometry )
        : elementGeometry_( elementGeometry ),
          hostLocalGeometry_( hostLocalGeometry )
      {}

      template< unsigned int numCorners >
      void calculate ( Coordinate (&corners)[ numCorners ] ) const
      {
        assert( numCorners == hostLocalGeometry_.corners() );
        for( unsigned int i = 0; i < numCorners; ++i )
          corners[ i ] = elementGeometry_.global( hostLocalGeometry_.corner( i ) );
      }

    private:
      const ElementGeometryImpl &elementGeometry_;
      HostLocalGeometry hostLocalGeometry_;
    };




    // CornerStorage
    // -------------

    template< class Topology, class Grid >
    class CornerStorage
    {
      typedef typename remove_const< Grid >::type::Traits Traits;

      typedef typename Traits::ctype ctype;

      static const int dimension = Traits::dimension;
      static const int mydimension = Topology::dimension;
      static const int codimension = dimension - mydimension;
      static const int dimensionworld = Traits::dimensionworld;

      typedef FieldVector< ctype, dimensionworld > Coordinate;

    public:
      static const unsigned int size = Topology::numCorners;

      template< class SubTopology >
      struct SubStorage
      {
        typedef CornerStorage< SubTopology, Grid > type;
      };

      template< bool fake >
      explicit
      CornerStorage ( const CoordVector< mydimension, Grid, fake > &coords )
      {
        coords.calculate( coords_ );
      }

      explicit CornerStorage ( const IntersectionCoordVector< Grid > &coords )
      {
        coords.calculate( coords_ );
      }

      template< class Mapping, unsigned int codim >
      explicit
      CornerStorage ( const GenericGeometry::SubMappingCoords< Mapping, codim > &coords )
      {
        for( unsigned int i = 0; i < size; ++i )
          coords_[ i ] = coords[ i ];
      }

      const Coordinate &operator[] ( unsigned int i ) const
      {
        return coords_[ i ];
      }

    private:
      Coordinate coords_[ size ];
    };

  }

}

#endif // #ifndef DUNE_GEOGRID_CORNERSTORAGE_HH
