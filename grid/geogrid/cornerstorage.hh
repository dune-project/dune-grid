// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOGRID_CORNERSTORAGE_HH
#define DUNE_GEOGRID_CORNERSTORAGE_HH

#include <dune/grid/genericgeometry/geometry.hh>

namespace Dune
{

  // GeometryGridCoordVector
  // -----------------------

  template< int mydim, class Grid, bool fake >
  class GeometryGridCoordVector;


  template< int mydim, class Grid >
  class GeometryGridCoordVector< mydim, Grid, false >
  {
    typedef typename remove_const< Grid > :: type :: Traits Traits;

    typedef typename Traits :: ctype ctype;

    static const int dimension = Traits :: dimension;
    static const int mydimension = mydim;
    static const int codimension = dimension - mydimension;
    static const int dimensionworld = Traits :: dimensionworld;

    typedef FieldVector< ctype, dimensionworld > Coordinate;

    typedef typename Traits :: HostGrid HostGrid;
    typedef typename Traits :: CoordFunction CoordFunction;

    typedef typename HostGrid :: template Codim< codimension > :: Geometry HostGeometry;

  private:
    const HostGeometry &hostGeometry_;
    const CoordFunction &coordFunction_;

  public:
    GeometryGridCoordVector ( const HostGeometry &hostGeometry,
                              const CoordFunction &coordFunction )
      : hostGeometry_( hostGeometry ),
        coordFunction_( coordFunction )
    {}

    template< unsigned int numCorners >
    void calculate ( Coordinate (&corners)[ numCorners ] ) const
    {
      assert( numCorners == hostGeometry_.corners() );
      for( unsigned int i = 0; i < numCorners; ++i )
        coordFunction_.evaluate( hostGeometry_[ i ], corners[ i ] );
    }
  };


  template< int mydim, class Grid >
  class GeometryGridCoordVector< mydim, Grid, true >
  {
    typedef typename remove_const< Grid > :: type :: Traits Traits;

    typedef typename Traits :: ctype ctype;

    static const int dimension = Traits :: dimension;
    static const int mydimension = mydim;
    static const int codimension = dimension - mydimension;
    static const int dimensionworld = Traits :: dimensionworld;

    typedef FieldVector< ctype, dimensionworld > Coordinate;

    typedef typename Traits :: HostGrid HostGrid;
    typedef typename Traits :: CoordFunction CoordFunction;

    typedef typename HostGrid :: template Codim< 0 > :: Geometry HostGeometry;

  private:
    const HostGeometry &hostGeometry_;
    const unsigned int subEntity_;
    const CoordFunction &coordFunction_;

  public:
    GeometryGridCoordVector ( const HostGeometry &hostGeometry,
                              const unsigned int subEntity,
                              const CoordFunction &coordFunction )
      : hostGeometry_( hostGeometry ),
        subEntity_( subEntity ),
        coordFunction_( coordFunction )
    {}

    template< unsigned int numCorners >
    void calculate ( Coordinate (&corners)[ numCorners ] ) const
    {
      const ReferenceElement< ctype, dimension > &refElement
        = ReferenceElements< ctype, dimension > :: general( hostGeometry_.type() );
      assert( numCorners == refElement.size( subEntity_, codimension, dimension ) );

      for( unsigned int i = 0; i < numCorners; ++i )
      {
        const int j = refElement.subEntity( subEntity_, codimension, i, dimension );
        coordFunction_.evaluate( hostGeometry_[ j ], corners[ i ] );
      }
    }
  };



  // GeometryGridCornerStorage
  // -------------------------

  template< class Topology, class Grid >
  class GeometryGridCornerStorage
  {
    typedef typename remove_const< Grid > :: type :: Traits Traits;

    typedef typename Traits :: ctype ctype;

    static const int dimension = Traits :: dimension;
    static const int mydimension = Topology :: dimension;
    static const int codimension = dimension - mydimension;
    static const int dimensionworld = Traits :: dimensionworld;

    typedef FieldVector< ctype, dimensionworld > Coordinate;

  public:
    static const unsigned int size = Topology :: numCorners;

  private:
    Coordinate coords_[ size ];

  public:
    template< bool fake >
    explicit
    GeometryGridCornerStorage ( const GeometryGridCoordVector< mydimension, Grid, fake > &coords )
    {
      coords.calculate( coords_ );
    }

    template< class Mapping, unsigned int codim >
    explicit
    GeometryGridCornerStorage ( const GenericGeometry :: SubMappingCoords< Mapping, codim > &coords )
    {
      for( unsigned int i = 0; i < size; ++i )
        coords_[ i ] = coords[ i ];
    }

    const Coordinate &operator[] ( unsigned int i ) const
    {
      return coords_[ i ];
    }
  };

}

#endif
