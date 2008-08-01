// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_CONVERSION_HH
#define DUNE_GENERICGEOMETRY_CONVERSION_HH

#include <limits>

#include <dune/common/static_assert.hh>
#include <dune/common/geometrytype.hh>

#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/genericgeometry/subtopologies.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // DuneGeometryType
    // ----------------

    template< class Topology, GeometryType :: BasicType defaultType >
    class DuneGeometryType;

    template< GeometryType :: BasicType defaultType >
    class DuneGeometryType< Point, defaultType >
    {
      dune_static_assert( (defaultType == GeometryType :: simplex)
                          || (defaultType == GeometryType :: cube),
                          "defaultType may only be a simplex or a cube." );

    public:
      enum { dimension = 0 };
      enum { basicType = defaultType };

      static GeometryType type ()
      {
        return GeometryType( (GeometryType :: BasicType)basicType, dimension );
      }
    };

    template< class BaseTopology, GeometryType :: BasicType defaultType >
    class DuneGeometryType< Prism< BaseTopology >, defaultType >
    {
      typedef DuneGeometryType< BaseTopology, defaultType > DuneBaseGeometryType;

      dune_static_assert( ((int)defaultType == (int)GeometryType :: simplex)
                          || ((int)defaultType == (int)GeometryType :: cube),
                          "defaultType may only be a simplex or a cube." );

      dune_static_assert( ((int)DuneBaseGeometryType :: basicType == (int)GeometryType :: simplex)
                          || ((int)DuneBaseGeometryType :: basicType == (int)GeometryType :: cube),
                          "Only prisms over simplices or cubes can be converted." );

    public:
      enum { dimension = DuneBaseGeometryType :: dimension + 1 };
      enum
      {
        basicType = (dimension == 1)
                    ? defaultType
                    : ((dimension == 2)
                       || ((int)DuneBaseGeometryType :: basicType == (int)GeometryType :: cube))
                    ? GeometryType :: cube
                    : GeometryType :: prism
      };

      static GeometryType type ()
      {
        return GeometryType( (GeometryType :: BasicType)basicType, dimension );
      }
    };

    template< class BaseTopology, GeometryType :: BasicType defaultType >
    class DuneGeometryType< Pyramid< BaseTopology >, defaultType >
    {
      typedef DuneGeometryType< BaseTopology, defaultType > DuneBaseGeometryType;

      dune_static_assert( ((int)defaultType == (int)GeometryType :: simplex)
                          || ((int)defaultType == (int)GeometryType :: cube),
                          "defaultType may only be a simplex or a cube." );

      dune_static_assert( ((int)DuneBaseGeometryType :: basicType == (int)GeometryType :: simplex)
                          || ((int)DuneBaseGeometryType :: basicType == (int)GeometryType :: cube),
                          "Only pyramids over simplices or cubes can be converted." );

    public:
      enum { dimension = DuneBaseGeometryType :: dimension + 1 };
      enum
      {
        basicType = (dimension == 1)
                    ? defaultType
                    : ((dimension == 2)
                       || ((int)DuneBaseGeometryType :: basicType == (int)GeometryType :: simplex))
                    ? GeometryType :: simplex
                    : GeometryType :: pyramid
      };

      static GeometryType type ()
      {
        return GeometryType( (GeometryType :: BasicType)basicType, dimension );
      }
    };



    // DuneGeometryTypeProvider
    // ------------------------

    template< unsigned int dim, GeometryType :: BasicType defaultType >
    struct DuneGeometryTypeProvider
    {
      enum { dimension = dim };
      enum { numTopologies = (1 << dimension) };

    private:
      template< int i >
      struct Builder;

      GeometryType types_[ numTopologies ];

      DuneGeometryTypeProvider ()
      {
        ForLoop< Builder, 0, (1 << dim)-1 > :: apply( types_ );
      }

      static const DuneGeometryTypeProvider &instance ()
      {
        static DuneGeometryTypeProvider inst;
        return inst;
      }

    public:
      static const GeometryType &type ( unsigned int topologyId )
      {
        assert( topologyId < numTopologies );
        return instance().types_[ topologyId ];
      }
    };


    template< unsigned int dim, GeometryType :: BasicType defaultType >
    template< int i >
    struct DuneGeometryTypeProvider< dim, defaultType > :: Builder
    {
      typedef typename GenericGeometry :: Topology< i, dimension > :: type Topology;
      typedef GenericGeometry :: DuneGeometryType< Topology, defaultType >
      DuneGeometryType;

      static void apply ( GeometryType (&types)[ numTopologies ] )
      {
        types[ i ] = DuneGeometryType :: type();
      }
    };



    // Convert
    // -------

    template< GeometryType :: BasicType type, unsigned int dim >
    struct Convert;

    template< unsigned int dim >
    struct Convert< GeometryType :: simplex, dim >
    {
      typedef Pyramid
      < typename Convert< GeometryType :: simplex, dim-1 > :: type >
      type;

      template< unsigned int codim >
      static unsigned int map ( unsigned int i )
      {
        static unsigned int tetra_edge[ 6 ] = { 0, 2, 1, 3, 4, 5 };

        if( (dim == 2) && (codim == 1) )
        {
          return 2 - i;
        }

        if( dim == 3 )
        {
          if( codim == 1 )
            return 3 - i;
          if( codim == 2 )
            return tetra_edge[ i ];
        }

        return i;
      }
    };

    template<>
    struct Convert< GeometryType :: simplex, 0 >
    {
      typedef Point type;

      template< unsigned int codim >
      static unsigned int map ( unsigned int i )
      {
        return i;
      }
    };

    template< unsigned int dim >
    struct Convert< GeometryType :: cube, dim >
    {
      typedef Prism< typename Convert< GeometryType :: cube, dim-1 > :: type >
      type;

      template< unsigned int codim >
      static unsigned int map ( unsigned int i )
      {
        static unsigned int edge_3d[ 12 ] = { 0, 1, 2, 3, 4, 5, 8, 9, 6, 7, 10, 11 };

        if( dim > 3 )
          DUNE_THROW( NotImplemented, "Cube renumbering is not implemented for dim > 3." );

        if( (dim == 3) && (codim == 2) )
          return edge_3d[ i ];

        return i;
      }
    };

    template<>
    struct Convert< GeometryType :: cube, 0 >
    {
      typedef Point type;

      template< unsigned int codim >
      static unsigned int map ( unsigned int i )
      {
        return i;
      }
    };

    template< unsigned int dim >
    struct Convert< GeometryType :: prism, dim >
    {
      typedef Prism
      < typename Convert< GeometryType :: simplex, dim-1 > :: type >
      type;

      template< unsigned int codim >
      static unsigned int map ( unsigned int i )
      {
        DUNE_THROW( NotImplemented, "Cube renumbering is not implemented for dim > 3." );
      }

    private:
      dune_static_assert( dim >= 3, "Dune prisms must be at least 3-dimensional." );
    };

    template< unsigned int dim >
    struct Convert< GeometryType :: pyramid, dim >
    {
      typedef Pyramid
      < typename Convert< GeometryType :: cube, dim-1 > :: type >
      type;

      // Note that we map dune numbering into the generic one
      // this is only important for pyramids
      template< unsigned int codim >
      static unsigned int map ( unsigned int i )
      {
        static unsigned int vertex_3d[ 5 ] = { 0, 1, 3, 2, 4 };
        static unsigned int edge_3d[ 8 ] = { 2, 1, 3, 0, 4, 5, 7, 6 };
        static unsigned int face_3d[ 5 ] = { 0, 3, 2, 4, 1 };

        if( dim > 3 )
          DUNE_THROW( NotImplemented, "Cube renumbering is not implemented for dim > 3." );

        if( codim == dim )
          return vertex_3d[ i ];
        if( codim == 2 )
          return edge_3d[ i ];
        if( codim == 1 )
          return face_3d[ i ];

        return i;
      }

    private:
      dune_static_assert( dim >= 3, "Dune pyramids must be at least 3-dimensional." );
    };



    // MapNumbering
    // ------------

    template< class Topology >
    struct MapNumbering;

    struct MapNumberingIdentical
    {
      template< unsigned int codim >
      static unsigned int dune2generic ( unsigned int i )
      {
        return i;
      }
    };

    // MapNumbering for Point
    template<>
    struct MapNumbering< Point >
      : public MapNumberingIdentical
    {};

    // MapNumbering for Line
    template<>
    struct MapNumbering< Prism< Point > >
      : public MapNumberingIdentical
    {};

    template<>
    struct MapNumbering< Pyramid< Point > >
      : public MapNumberingIdentical
    {};

    // MapNumbering for Triangle
    struct MapNumberingTriangle
    {
      template< unsigned int codim >
      static unsigned int dune2generic ( unsigned int i )
      {
        return (codim == 1 ? 2 - i : i);
      }
    };

    template<>
    struct MapNumbering< Pyramid< Pyramid< Point > > >
      : public MapNumberingTriangle
    {};

    template<>
    struct MapNumbering< Pyramid< Prism< Point > > >
      : public MapNumberingTriangle
    {};


    // MapNumbering for Quadrilateral
    template<>
    struct MapNumbering< Prism< Pyramid< Point > > >
      : public MapNumberingIdentical
    {};

    template<>
    struct MapNumbering< Prism< Prism< Point > > >
      : public MapNumberingIdentical
    {};

    // MapNumbering for Tetrahedron
    struct MapNumberingTetrahedron
    {
      template< unsigned int codim >
      static unsigned int dune2generic ( unsigned int i )
      {
        static unsigned int edge[ 6 ] = { 0, 2, 1, 3, 4, 5 };
        return (codim == 1 ? 3 - i : (codim == 2 ? edge[ i ] : i));
      }
    };

    template<>
    struct MapNumbering< Pyramid< Pyramid< Pyramid< Point > > > >
      : public MapNumberingTetrahedron
    {};

    template<>
    struct MapNumbering< Pyramid< Pyramid< Prism< Point > > > >
      : public MapNumberingTetrahedron
    {};

    // MapNumbering for Cube
    struct MapNumberingCube
    {
      template< unsigned int codim >
      static unsigned int dune2generic ( unsigned int i )
      {
        static unsigned int edge[ 12 ] = { 0, 1, 2, 3, 4, 5, 8, 9, 6, 7, 10, 11 };
        return (codim == 2 ? edge[ i ] : i);
      }
    };

    template<>
    struct MapNumbering< Prism< Prism< Pyramid< Point > > > >
      : public MapNumberingCube
    {};

    template<>
    struct MapNumbering< Prism< Prism< Prism< Point > > > >
      : public MapNumberingCube
    {};

    // MapNumbering for Pyramid
    struct MapNumberingPyramid
    {
      template< unsigned int codim >
      static unsigned int dune2generic ( unsigned int i )
      {
        static unsigned int vertex[ 5 ] = { 0, 1, 3, 2, 4 };
        static unsigned int edge[ 8 ] = { 2, 1, 3, 0, 4, 5, 7, 6 };
        static unsigned int face[ 5 ] = { 0, 3, 2, 4, 1 };

        if( codim == 3 )
          return vertex[ i ];
        else if( codim == 2 )
          return edge[ i ];
        else if( codim == 1 )
          return face[ i ];
        else
          return i;
      }
    };

    template<>
    struct MapNumbering< Pyramid< Prism< Pyramid< Point > > > >
      : public MapNumberingPyramid
    {};

    template<>
    struct MapNumbering< Pyramid< Prism< Prism< Point > > > >
      : public MapNumberingPyramid
    {};

    // MapNumbering for Prism
    struct MapNumberingPrism
    {
      template< unsigned int codim >
      static unsigned int dune2generic ( unsigned int i )
      {
        return (codim == 0 ? i : std :: numeric_limits< unsigned int > :: max());
      }
    };

    template<>
    struct MapNumbering< Prism< Pyramid< Pyramid< Point > > > >
      : public MapNumberingPrism
    {};

    template<>
    struct MapNumbering< Prism< Pyramid< Prism< Point > > > >
      : public MapNumberingPrism
    {};



    // MapNumberingProvider
    // --------------------

    template< unsigned int dim >
    class MapNumberingProvider
    {
      enum { dimension = dim };
      enum { numTopologies = (1 << dimension) };

    private:
      template< int i >
      struct Builder;

      typedef std :: vector< unsigned int > Map;

      Map dune2generic_[ numTopologies ][ dimension+1 ];

      MapNumberingProvider ()
      {
        ForLoop< Builder, 0, (1 << dim)-1 > :: apply( dune2generic_ );
      }

      static const MapNumberingProvider &instance ()
      {
        static MapNumberingProvider inst;
        return inst;
      }

    public:
      template< unsigned int codim >
      static unsigned int
      dune2generic ( unsigned int topologyId, unsigned int i )
      {
        assert( topologyId < numTopologies );
        return instance().dune2generic_[ topologyId ][ codim ][ i ];
      }
    };


    template< unsigned int dim >
    template< int i >
    struct MapNumberingProvider< dim > :: Builder
    {
      typedef typename GenericGeometry :: Topology< i, dimension > :: type Topology;
      typedef GenericGeometry :: MapNumbering< Topology > MapNumbering;

      template< int codim >
      struct Codim;

      static void apply ( Map (&dune2generic)[ numTopologies ][ dimension+1 ] )
      {
        ForLoop< Codim, 0, dimension > :: apply( dune2generic[ i ] );
      }
    };

    template< unsigned int dim >
    template< int i >
    template< int codim >
    struct MapNumberingProvider< dim > :: Builder< i > :: Codim
    {
      static void apply ( Map (&dune2generic)[ dimension+1 ] )
      {
        const unsigned int size = Size< Topology, codim > :: value;
        Map &d2g = dune2generic[ codim ];

        d2g.resize( size );
        for( unsigned int j = 0; j < size; ++j )
          d2g[ j ] = MapNumbering :: template dune2generic< codim >( j );
      }
    };

  }

}

#endif
