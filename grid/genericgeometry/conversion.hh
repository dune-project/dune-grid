// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_CONVERSION_HH
#define DUNE_GENERICGEOMETRY_CONVERSION_HH

#include <dune/common/static_assert.hh>
#include <dune/common/geometrytype.hh>

#include <dune/grid/genericgeometry/topologytypes.hh>

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

  }

}

#endif
