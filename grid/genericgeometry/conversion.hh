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
    };

    template< class BaseTopology, GeometryType :: BasicType defaultType >
    class DuneGeometryType< Prism< BaseTopology >, defaultType >
    {
      typedef DuneGeometryType< BaseTopology, defaultType > DuneBaseGeometryType;

      dune_static_assert( (defaultType == GeometryType :: simplex)
                          || (defaultType == GeometryType :: cube),
                          "defaultType may only be a simplex or a cube." );

      dune_static_assert( (DuneBaseGeometryType :: basicType == GeometryType :: simplex)
                          || (DuneBaseGeometryType :: basicType == GeometryType :: cube),
                          "Only prisms over simplices or cubes can be converted." );

    public:
      enum { dimension = DuneBaseGeometryType :: dimension + 1 };
      enum
      {
        basicType = (dimension == 1)
                    ? defaultType
                    : ((dimension == 2)
                       || (DuneBaseGeometryType :: basicType == GeometryType :: cube))
                    ? GeometryType :: cube
                    : GeometryType :: prism
      };
    };

    template< class BaseTopology, GeometryType :: BasicType defaultType >
    class DuneGeometryType< Pyramid< BaseTopology >, defaultType >
    {
      typedef DuneGeometryType< BaseTopology, defaultType > DuneBaseGeometryType;

      dune_static_assert( (defaultType == GeometryType :: simplex)
                          || (defaultType == GeometryType :: cube),
                          "defaultType may only be a simplex or a cube." );

      dune_static_assert( (DuneBaseGeometryType :: basicType == GeometryType :: simplex)
                          || (DuneBaseGeometryType :: basicType == GeometryType :: cube),
                          "Only pyramids over simplices or cubes can be converted." );

    public:
      enum { dimension = DuneBaseGeometryType :: dimension + 1 };
      enum
      {
        basicType = (dimension == 1)
                    ? defaultType
                    : ((dimension == 2)
                       || (DuneBaseGeometryType :: basicType == GeometryType :: simplex))
                    ? GeometryType :: simplex
                    : GeometryType :: pyramid
      };
    };



    template< GeometryType :: BasicType type, unsigned int dim >
    struct Convert;

    template< unsigned int dim >
    struct Convert< GeometryType :: simplex, dim >
    {
      typedef Pyramid
      < typename Convert< GeometryType :: simplex, dim-1 > :: type >
      type;
    };

    template<>
    struct Convert< GeometryType :: simplex, 0 >
    {
      typedef Point type;
    };

    template< unsigned int dim >
    struct Convert< GeometryType :: cube, dim >
    {
      typedef Prism< typename Convert< GeometryType :: cube, dim-1 > :: type >
      type;
    };

    template<>
    struct Convert< GeometryType :: cube, 0 >
    {
      typedef Point type;
    };

    template< unsigned int dim >
    struct Convert< GeometryType :: prism, dim >
    {
      typedef Prism
      < typename Convert< GeometryType :: simplex, dim-1 > :: type >
      type;

    private:
      dune_static_assert( dim >= 3, "Dune prisms must be at least 3-dimensional." );
    };

    template< unsigned int dim >
    struct Convert< GeometryType :: pyramid, dim >
    {
      typedef Pyramid
      < typename Convert< GeometryType :: cube, dim-1 > :: type >
      type;

    private:
      dune_static_assert( dim >= 3, "Dune pyramids must be at least 3-dimensional." );
    };

  }

}

#endif
