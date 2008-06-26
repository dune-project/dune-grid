// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_GEOMETRYTYPES_HH
#define DUNE_GENERICGEOMETRY_GEOMETRYTYPES_HH

#include <dune/common/static_assert.hh>
#include <dune/common/geometrytype.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< GeometryType :: BasicType basicType >
    class Point
    {
      dune_static_assert( (basicType == GeometryType :: simplex)
                          || (basicType == GeometryType :: cube),
                          "A point can only be a simplex or a cube." );

    protected:
      enum { basicDuneType = basicType };

    public:
      enum { dimension = 0 };
      enum { duneType = basicDuneType };
    };


    template< class BaseGeometry >
    class Prism
    {
      dune_static_assert( (BaseGeometry :: duneType == GeometryType :: simplex)
                          || (BaseGeometry :: duneType == GeometryType :: cube),
                          "Prisms can only be generated over simplices or cubes." );

    protected:
      enum { basicDuneType = BaseGeometry :: basicDuneType };


    public:
      enum { dimension = BaseGeometry :: dimension + 1 };
      enum { duneType = (dimension == 1)
                        ? basicDuneType
                        : ((dimension == 2)
                           || (BaseGeometry :: duneType == GeometryType :: cube))
                        ? GeometryType :: cube
                        : GeometryType :: prism };
    };


    template< class BaseGeometry >
    class Pyramid
    {
      dune_static_assert( (BaseGeometry :: duneType == GeometryType :: simplex)
                          || (BaseGeometry :: duneType == GeometryType :: cube),
                          "Pyramids can only be generated over simplices or cubes." );

    protected:
      enum { basicDuneType = BaseGeometry :: basicDuneType };

    public:
      enum { dimension = BaseGeometry :: dimension + 1 };
      enum { duneType = (dimension == 1)
                        ? basicDuneType
                        : ((dimension == 2)
                           || (BaseGeometry :: duneType == GeometryType :: simplex))
                        ? GeometryType :: simplex
                        : GeometryType :: pyramid };
    };



    template< GeometryType :: BasicType type, unsigned int dim,
        GeometryType :: BasicType basicType = GeometryType :: simplex >
    struct Convert;

    template< unsigned int dim, GeometryType :: BasicType basicType >
    struct Convert< GeometryType :: simplex, dim, basicType >
    {
    private:
      typedef typename Convert< GeometryType :: simplex, dim-1, basicType > :: type
      base_type;

    public:
      typedef Pyramid< base_type > type;
    };

    template< GeometryType :: BasicType basicType >
    struct Convert< GeometryType :: simplex, 0, basicType >
    {
      typedef Point< basicType > type;
    };

    template< unsigned int dim, GeometryType :: BasicType basicType >
    struct Convert< GeometryType :: cube, dim, basicType >
    {
    private:
      typedef typename Convert< GeometryType :: cube, dim-1, basicType > :: type
      base_type;

    public:
      typedef Prism< base_type > type;
    };

    template< GeometryType :: BasicType basicType >
    struct Convert< GeometryType :: cube, 0, basicType >
    {
      typedef Point< basicType > type;
    };

    template< unsigned int dim, GeometryType :: BasicType basicType >
    struct Convert< GeometryType :: prism, dim, basicType >
    {
    private:
      dune_static_assert( dim >= 3, "Prisms must be at least 3-dimensional." );

      typedef typename Convert< GeometryType :: simplex, dim-1, basicType > :: type
      base_type;

    public:
      typedef Prism< base_type > type;
    };

    template< unsigned int dim, GeometryType :: BasicType basicType >
    struct Convert< GeometryType :: pyramid, dim, basicType >
    {
    private:
      dune_static_assert( dim >= 3, "Pyramids must be at least 3-dimensional." );

      typedef typename Convert< GeometryType :: cube, dim-1, basicType > :: type
      base_type;

    public:
      typedef Pyramid< base_type > type;
    };

  }

}

#endif
