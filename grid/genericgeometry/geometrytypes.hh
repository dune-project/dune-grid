// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_GEOMETRYTYPES_HH
#define DUNE_GENERICGEOMETRY_GEOMETRYTYPES_HH

#include <string>

namespace Dune
{

  namespace GenericGeometry
  {

    struct Point
    {
      enum { dimension = 0 };
      enum { numCorners = 1 };

      static std :: string name ()
      {
        return "p";
      }
    };


    template< class BaseGeometry >
    struct Prism
    {
      enum { dimension = BaseGeometry :: dimension + 1 };
      enum { numCorners = 2 * BaseGeometry :: numCorners };

      static std :: string name ()
      {
        return BaseGeometry :: name() + "'";
      }
    };


    template< class BaseGeometry >
    struct Pyramid
    {
      enum { dimension = BaseGeometry :: dimension + 1 };
      enum { numCorners = BaseGeometry :: numCorners + 1 };

      static std :: string name ()
      {
        return BaseGeometry :: name() + "°";
      }
    };

  }

}

#endif
