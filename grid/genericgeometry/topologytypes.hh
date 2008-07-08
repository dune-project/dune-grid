// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_TOPOLOGYTYPES_HH
#define DUNE_GENERICGEOMETRY_TOPOLOGYTYPES_HH

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


    template< class BaseTopology >
    struct Prism
    {
      enum { dimension = BaseTopology :: dimension + 1 };
      enum { numCorners = 2 * BaseTopology :: numCorners };

      static std :: string name ()
      {
        return BaseTopology :: name() + "'";
      }
    };


    template< class BaseTopology >
    struct Pyramid
    {
      enum { dimension = BaseTopology :: dimension + 1 };
      enum { numCorners = BaseTopology :: numCorners + 1 };

      static std :: string name ()
      {
        return BaseTopology :: name() + "°";
      }
    };



    template< class Topology >
    struct BaseTopology;

    template< class Base >
    struct BaseTopology< Prism< Base > >
    {
      typedef Base type;
    };

    template< class Base >
    struct BaseTopology< Pyramid< Base > >
    {
      typedef Base type;
    };

  }

}

#endif
