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

      enum { id = 0 };

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

      enum { id = BaseTopology :: id + (1 << (dimension-1)) };

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

      enum { id = BaseTopology :: id };

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



    template< class Topology >
    struct IsHybrid
    {
      // Only cubes and simplices are nonhybrid
      enum { value = ((Topology :: id | 1) == (1 << Topology :: dimension) - 1)
                     || ((Topology :: id >> 1) == 0) };
    };



    template< unsigned int id, unsigned int dim >
    class Topology
    {
      enum { dimension = dim };

      dune_static_assert( (id < (1 << dimension)), "id too large." );

      enum { isPrism = ((id >> (dimension-1)) != 0) };
      typedef typename Topology< (id & ~(1 << (dimension-1))), dimension-1 > :: type
      BaseTopology;

      template< bool >
      struct Prism
      {
        typedef GenericGeometry :: Prism< BaseTopology > type;
      };

      template< bool >
      struct Pyramid
      {
        typedef GenericGeometry :: Pyramid< BaseTopology > type;
      };

    public:
      typedef typename ProtectedIf< isPrism, Prism, Pyramid > :: type type;
    };

    template< unsigned int id >
    class Topology< id, 0 >
    {
      enum { dimension = 0 };

      dune_static_assert( (id < (1 << dimension)), "id too large." );

    public:
      typedef Point type;
    };

  }

}

#endif
