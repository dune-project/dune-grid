// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_TOPOLOGYTYPES_HH
#define DUNE_GENERICGEOMETRY_TOPOLOGYTYPES_HH

#include <string>

#include <dune/common/static_assert.hh>
#include <dune/grid/genericgeometry/misc.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    struct Point
    {
      static const unsigned int dimension = 0;
      static const unsigned int numCorners = 1;

      static const unsigned int id = 0;

      static std :: string name ()
      {
        return "p";
      }
    };


    template< class BaseTopology >
    struct Prism
    {
      static const unsigned int dimension = BaseTopology :: dimension + 1;
      static const unsigned int numCorners = 2 * BaseTopology :: numCorners;

      static const unsigned int id = BaseTopology :: id + (1 << (dimension-1));

      static std :: string name ()
      {
        return BaseTopology :: name() + "'";
      }
    };


    template< class BaseTopology >
    struct Pyramid
    {
      static const unsigned int dimension = BaseTopology :: dimension + 1;
      static const unsigned int numCorners = BaseTopology :: numCorners + 1;

      static const unsigned int id = BaseTopology :: id;

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
    struct IsSimplex
    {
      static const bool value = ((Topology::id >> 1) == 0);
    };

    template< class Topology >
    struct IsCube
    {
      static const bool value = ((Topology::id | 1) == (1 << Topology::dimension) - 1);
    };

    template< class Topology >
    struct IsHybrid
    {
      static const bool value
        = !(IsSimplex< Topology >::value || IsCube< Topology >::value);
    };



    template< unsigned int id, unsigned int dim >
    class Topology
    {
      static const unsigned int dimension = dim;

      dune_static_assert( (id < (1 << dimension)), "id too large." );

      static const bool isPrism = ((id >> (dimension-1)) != 0);

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
      static const unsigned int dimension = 0;

      dune_static_assert( (id < (1 << dimension)), "id too large." );

    public:
      typedef Point type;
    };

  }

}

#endif
