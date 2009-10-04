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
        return BaseTopology :: name() + "l";
        // return BaseTopology :: name() + "'";
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
        return BaseTopology :: name() + "o";
        // return BaseTopology :: name() + "°";
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

    template< class Topology >
    struct IsGeneralizedPrism
    {
      static const bool value = false;
    };

    template< class BaseTopology >
    struct IsGeneralizedPrism< Prism< BaseTopology > >
    {
      static const bool value
        = (IsGeneralizedPrism< BaseTopology >::value || IsSimplex< BaseTopology >::value);
    };



    // SimplexTopology
    // ---------------

    template< unsigned int dim >
    struct SimplexTopology
    {
      typedef Pyramid< typename SimplexTopology< dim-1 >::type > type;
    };

    template<>
    struct SimplexTopology< 0 >
    {
      typedef Point type;
    };



    // CubeTopology
    // ------------

    template< unsigned int dim >
    struct CubeTopology
    {
      typedef Prism< typename CubeTopology< dim-1 >::type > type;
    };

    template<>
    struct CubeTopology< 0 >
    {
      typedef Point type;
    };



    // PyramidTopoogy
    // --------------

    template< unsigned int dim >
    struct PyramidTopology
    {
      typedef Pyramid< typename CubeTopology< dim-1 >::type > type;
    };



    // PrismTopoogy
    // --------------

    template< unsigned int dim >
    struct PrismTopology
    {
      typedef Prism< typename SimplexTopology< dim-1 >::type > type;
    };



    // Topology
    // --------

    template< unsigned int id, unsigned int dim >
    class Topology
    {
      static const unsigned int dimension = dim;

      dune_static_assert( (id < (1 << dimension)), "id too large." );

      static const bool isPrism = ((id >> (dimension-1)) != 0);

      typedef typename Topology< (id & ~(1 << (dimension-1))), dimension-1 >::type
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
      typedef typename SelectType< isPrism, Prism<true>, Pyramid<false> >::Type::type type;
    };

    template< unsigned int id >
    class Topology< id, 0 >
    {
      static const unsigned int dimension = 0;

      dune_static_assert( (id < (1 << dimension)), "id too large." );

    public:
      typedef Point type;
    };



    // IfTopology
    // ----------

    template< template< class > class Operation, int dim, class Topology = Point >
    class IfTopology
    {
      typedef IfTopology< Operation, dim-1, Prism< Topology > > IfPrism;
      typedef IfTopology< Operation, dim-1, Pyramid< Topology > > IfPyramid;

    public:
      static void apply ( const unsigned int topologyId )
      {
        if( topologyId & 1 )
          IfPrism::apply( topologyId >> 1 );
        else
          IfPyramid::apply( topologyId >> 1 );
      }

      template< class T1 >
      static void apply ( const unsigned int topologyId, T1 &p1 )
      {
        if( topologyId & 1 )
          IfPrism::apply( topologyId >> 1, p1 );
        else
          IfPyramid::apply( topologyId >> 1, p1 );
      }

      template< class T1, class T2 >
      static void apply ( const unsigned int topologyId, T1 &p1, T2 &p2 )
      {
        if( topologyId & 1 )
          IfPrism::apply( topologyId >> 1, p1, p2 );
        else
          IfPyramid::apply( topologyId >> 1, p1, p2 );
      }

      template< class T1, class T2, class T3 >
      static void apply ( const unsigned int topologyId, T1 &p1, T2 &p2, T3 &p3 )
      {
        if( topologyId & 1 )
          IfPrism::apply( topologyId >> 1, p1, p2, p3 );
        else
          IfPyramid::apply( topologyId >> 1, p1, p2, p3 );
      }

      template< class T1, class T2, class T3, class T4 >
      static void apply ( const unsigned int topologyId, T1 &p1, T2 &p2, T3 &p3, T4 &p4 )
      {
        if( topologyId & 1 )
          IfPrism::apply( topologyId >> 1, p1, p2, p3, p4 );
        else
          IfPyramid::apply( topologyId >> 1, p1, p2, p3, p4 );
      }
    };

    template< template< class > class Operation, class Topology >
    struct IfTopology< Operation, 0, Topology >
    {
      static void apply ( const unsigned int topologyId )
      {
        Operation< Topology >::apply();
      }

      template< class T1 >
      static void apply ( const unsigned int topologyId, T1 &p1 )
      {
        Operation< Topology >::apply( p1 );
      }

      template< class T1, class T2 >
      static void apply ( const unsigned int topologyId, T1 &p1, T2 &p2 )
      {
        Operation< Topology >::apply( p1, p2 );
      }

      template< class T1, class T2, class T3 >
      static void apply ( const unsigned int topologyId, T1 &p1, T2 &p2, T3 &p3 )
      {
        Operation< Topology >::apply( p1, p2, p3 );
      }

      template< class T1, class T2, class T3, class T4 >
      static void apply ( const unsigned int topologyId, T1 &p1, T2 &p2, T3 &p3, T4 &p4 )
      {
        Operation< Topology >::apply( p1, p2, p3, p4 );
      }
    };

  }

}

#endif
