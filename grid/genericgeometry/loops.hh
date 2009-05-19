// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
namespace Dune
{
  namespace GenericGeometry
  {

    template< template< class,unsigned int > class Operation, class Topology, unsigned int dimension, unsigned int codimension >
    class ForSubTopologies;

    template< template< class,unsigned int > class Operation, class Topology >
    struct ForSubTopology
    {
      static void apply ( )
      {
        ForLoop< Helper, 0, Topology::dimension >::apply( );
      }
    private:
      template <int dim>
      struct Helper
      {
        static void apply ( )
        {
          ForSubTopologies< Operation, Topology, dim, Topology::dimension-dim >::template apply<0>( );
        }
      };
    };

    // ForSubTopologies on Pyramids:
    // --------------------

    // codimension = dimension (start)
    template< template< class,unsigned int > class Operation, class BaseTopology, unsigned int codimension  >
    class ForSubTopologies< Operation, Pyramid< BaseTopology>, 0, codimension >
    {
      typedef Pyramid< BaseTopology > Topology;
      static const unsigned int dimension = 0;
    public:
      template <unsigned int start>
      static void apply ( )
      {
        // over all points in base
        ForSubTopologies< Operation, BaseTopology, dimension, codimension-1 >::template apply<start>( );
        // call on the top point
        const unsigned int numVtx = Size< BaseTopology, codimension-1 >::value;
        Operation< Point,start+numVtx >::apply( );
      }
      template < template <class B> class T,unsigned int start, unsigned int size >
      static void apply ( )
      {
        ForLoop< ForSubTopologies< Operation,Point,dimension,0 >::template
            Helper< T >::template Apply,start,start+size-1 >::apply( );
      }
    };

    // codimension not equal to 0 or dimension
    template< template< class,unsigned int > class Operation, class BaseTopology, unsigned int dimension, unsigned int codimension >
    class ForSubTopologies< Operation, Pyramid< BaseTopology >, dimension, codimension >
    {
      typedef Pyramid< BaseTopology > Topology;
    public:
      template <unsigned int start>
      static void apply ( )
      {
        // first call same codimension on base
        ForSubTopologies< Operation, BaseTopology, dimension, codimension-1 >::template apply<start>( );
        // now construct upwards going topologies
        const unsigned int nextStart = start+Size< BaseTopology, codimension-1 >::value;
        const unsigned int numSub = Size< BaseTopology, codimension >::value;
        ForSubTopologies< Operation, BaseTopology, dimension-1, codimension >::template apply<HelperBase::template Apply,nextStart,numSub>( );
      }
      template < template <class B> class T,unsigned int start,unsigned int size >
      static void apply ( )
      {
        const unsigned int nextStart = start;
        const unsigned int nextSize = size;
        ForSubTopologies< Operation, BaseTopology,dimension-1,codimension >::template apply< Helper<T>::template Apply,nextStart,nextSize >( );
      }
    private:
      struct HelperBase
      {
        template <class T>
        struct Apply
        {
          typedef Pyramid<T> type;
        };
      };
      template <template <class> class TopT>
      struct Helper
      {
        template <class T>
        struct Apply
        {
          typedef typename TopT< Pyramid<T> >::type type;
        };
      };
    };

    // codimension 0 (final)
    template< template< class,unsigned int > class Operation, class BaseTopology, unsigned int dimension >
    class ForSubTopologies< Operation, Pyramid<BaseTopology>, dimension, 0 >
    {
      typedef Pyramid< BaseTopology > Topology;
    public:
      template <unsigned int start>
      static void apply ( )
      {
        Operation< Topology,start > :: apply ( );
      }
    };


    // ForSubTopologies on Prism:
    // --------------------

    // codimension = dimension (start)
    template< template< class,unsigned int > class Operation, class BaseTopology, unsigned int codimension  >
    class ForSubTopologies< Operation, Prism< BaseTopology>, 0, codimension >
    {
      typedef Prism< BaseTopology > Topology;
      static const unsigned int dimension = 0;
    public:
      template <unsigned int start>
      static void apply ( )
      {
        // over all points in base
        ForSubTopologies< Operation, BaseTopology, dimension, codimension-1 >::template apply<start>( );
        // call on all top point
        const unsigned int nextStart = start+Size< BaseTopology, codimension-1 >::value;
        ForSubTopologies< Operation, BaseTopology, dimension, codimension-1 >::template apply<nextStart>( );
      }
      template < template <class B> class T,unsigned int start, unsigned int size >
      static void apply ( )
      {
        ForLoop< ForSubTopologies< Operation,Point,dimension,0 >::template
            Helper< T >::template Apply,start,start+size-1 >::apply( );
      }
    };

    // codimension not equal to 0 or dimension
    template< template< class,unsigned int > class Operation, class BaseTopology, unsigned int dimension, unsigned int codimension >
    class ForSubTopologies< Operation, Prism< BaseTopology >, dimension, codimension >
    {
      typedef Prism< BaseTopology > Topology;
    public:
      template <unsigned int start>
      static void apply ( )
      {
        // firstconstruct upwards going topologies
        const unsigned int numSub = Size< BaseTopology, codimension >::value;
        ForSubTopologies< Operation, BaseTopology, dimension-1, codimension >::template apply<HelperBase::template Apply,start,numSub>( );
        // now call same codimension on bottom base
        const unsigned int baseSize = Size< BaseTopology, codimension-1 >::value;
        ForSubTopologies< Operation, BaseTopology, dimension, codimension-1 >::template apply<start+numSub>( );
        // now call same codimension on top base
        ForSubTopologies< Operation, BaseTopology, dimension, codimension-1 >::template apply<start+numSub+baseSize>( );
      }
      template < template <class B> class T,unsigned int start,unsigned int size >
      static void apply ( )
      {
        const unsigned int nextStart = start;
        const unsigned int nextSize = size;
        ForSubTopologies< Operation, BaseTopology,dimension-1,codimension >::template apply< Helper<T>::template Apply,nextStart,nextSize >( );
      }
    private:
      struct HelperBase
      {
        template <class T>
        struct Apply
        {
          typedef Prism<T> type;
        };
      };
      template <template <class> class TopT>
      struct Helper
      {
        template <class T>
        struct Apply
        {
          typedef typename TopT< Prism<T> >::type type;
        };
      };
    };

    // codimension 0 (final)
    template< template< class,unsigned int > class Operation, class BaseTopology, unsigned int dimension >
    class ForSubTopologies< Operation, Prism<BaseTopology>, dimension, 0 >
    {
      typedef Prism< BaseTopology > Topology;
    public:
      template <unsigned int start>
      static void apply ( )
      {
        Operation< Topology,start > :: apply ( );
      }
    };


    // ForSubTopologies on Points:
    // --------------------
    template< template< class,unsigned int > class Operation >
    class ForSubTopologies< Operation, Point, 0, 0 >
    {
      typedef Point Topology;
    public:
      template <unsigned int start>
      static void apply ( )
      {
        Operation< Topology,start > :: apply ( );
      }
      template < template <class B> class T,unsigned int start,unsigned int size >
      static void apply ( )
      {
        ForLoop<Helper< T >::template Apply,start,start+size-1>::apply( );
      }
      template <template <class> class T>
      struct Helper
      {
        template <int i>
        struct Apply
        {
          static void apply ( )
          {
            Operation<typename T<Point>::type,i>::apply( );
          }
        };
      };
    };

  }
}
