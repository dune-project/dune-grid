// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_SUBTOPOLOGIES_HH
#define DUNE_GENERICGEOMETRY_SUBTOPOLOGIES_HH

#include <vector>

#include <dune/common/static_assert.hh>

#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/genericgeometry/codimtable.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< class Topology, unsigned int codim >
    struct Size;

    template< class Topology, unsigned int codim, unsigned int i >
    struct SubTopology;

    template< class Topology, unsigned int codim, unsigned int i,
        unsigned int subcodim, unsigned int j >
    class SubTopologyNumber;

    template< class Topology >
    class SubTopologyNumbering;



    // Size
    // ----

    template< class Topology, unsigned int dim, unsigned int codim >
    class SizeImpl;

    template< unsigned int dim, unsigned int codim >
    class SizeImpl< Point, dim, codim >
    {
      typedef Point Topology;
      dune_static_assert( (dim == Topology :: dimension), "Wrong dimension" );
      dune_static_assert( (codim <= dim), "Invalid codimension" );

    public:
      enum { value = 1 };
    };

    template< class BaseTopology, unsigned int dim, unsigned int codim >
    class SizeImpl< Prism< BaseTopology >, dim, codim >
    {
      typedef Prism< BaseTopology > Topology;
      dune_static_assert( (dim == Topology :: dimension), "Wrong dimension" );
      dune_static_assert( (codim <= dim), "Invalid codimension" );

      enum { m = Size< BaseTopology, codim-1 > :: value };
      enum { n = Size< BaseTopology, codim > :: value };

    public:
      enum { value = n + 2*m };
    };

    template< class BaseTopology, unsigned int dim >
    class SizeImpl< Prism< BaseTopology >, dim, 0 >
    {
      typedef Prism< BaseTopology > Topology;
      dune_static_assert( (dim == Topology :: dimension), "Wrong dimension" );

    public:
      enum { value = 1 };
    };

    template< class BaseTopology, unsigned int dim >
    class SizeImpl< Prism< BaseTopology >, dim, dim >
    {
      typedef Prism< BaseTopology > Topology;
      dune_static_assert( (dim == Topology :: dimension), "Wrong dimension" );

      enum { m = Size< BaseTopology, dim-1 > :: value };

    public:
      enum { value = 2*m };
    };

    template< class BaseTopology, unsigned int dim, unsigned int codim >
    struct SizeImpl< Pyramid< BaseTopology >, dim, codim >
    {
      typedef Pyramid< BaseTopology > Topology;
      dune_static_assert( (dim == Topology :: dimension), "Wrong dimension" );
      dune_static_assert( (codim <= dim), "Invalid codimension" );

      enum { m = Size< BaseTopology, codim-1 > :: value };
      enum { n = Size< BaseTopology, codim > :: value };

    public:
      enum { value = m+n };
    };

    template< class BaseTopology, unsigned int dim >
    class SizeImpl< Pyramid< BaseTopology >, dim, 0 >
    {
      typedef Pyramid< BaseTopology > Topology;
      dune_static_assert( (dim == Topology :: dimension), "Wrong dimension" );

    public:
      enum { value = 1 };
    };

    template< class BaseTopology, unsigned int dim >
    class SizeImpl< Pyramid< BaseTopology >, dim, dim >
    {
      typedef Pyramid< BaseTopology > Topology;
      dune_static_assert( (dim == Topology :: dimension), "Wrong dimension" );

      enum { m = Size< BaseTopology, dim-1 > :: value };

    public:
      enum { value = m+1 };
    };


    template< class Topology, unsigned int codim >
    struct Size
    {
      enum { value = SizeImpl< Topology, Topology :: dimension, codim > :: value };
    };



    // SubTopology
    // -----------

    template< class Topology, unsigned int dim, unsigned int codim, unsigned int i >
    class SubTopologyImpl;

    template< unsigned int dim, unsigned int codim, unsigned int i >
    class SubTopologyImpl< Point, dim, codim, i >
    {
      typedef Point Topology;
      dune_static_assert( (dim == Topology :: dimension), "Wrong dimension" );
      dune_static_assert( (codim <= dim), "Invalid codimension" );
      dune_static_assert( (i < Size< Topology, codim > :: value),
                          "Invalid subentity index" );

    public:
      typedef Topology type;
    };

    template< class BaseTopology, unsigned int dim, unsigned int codim, unsigned int i >
    class SubTopologyImpl< Prism< BaseTopology >, dim, codim, i >
    {
      typedef Prism< BaseTopology > Topology;
      dune_static_assert( (dim == Topology :: dimension), "Wrong dimension" );
      dune_static_assert( (codim <= dim), "Invalid codimension" );
      dune_static_assert( (i < Size< Topology, codim > :: value),
                          "Invalid subentity index" );

      enum { m = Size< BaseTopology, codim-1 > :: value };
      enum { n = Size< BaseTopology, codim > :: value };

      enum { s = (i < n+m ? 0 : 1) };

      template< bool >
      struct PrismSub
      {
        typedef Prism< typename SubTopology< BaseTopology, codim, i > :: type > type;
      };

      template< bool >
      struct BaseSub
      {
        typedef typename SubTopology< BaseTopology, codim-1, i-(n+s*m) > :: type type;
      };

    public:
      typedef typename ProtectedIf< (i < n), PrismSub, BaseSub > :: type type;
    };

    template< class BaseTopology, unsigned int dim, unsigned int i >
    class SubTopologyImpl< Prism< BaseTopology >, dim, 0, i >
    {
      typedef Prism< BaseTopology > Topology;
      dune_static_assert( (dim == Topology :: dimension), "Wrong dimension" );
      dune_static_assert( (i < Size< Topology, 0 > :: value),
                          "Invalid subentity index" );
    public:
      typedef Topology type;
    };

    template< class BaseTopology, unsigned int dim, unsigned int i >
    class SubTopologyImpl< Prism< BaseTopology >, dim, dim, i >
    {
      typedef Prism< BaseTopology > Topology;
      dune_static_assert( (dim == Topology :: dimension), "Wrong dimension" );
      dune_static_assert( (i < Size< Topology, dim > :: value),
                          "Invalid subentity index" );
    public:
      typedef Point type;
    };

    template< class BaseTopology, unsigned int dim, unsigned int codim, unsigned int i >
    class SubTopologyImpl< Pyramid< BaseTopology >, dim, codim, i >
    {
      typedef Pyramid< BaseTopology > Topology;
      dune_static_assert( (dim == Topology :: dimension), "Wrong dimension" );
      dune_static_assert( (codim <= dim), "Invalid codimension" );
      dune_static_assert( (i < Size< Topology, codim > :: value),
                          "Invalid subentity index" );

      enum { m = Size< BaseTopology, codim-1 > :: value };

      template< bool >
      struct BaseSub
      {
        typedef typename SubTopology< BaseTopology, codim-1, i > :: type type;
      };

      template< bool >
      struct PyramidSub
      {
        typedef Pyramid< typename SubTopology< BaseTopology, codim, i-m > :: type > type;
      };

    public:
      typedef typename ProtectedIf< (i < m), BaseSub, PyramidSub > :: type type;
    };

    template< class BaseTopology, unsigned int dim, unsigned int i >
    class SubTopologyImpl< Pyramid< BaseTopology >, dim, 0, i >
    {
      typedef Pyramid< BaseTopology > Topology;
      dune_static_assert( (dim == Topology :: dimension), "Wrong dimension" );
      dune_static_assert( (i < Size< Topology, 0 > :: value),
                          "Invalid subentity index" );

    public:
      typedef Topology type;
    };

    template< class BaseTopology, unsigned int dim, unsigned int i >
    class SubTopologyImpl< Pyramid< BaseTopology >, dim, dim, i >
    {
      typedef Pyramid< BaseTopology > Topology;
      dune_static_assert( (dim == Topology :: dimension), "Wrong dimension" );
      dune_static_assert( (i < Size< Topology, dim > :: value),
                          "Invalid subentity index" );

    public:
      typedef Point type;
    };

    template< class Topology, unsigned int codim, unsigned int i >
    struct SubTopology
    {
      typedef typename SubTopologyImpl< Topology, Topology :: dimension, codim, i > :: type type;
    };



    // SubTopologyNumber
    // -----------------

    template< class Topology, unsigned int codim, unsigned int i,
        unsigned int subdim, unsigned int subcodim, unsigned int j >
    struct SubTopologyNumberHelper;

    template< class BaseTopology, unsigned int codim, unsigned int i,
        unsigned int subdim, unsigned int subcodim, unsigned int j >
    struct SubTopologyNumberHelper< Prism< BaseTopology >, codim, i, subdim, subcodim, j >
    {
      typedef Prism< BaseTopology > Topology;
      typedef typename GenericGeometry :: SubTopology< Topology, codim, i > :: type
      SubTopology;

      enum { m = Size< BaseTopology, codim-1 > :: value };
      enum { n = Size< BaseTopology, codim > :: value };

      enum { mb = Size< BaseTopology, codim+subcodim-1 > :: value };
      enum { nb = Size< BaseTopology, codim+subcodim > :: value };

      enum { s = (i < n+m ? 0 : 1) };

      template< bool >
      struct PrismSub
      {
        typedef typename GenericGeometry :: BaseTopology< SubTopology > :: type
        SubBaseTopology;

        enum { ms = Size< SubBaseTopology, subcodim-1 > :: value };
        enum { ns = Size< SubBaseTopology, subcodim > :: value };

        enum { s = (j < ns+ms ? 0 : 1) };

        template< bool >
        struct PrismSubSub
        {
          enum { value = SubTopologyNumber< BaseTopology, codim, i, subcodim, j > :: value };
        };

        template< bool >
        struct BaseSubSub
        {
          enum
          {
            value = SubTopologyNumber< BaseTopology, codim, i, subcodim-1, j-(ns+s*ms) > :: value
                    + nb + s*mb
          };
        };

        enum { value = ProtectedIf< (j < ns), PrismSubSub, BaseSubSub > :: value };
      };

      template< bool >
      struct BaseSub
      {
        enum
        {
          value = SubTopologyNumber< BaseTopology, codim-1, i-(n+s*m), subcodim, j > :: value
                  + nb + s*mb
        };
      };

      enum { value = ProtectedIf< (i < n), PrismSub, BaseSub > :: value };
    };

    template< class BaseTopology, unsigned int codim, unsigned int i,
        unsigned int subdim, unsigned int j >
    struct SubTopologyNumberHelper< Prism< BaseTopology >, codim, i, subdim, 0, j >
    {
      enum { value = i };
    };

    template< class BaseTopology, unsigned int codim, unsigned int i,
        unsigned int subdim, unsigned int j >
    struct SubTopologyNumberHelper< Prism< BaseTopology >, codim, i, subdim, subdim, j >
    {
      typedef Prism< BaseTopology > Topology;

      typedef typename GenericGeometry :: SubTopology< Topology, codim, i > :: type
      SubTopology;

      enum { m = Size< BaseTopology, codim-1 > :: value };
      enum { n = Size< BaseTopology, codim > :: value };

      enum { mb = Size< BaseTopology, codim+subdim-1 > :: value };

      enum { s = (i < n+m ? 0 : 1) };

      template< bool >
      struct PrismSub
      {
        typedef typename GenericGeometry :: BaseTopology< SubTopology > :: type
        SubBaseTopology;

        enum { ms = Size< SubBaseTopology, subdim-1 > :: value };

        enum { s = (j < ms ? 0 : 1) };

        enum
        {
          value = SubTopologyNumber< BaseTopology, codim, i, subdim-1, j-s*ms > :: value
                  + s*mb
        };
      };

      template< bool >
      struct BaseSub
      {
        enum
        {
          value = SubTopologyNumber< BaseTopology, codim-1, i-(n+s*m), subdim, j > :: value
                  + s*mb
        };
      };

      enum { value = ProtectedIf< (i < n), PrismSub, BaseSub > :: value };
    };

    template< class BaseTopology, unsigned int codim, unsigned int i,
        unsigned int subdim, unsigned int subcodim, unsigned int j >
    struct SubTopologyNumberHelper< Pyramid< BaseTopology >, codim, i, subdim, subcodim, j >
    {
      typedef Pyramid< BaseTopology > Topology;

      typedef typename GenericGeometry :: SubTopology< Topology, codim, i > :: type
      SubTopology;

      enum { m = Size< BaseTopology, codim-1 > :: value };

      enum { mb = Size< BaseTopology, codim+subcodim-1 > :: value };

      template< bool >
      struct BaseSub
      {
        enum
        {
          value = SubTopologyNumber< BaseTopology, codim-1, i, subcodim, j > :: value
        };
      };

      template< bool >
      struct PyramidSub
      {
        typedef typename GenericGeometry :: BaseTopology< SubTopology > :: type
        SubBaseTopology;

        enum { ms = Size< SubBaseTopology, subcodim-1 > :: value };

        template< bool >
        struct BaseSubSub
        {
          enum { value = SubTopologyNumber< BaseTopology, codim, i-m, subcodim-1, j > :: value };
        };

        template< bool >
        struct PyramidSubSub
        {
          enum
          {
            value = SubTopologyNumber< BaseTopology, codim, i-m, subcodim, j-ms > :: value
                    + mb
          };
        };

        enum { value = ProtectedIf< (j < ms), BaseSubSub, PyramidSubSub > :: value };
      };

      enum { value = ProtectedIf< (i < m), BaseSub, PyramidSub > :: value };
    };

    template< class BaseTopology, unsigned int codim, unsigned int i,
        unsigned int subdim, unsigned int j >
    struct SubTopologyNumberHelper< Pyramid< BaseTopology >, codim, i, subdim, 0, j >
    {
      enum { value = i };
    };

    template< class BaseTopology, unsigned int codim, unsigned int i,
        unsigned int subdim, unsigned int j >
    struct SubTopologyNumberHelper< Pyramid< BaseTopology >, codim, i, subdim, subdim, j >
    {
      typedef Pyramid< BaseTopology > Topology;

      typedef typename GenericGeometry :: SubTopology< Topology, codim, i > :: type
      SubTopology;

      enum { m = Size< BaseTopology, codim-1 > :: value };

      enum { mb = Size< BaseTopology, codim+subdim-1 > :: value };

      template< bool >
      struct BaseSub
      {
        enum
        {
          value = SubTopologyNumber< BaseTopology, codim-1, i, subdim, j > :: value
        };
      };

      template< bool >
      struct PyramidSub
      {
        typedef typename GenericGeometry :: BaseTopology< SubTopology > :: type
        SubBaseTopology;

        enum { ms = Size< SubBaseTopology, subdim-1 > :: value };

        template< bool >
        struct BaseSubSub
        {
          enum { value = SubTopologyNumber< BaseTopology, codim, i-m, subdim-1, j > :: value };
        };

        template< bool >
        struct PyramidSubSub
        {
          enum { value = mb };
        };

        enum { value = ProtectedIf< (j < ms), BaseSubSub, PyramidSubSub > :: value };
      };

      enum { value = ProtectedIf< (i < m), BaseSub, PyramidSub > :: value };
    };


    template< class Topology, unsigned int codim, unsigned int i,
        unsigned int subcodim, unsigned int j >
    class SubTopologyNumber
    {
      dune_static_assert( (codim <= Topology :: dimension), "Invald codimension" );
      dune_static_assert( (i < Size< Topology, codim > :: value),
                          "Invalid subentity index" );

      typedef typename GenericGeometry :: SubTopology< Topology, codim, i > :: type
      SubTopology;
      dune_static_assert( (subcodim <= SubTopology :: dimension),
                          "Invald subcodimension" );
      dune_static_assert( (j < Size< SubTopology, subcodim > :: value),
                          "Invalid subsubentity index" );

      template< bool >
      struct Border
      {
        enum { value = (codim == 0 ? j : i) };
      };

      template< bool >
      struct Interior
      {
        enum
        {
          value = SubTopologyNumberHelper< Topology, codim, i, SubTopology :: dimension,
              subcodim, j > :: value
        };
      };

    public:
      enum
      {
        value = ProtectedIf< (codim == 0) || (codim == Topology :: dimension),
            Border, Interior > :: value
      };
    };



    // SubTopologyNumbering
    // --------------------

    template< class Topology >
    class SubTopologyNumbering
    {
      template< int codim >
      class Numbering;

      template< int codim >
      class CodimNumbering
        : public std :: vector< Numbering< codim > >
      {};

      template< int codim >
      struct Builder;

      CodimTable< CodimNumbering, Topology :: dimension > codimNumbering_;

    public:
      template< unsigned int codim, unsigned int subcodim >
      static unsigned int subEntity ( unsigned int i, unsigned int j )
      {
        Int2Type< codim > codimVariable;
        const CodimNumbering< codim > &numbering
          = instance().codimNumbering_[ codimVariable ];
        return numbering[ i ].template number< subcodim >( j );
      }

    private:
      SubTopologyNumbering ()
      {
        ForLoop< Builder, 0, Topology :: dimension > :: apply( *this );
      }

      static const SubTopologyNumbering &instance ()
      {
        static SubTopologyNumbering inst;
        return inst;
      }
    };



    template< class Topology >
    template< int codim >
    class SubTopologyNumbering< Topology > :: Numbering
    {
      std :: vector< unsigned int > numbering_[ Topology :: dimension - codim + 1 ];

    public:
      template< int subcodim >
      const unsigned int &number ( unsigned int j ) const
      {
        return numbering_[ subcodim ][ j ];
      }

      template< int subcodim >
      unsigned int &number ( unsigned int j )
      {
        return numbering_[ subcodim ][ j ];
      }

      template< int subcodim >
      void resize ( unsigned int size )
      {
        numbering_[ subcodim ].resize( size );
      }
    };



    template< class Topology >
    template< int codim >
    struct SubTopologyNumbering< Topology > :: Builder
    {
      typedef GenericGeometry :: SubTopologyNumbering< Topology > SubTopologyNumbering;
      typedef typename SubTopologyNumbering :: template CodimNumbering< codim >
      CodimNumbering;
      typedef typename SubTopologyNumbering :: template Numbering< codim > Numbering;

      template< int i >
      struct Sub
      {
        typedef typename GenericGeometry :: SubTopology< Topology, codim, i > :: type
        SubTopology;

        template< int subcodim >
        struct SubCodim
        {
          template< int j >
          struct SubSub
          {
            static void apply( Numbering &numbering )
            {
              numbering.template number< subcodim >( j )
                = SubTopologyNumber< Topology, codim, i, subcodim, j > :: value;
            }
          };

          static void apply ( Numbering &numbering )
          {
            enum { numSubSubs = Size< SubTopology, subcodim > :: value };
            numbering.template resize< subcodim >( numSubSubs );
            ForLoop< SubSub, 0, numSubSubs-1 > :: apply( numbering );
          }
        };

        static void apply ( CodimNumbering &numbering )
        {
          ForLoop< SubCodim, 0, SubTopology :: dimension > :: apply( numbering[ i ] );
        }
      };

      static void apply ( SubTopologyNumbering &numbering )
      {
        enum { numSubs = Size< Topology, codim > :: value };
        Int2Type< codim > codimVariable;
        CodimNumbering &codimNumbering = numbering.codimNumbering_[ codimVariable ];
        codimNumbering.resize( numSubs );
        ForLoop< Sub, 0, numSubs-1 > :: apply( codimNumbering );
      }
    };

  }

}

#endif
