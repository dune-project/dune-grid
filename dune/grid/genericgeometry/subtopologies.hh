// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_SUBTOPOLOGIES_HH
#define DUNE_GENERICGEOMETRY_SUBTOPOLOGIES_HH

#include <cassert>
#include <vector>

#include <dune/common/forloop.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/typetraits.hh>

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

    template< class Topology, unsigned int codim, unsigned int subcodim >
    class SubTopologySize;

    template< class Topology, unsigned int codim, unsigned int subcodim >
    class GenericSubTopologyNumbering;

    template< class Topology, unsigned int codim, unsigned int subcodim >
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
      typedef typename SelectType< (i < n), PrismSub<true>, BaseSub<false> > :: Type :: type type;
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
      typedef typename SelectType< (i < m), BaseSub<true>, PyramidSub<false> > :: Type :: type type;
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



    // SubTopologySize
    // ---------------

    template< class Topology, unsigned int codim, unsigned int subcodim >
    class SubTopologySize
    {
      template< int i >
      struct Builder;

      unsigned int size_[ Size< Topology, codim > :: value ];

      SubTopologySize ()
      {
        ForLoop< Builder, 0, Size< Topology, codim > :: value-1 >
        :: apply( *this );
      }

      SubTopologySize ( const SubTopologySize & );

      static const SubTopologySize &instance ()
      {
        static SubTopologySize inst;
        return inst;
      }

    public:
      static unsigned int size ( unsigned int i )
      {
        assert( (i < Size< Topology, codim > :: value) );
        return instance().size_[ i ];
      }
    };

    template< class Topology, unsigned int codim, unsigned int subcodim >
    template< int i >
    struct SubTopologySize< Topology, codim, subcodim > :: Builder
    {
      typedef GenericGeometry :: SubTopologySize< Topology, codim, subcodim >
      SubTopologySize;
      typedef typename GenericGeometry :: SubTopology< Topology, codim, i > :: type
      SubTopology;

      static void apply ( SubTopologySize &subTopologySize )
      {
        subTopologySize.size_[ i ] = Size< SubTopology, subcodim > :: value;
      }
    };



    // GenericSubTopologyNumbering
    // ---------------------------

    template< class Topology, unsigned int codim,
        unsigned int subdim, unsigned int subcodim >
    struct GenericSubTopologyNumberingHelper;

    template< class BaseTopology, unsigned int codim,
        unsigned int subdim, unsigned int subcodim >
    struct GenericSubTopologyNumberingHelper
    < Prism< BaseTopology >, codim, subdim, subcodim >
    {
      typedef Prism< BaseTopology > Topology;

      enum { m = Size< BaseTopology, codim-1 > :: value };
      enum { n = Size< BaseTopology, codim > :: value };

      enum { mb = Size< BaseTopology, codim+subcodim-1 > :: value };
      enum { nb = Size< BaseTopology, codim+subcodim > :: value };

      static unsigned int number ( unsigned int i, unsigned int j )
      {
        const unsigned int s = (i < n+m ? 0 : 1);
        if( i < n )
        {
          const unsigned int ms = SubTopologySize< BaseTopology, codim, subcodim-1 > :: size( i );
          const unsigned int ns = SubTopologySize< BaseTopology, codim, subcodim > :: size( i );
          const unsigned int ss = (j < ns+ms ? 0 : 1);
          if( j < ns )
            return GenericSubTopologyNumbering< BaseTopology, codim, subcodim >
                   :: number( i, j );
          else
            return GenericSubTopologyNumbering< BaseTopology, codim, subcodim-1 >
                   :: number( i, j-(ns+ss*ms) ) + nb + ss*mb;
        }
        else
          return GenericSubTopologyNumbering< BaseTopology, codim-1, subcodim >
                 :: number( i-(n+s*m), j ) + nb + s*mb;
      }
    };

    template< class BaseTopology, unsigned int codim, unsigned int subdim >
    struct GenericSubTopologyNumberingHelper
    < Prism< BaseTopology >, codim, subdim, 0 >
    {
      typedef Prism< BaseTopology > Topology;

      static unsigned int number ( unsigned int i, unsigned int j )
      {
        return i;
      }
    };

    template< class BaseTopology, unsigned int codim, unsigned int subdim >
    struct GenericSubTopologyNumberingHelper
    < Prism< BaseTopology >, codim, subdim, subdim >
    {
      typedef Prism< BaseTopology > Topology;

      enum { m = Size< BaseTopology, codim-1 > :: value };
      enum { n = Size< BaseTopology, codim > :: value };

      enum { mb = Size< BaseTopology, codim+subdim-1 > :: value };

      static unsigned int number ( unsigned int i, unsigned int j )
      {
        const unsigned int s = (i < n+m ? 0 : 1);
        if( i < n )
        {
          const unsigned int ms = SubTopologySize< BaseTopology, codim, subdim-1 > :: size( i );
          const unsigned int ss = (j < ms ? 0 : 1);
          return GenericSubTopologyNumbering< BaseTopology, codim, subdim-1 >
                 :: number( i, j-ss*ms ) + ss*mb;
        }
        else
          return GenericSubTopologyNumbering< BaseTopology, codim-1, subdim >
                 :: number( i-(n+s*m), j ) + s*mb;
      }
    };

    template< class BaseTopology, unsigned int codim,
        unsigned int subdim, unsigned int subcodim >
    struct GenericSubTopologyNumberingHelper
    < Pyramid< BaseTopology >, codim, subdim, subcodim >
    {
      typedef Pyramid< BaseTopology > Topology;

      enum { m = Size< BaseTopology, codim-1 > :: value };

      enum { mb = Size< BaseTopology, codim+subcodim-1 > :: value };

      static unsigned int number ( unsigned int i, unsigned int j )
      {
        if( i < m )
          return GenericSubTopologyNumbering< BaseTopology, codim-1, subcodim >
                 :: number( i, j );
        else
        {
          const unsigned int ms = SubTopologySize< BaseTopology, codim, subcodim-1 > :: size( i-m );
          if( j < ms )
            return GenericSubTopologyNumbering< BaseTopology, codim, subcodim-1 >
                   :: number( i-m, j );
          else
            return GenericSubTopologyNumbering< BaseTopology, codim, subcodim >
                   :: number( i-m, j-ms ) + mb;
        }
      }
    };

    template< class BaseTopology, unsigned int codim, unsigned int subdim >
    struct GenericSubTopologyNumberingHelper
    < Pyramid< BaseTopology >, codim, subdim, 0 >
    {
      typedef Pyramid< BaseTopology > Topology;

      static unsigned int number ( unsigned int i, unsigned int j )
      {
        return i;
      }
    };

    template< class BaseTopology, unsigned int codim, unsigned int subdim >
    struct GenericSubTopologyNumberingHelper
    < Pyramid< BaseTopology >, codim, subdim, subdim >
    {
      typedef Pyramid< BaseTopology > Topology;

      enum { m = Size< BaseTopology, codim-1 > :: value };

      enum { mb = Size< BaseTopology, codim+subdim-1 > :: value };

      static unsigned int number ( unsigned int i, unsigned int j )
      {
        if( i < m )
          return GenericSubTopologyNumbering< BaseTopology, codim-1, subdim >
                 :: number( i, j );
        else
        {
          const unsigned int ms = SubTopologySize< BaseTopology, codim, subdim-1 > :: size( i-m );
          if( j < ms )
            return GenericSubTopologyNumbering< BaseTopology, codim, subdim-1 >
                   :: number( i-m, j );
          else
            return mb;
        }
      }
    };

    template< class Topology, unsigned int codim, unsigned int subcodim >
    class GenericSubTopologyNumbering
    {
      dune_static_assert( (codim <= Topology :: dimension), "Invalid codimension" );
      dune_static_assert( (codim + subcodim <= Topology :: dimension),
                          "Invalid subcodimension" );

      template< bool >
      struct BorderCodim
      {
        static unsigned int number ( unsigned int i, unsigned int j )
        {
          return (codim == 0 ? j : i );
        }
      };

      template< bool >
      struct InnerCodim
      {
        static unsigned int number ( unsigned int i, unsigned int j )
        {
          return GenericSubTopologyNumberingHelper
                 < Topology, codim, Topology :: dimension - codim, subcodim >
                 :: number( i, j );
        }
      };

    public:
      static unsigned int number ( unsigned int i, unsigned int j )
      {
        assert( (j <= SubTopologySize< Topology, codim, subcodim > :: size( i )) );
        return SelectType
               < (codim == 0) || (codim == Topology :: dimension), BorderCodim<true>, InnerCodim<false> >
               :: Type :: number( i, j );
      }
    };



    // SubTopologyNumbering
    // --------------------

    template< class Topology, unsigned int codim, unsigned int subcodim >
    class SubTopologyNumbering
    {
      typedef GenericSubTopologyNumbering< Topology, codim, subcodim >
      GenericNumbering;

      std :: vector< unsigned int > numbering_[ Size< Topology, codim > :: value ];

    public:
      static unsigned int number ( unsigned int i, unsigned int j )
      {
        assert( (j <= SubTopologySize< Topology, codim, subcodim > :: size( i )) );
        return instance().numbering_[ i ][ j ];
      }

    private:
      SubTopologyNumbering ()
      {
        for( unsigned int i = 0; i < Size< Topology, codim > :: value; ++i )
        {
          const unsigned int size = SubTopologySize< Topology, codim, subcodim > :: size( i );
          numbering_[ i ].resize( size );
          for( unsigned int j = 0; j < size; ++j )
            numbering_[ i ][ j ] = GenericNumbering :: number( i, j );
        }
      }

      static const SubTopologyNumbering &instance ()
      {
        static SubTopologyNumbering inst;
        return inst;
      }
    };



    // IsCodimHybrid
    // -------------

    template< class Topology, unsigned int codim >
    struct IsCodimHybrid
    {
      static const bool value = (codim != 0) && IsHybrid< Topology >::value;
    };



    // SubTopologyMapper
    // -----------------

    template< class Topology >
    class SubTopologyMapper
    {
      static const unsigned int dimension = Topology::dimension;

      template< class A, class B >
      struct StaticSum
      {
        static const unsigned int value = A::value + B::value;
      };

      template< int codim >
      struct Size
      {
        static const unsigned int value = GenericGeometry::Size< Topology, codim >::value;
      };

      template< int codim >
      struct CalcOffset
      {
        static void apply ( unsigned int (&offsets)[ dimension+2 ] )
        {
          offsets[ codim+1 ] = offsets[ codim ] + Size< codim >::value;
        }
      };

    public:
      static const unsigned int staticSize = GenericForLoop< StaticSum, Size, 0, dimension >::value;

      SubTopologyMapper ()
      {
        offsets_[ 0 ] = 0;
        ForLoop< CalcOffset, 0, dimension >::apply( offsets_ );
        assert( size() == staticSize );
      };

      unsigned int operator() ( const unsigned int codim, const unsigned int subEntity ) const
      {
        const unsigned int offset = offsets_[ codim ];
        assert( offset + subEntity < offsets_[ codim+1 ] );
        return offset + subEntity;
      }

      unsigned int size () const
      {
        return offsets_[ dimension+1 ];
      }

    private:
      unsigned int offsets_[ dimension+2 ];
    };

  }

}

#endif // #ifndef DUNE_GENERICGEOMETRY_SUBTOPOLOGIES_HH
