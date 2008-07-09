// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_REFERENCEELEMENT_HH
#define DUNE_GENERICGEOMETRY_REFERENCEELEMENT_HH

#include <dune/common/fvector.hh>

#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< class Topology, class ctype >
    class ReferenceElement;



    // Corner
    // ------

    template< class Topology >
    class Corner;

    template<>
    class Corner< Point >
    {
      typedef Point Topology;

      template< class > friend class Corner;

      template< class ctype, int dim >
      static void evaluate_ ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        assert( i < Topology :: numCorners );
      }

    public:
      enum { numCorners = Topology :: numCorners };
      enum { dimension = Topoloty :: dimension };

      template< class ctype >
      static void evaluate ( unsigned int i, FieldVector< ctype, dimension > &x )
      {
        x = ctype( 0 );
        evaluate_( i, x );
      }
    };

    template< class BaseTopology >
    class Corner< Prism< BaseTopology > >
    {
      typedef Prism< BaseTopology > Topology;

      template< class > friend class Corner;

      template< class ctype, int dim >
      static void evaluate_ ( unsigned int i, FieldVector< ctype, dim > &x )
      {
        assert( i < Topology :: numCorners );
        const unsigned int j = i % BaseTopology :: numCorners;
        Corner< BaseTopology > :: evaluate_( j, x );
        if( j >= BaseTopology :: numCorners )
          x[ dimension - 1 ] = ctype( 1 );
      }

    public:
      enum { numCorners = Topology :: numCorners };
      enum { dimension = Topology :: dimension };

      template< class ctype >
      static void evaluate ( unsigned int i, FieldVector< ctype, dimension > &x )
      {
        x = ctype( 0 );
        evaluate_( i, x );
      }
    };

    template< class BaseTopology >
    class Corner< Pyramid< BaseTopology > >
    {
      typedef Prism< BaseTopology > Topology;

      template< class > friend class Corner;

      template< class ctype, int dim >
      static void evaluate_ ( unsigned int i, FieldVector< ctype, dim > &x )
      {
        assert( i < Topology :: numCorners );
        if( i < BaseTopology :: numCorners )
          Corner< BaseTopology > :: evaluate_( i, x );
        else
          x[ dimension - 1 ] = ctype( 1 );
      }

    public:
      enum { numCorners = Topology :: numCorners };
      enum { dimension = Topology :: dimension };

      template< class ctype >
      static void evaluate ( unsigned int i, FieldVector< ctype, dimension > &x )
      {
        x = ctype( 0 );
        evaluate_( i, x );
      }
    };



    // IntegrationOuterNormal
    // ----------------------

    template< class Topology >
    class IntegrationOuterNormal;

    template< class BaseTopology >
    class IntegrationOuterNormal< Prism< BaseTopology > >
    {
      typedef Prism< BaseTopology > Topology;

      enum { dimension = Topology :: dimension };

      template< class > friend class IntegrationOuterNormal;

      template< class ctype, int dim >
      static void evaluate_ ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        if( i >= Size< BaseTopology, 1 > :: value )
        {
          const unsigned int j = i - Size< BaseTopology, 1 > :: value;
          n[ dimension - 1 ] = (j == 0 ? ctype( -1 ) : ctype( 1 ));
        }
        else
          IntegrationOuterNormal< BaseTopology > :: evaluate_( i, n );
      }

    public:
      template< class ctype >
      static void evaluate ( unsigned int i, FieldVector< ctype, dimension > &n )
      {
        n = ctype( 0 );
        evaluate_( i, n );
      }
    };

    template<>
    class IntegrationOuterNormal< Prism< Point > >
    {
      typedef Prism< Point > Topology;

      enum { dimension = Topology :: dimension };

      template< class > friend class IntegrationOuterNormal;

      template< class ctype, int dim >
      static void evaluate_ ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        n[ dimension - 1 ] = (i == 0 ? ctype( -1 ) : ctype( 1 ));
      }

    public:
      template< class ctype >
      static void evaluate ( unsigned int i, FieldVector< ctype, dimension > &n )
      {
        n = ctype( 0 );
        evaluate_( i, n );
      }
    };

    template< class BaseTopology >
    class IntegrationOuterNormal< Pyramid< BaseTopology > >
    {
      typedef Pyramid< BaseTopology > Topology;

      enum { dimension = Topology :: dimension };

      template< class > friend class IntegrationOuterNormal;

      template< class ctype, int dim >
      static void evaluate_ ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        typedef SubTopologyNumbering< BaseTopology > Numbering;

        if( i < Size< BaseTopology, 1 > :: value )
        {
          const unsigned int j
            = Numbering :: template subEntity< 1, dimension-2 >( i, 0 );
          FieldVector< ctype, dim > x( ctype( 0 ) );
          Corner< BaseTopology > :: evaluate_( j, x );

          IntegrationOuterNormal< BaseTopology > :: evaluate_( i, n );
          n[ dimension - 1 ] = (x * n);
        }
        else
          n[ dimension - 1 ] = ctype( -1 );
      }

    public:
      template< class ctype >
      static void evaluate ( unsigned int i, FieldVector< ctype, dimension > &n )
      {
        n = ctype( 0 );
        evaluate_( i, n );
      }
    };

    template<>
    class IntegrationOuterNormal< Pyramid< Point > >
    {
      typedef Pyramid< Point > Topology;

      enum { dimension = Topology :: dimension };

      template< class > friend class IntegrationOuterNormal;

      template< class ctype, int dim >
      static void evaluate_ ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        n[ dimension - 1 ] = (i == 0 ? ctype( -1 ) : ctype( 1 ));
      }

    public:
      template< class ctype >
      static void evaluate ( unsigned int i, FieldVector< ctype, dimension > &n )
      {
        n = ctype( 0 );
        evaluate_( i, n );
      }
    };



    // Volume
    // ------

    template< class Topology >
    struct Volume;

    template<>
    struct Volume< Point >
    {
      template< class ctype >
      static ctype evaluate ()
      {
        return ctype( 1 );
      }
    };

    template< class BaseTopology >
    struct Volume< Prism< BaseTopology > >
    {
      template< class ctype >
      static ctype evaluate ()
      {
        return Volume< BaseTopology > :: template evaluate< ctype >();
      }
    };

    template< class BaseTopology >
    struct Volume< Pyramid< BaseTopology > >
    {
      template< class ctype >
      static ctype evaluate ()
      {
        const ctype baseVolume
          = Volume< BaseTopology > :: template evaluate< ctype >();

        return baseVolume / ctype( Faculty< dimension > :: value );
      }
    };



    // ReferenceElement
    // ----------------

    template< class Topology, class ctype >
    struct ReferenceElement
    {
      enum { numCorners = Topology :: numCorners };
      enum { dimension = Topology :: dimension };

      typedef FieldVector< ctype, dimension > CoordinateType;

      template< unsigned int codim, unsigned int subcodim >
      static unsigned int subNumbering ( unsigned int i, unsigned int j )
      {
        return SubTopologyNumbering< Topology, codim, subcodim > :: number( i, j );
      }

      static const CoordinateType &corner ( unsigned int i )
      {
        return instance().corners_[ i ];
      }

      template< class ctype >
      static const CoordinateType &
      integrationOuterNormal ( unsigned int i )
      {
        return instance().normals_[ i ];
      }

      static ctype volume ()
      {
        return Volume< Topology > :: template evaluate< ctype >();
      }

    private:
      enum { numFaces = Size< Topology, 1 > :: value };

      ReferenceElement ()
      {
        for( unsigned int i = 0; i < numCorners; ++i )
          Corners< Topology > :: evaluate( i, corners_[ i ] );
        for( unsigned int i = 0; i < numFaces; ++i )
          IntegrationOuterNormal< Topology > :: evaluate( i, normals_[ i ] );
      }

      static const ReferenceElement &instance ()
      {
        static ReferenceElement inst;
        return inst;
      }

      CoordinateType corners_[ numCorners ];
      CoordinateType normals_[ numFaces ];
    };

    template< class ctype >
    struct ReferenceElement< Point, ctype >
    {
      typedef Point Topology;

      enum { numCorners = Topology :: numCorners };
      enum { dimension = Topology :: dimension };

      typedef FieldVector< ctype, dimension > CoordinateType;

      template< unsigned int codim, unsigned int subcodim >
      static unsigned int subNumbering ( unsigned int i, unsigned int j )
      {
        return SubTopologyNumbering< Topology, codim, subcodim > :: number( i, j );
      }

      static const CoordinateType &corner ( unsigned int i )
      {
        return instance().corners_[ i ];
      }

      template< class ctype >
      static const CoordinateType &
      integrationOuterNormal ( unsigned int i )
      {
        abort();
      }

      static ctype volume ()
      {
        return Volume< Topology > :: template evaluate< ctype >();
      }

    private:
      ReferenceElement ()
      {
        for( unsigned int i = 0; i < numCorners; ++i )
          Corners< Topology > :: evaluate( i, corners_[ i ] );
      }

      static const ReferenceElement &instance ()
      {
        static ReferenceElement inst;
        return inst;
      }

      CoordinateType corners_[ numCorners ];
    };

  }

}

#endif
