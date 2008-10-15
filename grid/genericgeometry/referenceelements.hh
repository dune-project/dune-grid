// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_REFERENCEELEMENTS_HH
#define DUNE_GENERICGEOMETRY_REFERENCEELEMENTS_HH

#include <dune/common/fixedarray.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/genericgeometry/conversion.hh>
#include <dune/grid/genericgeometry/subtopologies.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< class Topology, class ctype >
    class ReferenceElement;



    // ReferenceDomain
    // ---------------

    template< class Topology >
    struct ReferenceDomainBase;

    template<>
    struct ReferenceDomainBase< Point >
    {
      typedef Point Topology;

      template< class > friend class ReferenceDomain;
      template< class > friend class ReferenceDomainBase;

      enum { numNormals = 0 };

      template< class ctype, int dim >
      static void corner ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        assert( i < Topology :: numCorners );
      }

      template< class ctype, int dim >
      static bool
      checkInside ( const FieldVector< ctype, dim > &x, ctype factor )
      {
        return true;
      }

      template< class ctype, int dim >
      static void
      integrationOuterNormal ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        assert( i < numNormals );
      }

      template< class ctype >
      static ctype volume ()
      {
        return ctype( 1 );
      }
    };

    template< class BaseTopology >
    struct ReferenceDomainBase< Prism< BaseTopology > >
    {
      typedef Prism< BaseTopology > Topology;

      template< class > friend class ReferenceDomain;
      template< class > friend class ReferenceDomainBase;

      enum { numNormals = Size< Topology, 1 > :: value };

      static const unsigned int dimension = Topology :: dimension;
      static const unsigned int myindex = dimension - 1;

      template< bool >
      struct MultiDimensional
      {
        template< class ctype, int dim >
        static void
        integrationOuterNormal ( unsigned int i, FieldVector< ctype, dim > &n )
        {
          if( i >= Size< BaseTopology, 1 > :: value )
          {
            const unsigned int j = i - Size< BaseTopology, 1 > :: value;
            n[ myindex ] = (j == 0 ? ctype( -1 ) : ctype( 1 ));
          }
          else
            ReferenceDomainBase< BaseTopology > :: integrationOuterNormal( i, n );
        }
      };

      template< bool >
      struct OneDimensional
      {
        template< class ctype, int dim >
        static void
        integrationOuterNormal ( unsigned int i, FieldVector< ctype, dim > &n )
        {
          n[ myindex ] = (i > 0) ? ctype( 1 ) : ctype( -1 );
        }
      };


      template< class ctype, int dim >
      static void corner ( unsigned int i, FieldVector< ctype, dim > &x )
      {
        assert( i < Topology :: numCorners );
        const unsigned int j = i % BaseTopology :: numCorners;
        ReferenceDomainBase< BaseTopology > :: corner( j, x );
        if( i >= BaseTopology :: numCorners )
          x[ myindex ] = ctype( 1 );
      }

      template< class ctype, int dim >
      static bool
      checkInside ( const FieldVector< ctype, dim > &x, ctype factor )
      {
        const ctype xn = x[ myindex ];
        const ctype cxn = factor - xn;
        return (xn > -1e-12) && (cxn > -1e-12)
               && ReferenceDomainBase< BaseTopology > :: checkInside( x, factor );
      }

      template< class ctype, int dim >
      static void
      integrationOuterNormal ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        ProtectedIf< (dimension > 1), MultiDimensional, OneDimensional >
        :: integrationOuterNormal( i, n );
      }

      template< class ctype >
      static ctype volume ()
      {
        return ReferenceDomainBase< BaseTopology > :: template volume< ctype >();
      }
    };

    template< class BaseTopology >
    struct ReferenceDomainBase< Pyramid< BaseTopology > >
    {
      typedef Pyramid< BaseTopology > Topology;

      template< class > friend class ReferenceDomain;
      template< class > friend class ReferenceDomainBase;

      enum { numNormals = Size< Topology, 1 > :: value };

      static const unsigned int dimension = Topology :: dimension;
      static const unsigned int myindex = dimension - 1;

      template< bool >
      struct MultiDimensional
      {
        template< class ctype, int dim >
        static void
        integrationOuterNormal ( unsigned int i, FieldVector< ctype, dim > &n )
        {
          typedef SubTopologyNumbering< BaseTopology,1,dimension-2 > Numbering;
          if( i > 0 )
          {
            const unsigned int j = Numbering :: number( i-1, 0 );
            FieldVector< ctype, dim > x( ctype( 0 ) );
            ReferenceDomainBase< BaseTopology > :: corner( j, x );

            ReferenceDomainBase< BaseTopology > :: integrationOuterNormal ( i-1, n );
            n[ myindex ] = (x * n);
          }
          else
            n[ myindex ] = ctype( -1 );
        }
      };

      template< bool >
      struct OneDimensional
      {
        template< class ctype, int dim >
        static void
        integrationOuterNormal ( unsigned int i, FieldVector< ctype, dim > &n )
        {
          n[ myindex ] = (i > 0) ? ctype( 1 ) : ctype( -1 );
        }
      };

      template< class ctype, int dim >
      static void corner ( unsigned int i, FieldVector< ctype, dim > &x )
      {
        assert( i < Topology :: numCorners );
        if( i < BaseTopology :: numCorners )
          ReferenceDomainBase< BaseTopology > :: corner( i, x );
        else
          x[ myindex ] = ctype( 1 );
      }

      template< class ctype, int dim >
      static bool
      checkInside ( const FieldVector< ctype, dim > &x, ctype factor )
      {
        const ctype xn = x[ myindex ];
        const ctype cxn = factor - xn;
        return (xn > -1e-12) && (cxn > -1e-12)
               && ReferenceDomainBase< BaseTopology > :: checkInside( x, factor * cxn );
      }

      template< class ctype, int dim >
      static void
      integrationOuterNormal ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        ProtectedIf< (dimension > 1), MultiDimensional, OneDimensional >
        :: integrationOuterNormal( i, n );
      }

      template< class ctype >
      static ctype volume ()
      {
        const ctype baseVolume
          = ReferenceDomainBase< BaseTopology > :: template volume< ctype >();
        return baseVolume / ctype( dimension );
      }
    };



    template< class Topology >
    struct ReferenceDomain
    {
      static const unsigned int numCorners = Topology :: numCorners;
      static const unsigned int dimension = Topology :: dimension;
      static const unsigned int numNormals
        = ReferenceDomainBase< Topology > :: numNormals;

      template< class ctype >
      static void corner ( unsigned int i, FieldVector< ctype, dimension > &x )
      {
        x = ctype( 0 );
        ReferenceDomainBase< Topology > :: corner( i, x );
      }

      template< class ctype >
      static bool checkInside ( const FieldVector< ctype, dimension > &x )
      {
        return ReferenceDomainBase< Topology > :: checkInside( x, ctype( 1 ) );
      }

      template< class ctype >
      static void
      integrationOuterNormal ( unsigned int i, FieldVector< ctype, dimension > &n )
      {
        n = ctype( 0 );
        return ReferenceDomainBase< Topology > :: integrationOuterNormal( i, n );
      }

      template< class ctype >
      static ctype volume ()
      {
        return ReferenceDomainBase< Topology > :: template volume< ctype >();
      }
    };



    // ReferenceElement
    // ----------------

    template< class Topology, class ctype >
    struct ReferenceElement
    {
      static const unsigned int topologyId = Topology :: id;
      static const unsigned int dimension = Topology :: dimension;

      static const unsigned int numCorners = Topology :: numCorners;
      static const unsigned int numNormals = ReferenceDomain< Topology > :: numNormals;

      typedef FieldVector< ctype, dimension > CoordinateType;

      template< unsigned int codim >
      struct Codim
      {
        enum { size = Size< Topology, codim > :: value };
      };

      template< unsigned int codim, unsigned int subcodim >
      static unsigned int subNumbering ( unsigned int i, unsigned int j )
      {
        return SubTopologyNumbering< Topology, codim, subcodim > :: number( i, j );
      }

      template< unsigned int codim, unsigned int subcodim >
      static unsigned int size ( unsigned int i )
      {
        return SubTopologySize< Topology, codim, subcodim > :: size( i );
      }

      template< unsigned int codim >
      static const FieldVector< ctype, dimension > &
      baryCenter ( unsigned int i )
      {
        Int2Type< codim > codimVariable;
        return instance().baryCenters_[ codimVariable ][ i ];
      }

      static const CoordinateType &corner ( unsigned int i )
      {
        assert( i < numCorners );
        return instance().corners_[ i ];
      }

      static bool checkInside ( const CoordinateType &x )
      {
        return ReferenceDomain< Topology > :: checkInside( x );
      }

      static const CoordinateType &
      integrationOuterNormal ( unsigned int i )
      {
        assert( i < numNormals );
        return instance().normals_[ i ];
      }

      static ctype volume ()
      {
        return ReferenceDomain< Topology > :: template volume< ctype >();
      }

      static const ReferenceElement &instance ()
      {
        static ReferenceElement inst;
        return inst;
      }

    private:
      template< int codim >
      class BaryCenterArray;

      ReferenceElement ()
      {
        for( unsigned int i = 0; i < numCorners; ++i )
          ReferenceDomain< Topology > :: corner( i, corners_[ i ] );
        for( unsigned int i = 0; i < numNormals; ++i )
          ReferenceDomain< Topology > :: integrationOuterNormal( i, normals_[ i ] );
      }

      Dune::array< CoordinateType, numCorners > corners_;
      CodimTable< BaryCenterArray, dimension > baryCenters_;
      Dune::array< CoordinateType, numNormals > normals_;
    };



    template< class Topology, class ctype >
    template< int codim >
    class ReferenceElement< Topology, ctype > :: BaryCenterArray
    {
      enum { Size = GenericGeometry :: Size< Topology, codim > :: value };

      typedef FieldVector< ctype, dimension > CoordinateType;

      template< int i >
      struct Builder;

      CoordinateType baryCenters_[ Size ];

    public:
      BaryCenterArray ()
      {
        ForLoop< Builder, 0, Size-1 > :: apply( baryCenters_ );
      }

      const CoordinateType &operator[] ( unsigned int i ) const
      {
        assert( i < Size );
        return baryCenters_[ i ];
      }

      static unsigned int size ()
      {
        return Size;
      }
    };

    template< class Topology, class ctype >
    template< int codim >
    template< int i >
    struct ReferenceElement< Topology, ctype > :: BaryCenterArray< codim > :: Builder
    {
      static void apply ( CoordinateType (&baryCenters)[ Size ] )
      {
        typedef SubTopologyNumbering< Topology, codim, dimension - codim > Numbering;
        typedef SubTopologySize< Topology, codim, dimension - codim > Size;

        CoordinateType &x = baryCenters[ i ];
        x = 0;
        const unsigned int numCorners = Size :: size( i );
        for( unsigned int k = 0; k < numCorners; ++k )
        {
          unsigned int j = Numbering :: number( i, k );

          CoordinateType y;
          ReferenceDomain< Topology > :: corner( j, y );
          x += y;
        }
        x *= ctype( 1 ) / ctype( numCorners );
      }
    };

  }

}

#endif
