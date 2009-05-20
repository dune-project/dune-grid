// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_REFERENCEELEMENTS_HH
#define DUNE_GENERICGEOMETRY_REFERENCEELEMENTS_HH

#include <dune/common/fixedarray.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/genericgeometry/referencedomain.hh>

namespace Dune
{

  namespace GenericGeometry
  {

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
        return ReferenceDomain< Topology >::checkInside( x );
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
