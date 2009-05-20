// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_REFERENCEDOMAIN_HH
#define DUNE_GENERICGEOMETRY_REFERENCEDOMAIN_HH

#include <dune/common/fixedarray.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/genericgeometry/subtopologies.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class Topology >
    struct ReferenceDomain;



    // ReferenceDomain
    // ---------------

    template< class Topology >
    struct ReferenceDomainBase;

    /** \cond */
    template<>
    class ReferenceDomainBase< Point >
    {
      typedef Point Topology;

      friend class ReferenceDomain< Topology >;
      friend class ReferenceDomainBase< Prism< Topology > >;
      friend class ReferenceDomainBase< Pyramid< Topology > >;

      static const unsigned int numNormals = 0;

      template< class ctype, int dim >
      static void corner ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        assert( i < Topology::numCorners );
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
    class ReferenceDomainBase< Prism< BaseTopology > >
    {
      typedef Prism< BaseTopology > Topology;

      friend class ReferenceDomain< Topology >;
      friend class ReferenceDomainBase< Prism< Topology > >;
      friend class ReferenceDomainBase< Pyramid< Topology > >;

      static const unsigned int numNormals = Size< Topology, 1 >::value;

      static const unsigned int dimension = Topology::dimension;
      static const unsigned int myindex = dimension - 1;

      template< class ctype, int dim >
      static void corner ( unsigned int i, FieldVector< ctype, dim > &x )
      {
        assert( i < Topology::numCorners );
        const unsigned int j = i % BaseTopology::numCorners;
        ReferenceDomainBase< BaseTopology >::corner( j, x );
        if( i >= BaseTopology::numCorners )
          x[ myindex ] = ctype( 1 );
      }

      template< class ctype, int dim >
      static bool
      checkInside ( const FieldVector< ctype, dim > &x, ctype factor )
      {
        const ctype xn = x[ myindex ];
        const ctype cxn = factor - xn;
        return (xn > -1e-12) && (cxn > -1e-12)
               && ReferenceDomainBase< BaseTopology >::checkInside( x, factor );
      }

      template< class ctype, int dim >
      static void
      integrationOuterNormal ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        typedef ReferenceDomainBase< BaseTopology > BaseReferenceDomain;

        if( i >= BaseReferenceDomain::numNormals )
        {
          const unsigned int j = i - BaseReferenceDomain::numNormals;
          n[ myindex ] = (j == 0 ? ctype( -1 ) : ctype( 1 ));
        }
        else
          BaseReferenceDomain::integrationOuterNormal( i, n );
      }

      template< class ctype >
      static ctype volume ()
      {
        typedef ReferenceDomainBase< BaseTopology > BaseReferenceDomain;
        return BaseReferenceDomain::template volume< ctype >();
      }
    };


    template< class BaseTopology >
    class ReferenceDomainBase< Pyramid< BaseTopology > >
    {
      typedef Pyramid< BaseTopology > Topology;

      friend class ReferenceDomain< Topology >;
      friend class ReferenceDomainBase< Prism< Topology > >;
      friend class ReferenceDomainBase< Pyramid< Topology > >;

      static const unsigned int numNormals = Size< Topology, 1 >::value;

      static const unsigned int dimension = Topology::dimension;
      static const unsigned int myindex = dimension - 1;

      template< bool >
      struct MultiDimensional
      {
        template< class ctype, int dim >
        static void
        integrationOuterNormal ( unsigned int i, FieldVector< ctype, dim > &n )
        {
          multiDimensionalIntegrationOuterNormal( i, n );
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
        assert( i < Topology::numCorners );
        if( i < BaseTopology::numCorners )
          ReferenceDomainBase< BaseTopology >::corner( i, x );
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
               && ReferenceDomainBase< BaseTopology >::checkInside( x, cxn );
      }

      template< class ctype, int dim >
      static void
      multiDimensionalIntegrationOuterNormal ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        typedef ReferenceDomainBase< BaseTopology > BaseReferenceDomain;
        typedef SubTopologyNumbering< BaseTopology, 1, dimension-2 > Numbering;

        if( i > 0 )
        {
          const unsigned int j = Numbering::number( i-1, 0 );
          FieldVector< ctype, dim > x( ctype( 0 ) );
          BaseReferenceDomain::corner( j, x );

          BaseReferenceDomain::integrationOuterNormal ( i-1, n );
          n[ myindex ] = (x * n);
        }
        else
          n[ myindex ] = ctype( -1 );
      }

      template< class ctype, int dim >
      static void
      integrationOuterNormal ( unsigned int i, FieldVector< ctype, dim > &n )
      {
        ProtectedIf< (dimension > 1), MultiDimensional, OneDimensional >
        ::integrationOuterNormal( i, n );
      }

      template< class ctype >
      static ctype volume ()
      {
        typedef ReferenceDomainBase< BaseTopology > BaseReferenceDomain;
        const ctype baseVolume = BaseReferenceDomain::template volume< ctype >();
        return baseVolume / ctype( (unsigned int)(dimension) ); // linker problem when using dimension directly
      }
    };
    /** \endcond */



    // ReferenceDomain
    // ---------------

    template< class Topology >
    struct ReferenceDomain
    {
      static const unsigned int numCorners = Topology::numCorners;
      static const unsigned int dimension = Topology::dimension;

      static const unsigned int numNormals
        = ReferenceDomainBase< Topology >::numNormals;

      template< class ctype >
      static void corner ( unsigned int i, FieldVector< ctype, dimension > &x )
      {
        x = ctype( 0 );
        ReferenceDomainBase< Topology >::corner( i, x );
      }

      template< class ctype >
      static bool checkInside ( const FieldVector< ctype, dimension > &x )
      {
        return ReferenceDomainBase< Topology >::checkInside( x, ctype( 1 ) );
      }

      template< class ctype >
      static void
      integrationOuterNormal ( unsigned int i, FieldVector< ctype, dimension > &n )
      {
        n = ctype( 0 );
        return ReferenceDomainBase< Topology >::integrationOuterNormal( i, n );
      }

      template< class ctype >
      static ctype volume ()
      {
        return ReferenceDomainBase< Topology >::template volume< ctype >();
      }
    };

  }

}

#endif
