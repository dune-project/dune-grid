// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_CORNERMAPPING_HH
#define DUNE_GENERICGEOMETRY_CORNERMAPPING_HH

#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>
#include <dune/grid/genericgeometry/referenceelements.hh>
#include <dune/grid/genericgeometry/matrix.hh>
#include <dune/grid/genericgeometry/submapping.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    /* Mapping have two template arguments:
     * Topology:    the generic geometry describing
     *              the domain.
     * CoordTraits: a traits class describing the
     *              vector and derivative types needed
     *              to describe the range of the mapping.
     *              This class fixes the local (domain) vector type
     *              and the range (world) vector type. Note
     *              that the dimension of the local coordinates (dimG)
     *              must be greater or equal to Topology::dimension;
     *              if dimG > Topology::dimension then only the
     *              first components are used.
     *              The dimension of the global coordinates (dimW)
     *              must be greater or equal to dimG.
     *              Matrix types both for the jacobian (dimW x dimG)
     *              its transpose (dimG x dimW) and
     *              a square matrix of dimension (dimG x dimG).
     * struct CoordTraits {
     *   enum {dimW = };              // world dimension
     *   enum {dimG = };              // grid dimension
     *   typedef ... FieldType;
     *   // general vector and matrix types
     *   template <int dim> struct Vector { typedef ... Type; };
     *   template <int dimR,dimC> struct Matrix { typedef ... Type; };
     *
     *   // Vector of global vectors denoting the edges of the range
     *   // domain, used to construct a mapping together with an offset.
     *   // Corners used are
     *   // p[offset],...,p[offset+Topology::numCorners]
     *   // Assumption: coord_vector[i] is of type const Vector<dimW>&
     *   typedef ... coord_vector;
     *
     *   // mapping is of the form Ax+b (used untested)
     *   enum {affine = 0|1};
     * };
     */



    // GenericMapping
    // --------------

    template< class Topology, class Traits, unsigned int offset = 0 >
    class GenericMapping;

    template< class Traits, unsigned int offset >
    class GenericMapping < Point, Traits, offset >
    {
      typedef Point Topology;

    public:
      enum { dim = Topology :: dimension };

      enum { dimW = Traits :: dimW };
      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianTransposedType JacobianTransposedType;

      enum { alwaysAffine = true };

      static bool isZero ( const FieldType &a )
      {
        return std::abs(a)<1e-12;
      }

      template< unsigned int numCorners >
      static const GlobalCoordType &
      origin ( const GlobalCoordType *const (&coords)[ numCorners ] )
      {
        dune_static_assert( (offset < numCorners), "Invalid offset." );
        return *(coords[ offset ]);
      }

      template< unsigned int numCorners >
      static void phi_set ( const GlobalCoordType *const (&coords)[ numCorners ],
                            const LocalCoordType &x,
                            const FieldType &factor,
                            GlobalCoordType &p )
      {
        const GlobalCoordType &y = origin( coords );
        for( unsigned int i = 0; i < dimW; ++i )
          p[ i ] = factor * y[ i ];
      }

      template< unsigned int numCorners >
      static void phi_add ( const GlobalCoordType *const (&coords)[ numCorners ],
                            const LocalCoordType &x,
                            const FieldType &factor,
                            GlobalCoordType &p )
      {
        const GlobalCoordType &y = origin( coords );
        for( unsigned int i = 0; i < dimW; ++i )
          p[ i ] += factor * y[ i ];
      }

      template< unsigned int numCorners >
      static bool Dphi_set ( const GlobalCoordType *const (&coords)[ numCorners ],
                             const LocalCoordType &x,
                             const FieldType &factor,
                             JacobianTransposedType &J )
      {
        return true;
      }

      template< unsigned int numCorners >
      static bool Dphi_add ( const GlobalCoordType *const (&coords)[ numCorners ],
                             const LocalCoordType &x,
                             const FieldType &factor,
                             JacobianTransposedType &J )
      {
        return true;
      }


      static bool inDomain ( const LocalCoordType &x, FieldType factor )
      {
        return isZero( x[ 0 ] );
      }

    };


    template< class BaseTopology, class Traits, unsigned int offset >
    class GenericMapping< Prism< BaseTopology >, Traits, offset >
    {
      typedef Prism< BaseTopology > Topology;

      typedef GenericMapping< BaseTopology, Traits, offset > BottomMapping;
      typedef GenericMapping< BaseTopology, Traits, offset + BaseTopology :: numCorners >
      TopMapping;
    public:
      enum { dim  = Topology :: dimension };

      enum { dimW = Traits :: dimW };
      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianTransposedType JacobianTransposedType;

      enum { alwaysAffine = ((dim < 2) || Traits :: affine) };

      template< unsigned int numCorners >
      static const GlobalCoordType &
      origin ( const GlobalCoordType *const (&coords)[ numCorners ] )
      {
        return BottomMapping :: origin( coords );
      }

      template< unsigned int numCorners >
      static void phi_set ( const GlobalCoordType *const (&coords)[ numCorners ],
                            const LocalCoordType &x,
                            const FieldType &factor,
                            GlobalCoordType &p )
      {
        const FieldType xn = x[ dim-1 ];
        const FieldType cxn = FieldType( 1 ) - xn;
        BottomMapping :: phi_set( coords, x, factor * cxn, p );
        TopMapping :: phi_add( coords, x, factor * xn, p );
      }

      template< unsigned int numCorners >
      static void phi_add ( const GlobalCoordType *const (&coords)[ numCorners ],
                            const LocalCoordType &x,
                            const FieldType &factor,
                            GlobalCoordType &p )
      {
        const FieldType xn = x[ dim-1 ];
        const FieldType cxn = FieldType( 1 ) - xn;
        BottomMapping :: phi_add( coords, x, factor * cxn, p );
        TopMapping :: phi_add( coords, x, factor * xn, p );
      }

      template< unsigned int numCorners >
      static bool Dphi_set ( const GlobalCoordType *const (&coords)[ numCorners ],
                             const LocalCoordType &x,
                             const FieldType &factor,
                             JacobianTransposedType &J )
      {
        const FieldType xn = x[ dim-1 ];
        bool affine = true;
        if( alwaysAffine )
        {
          const FieldType cxn = FieldType( 1 ) - xn;
          BottomMapping :: Dphi_set( coords, x, factor * cxn, J );
          TopMapping :: Dphi_add( coords, x, factor * xn, J );
        }
        else
        {
          JacobianTransposedType Jtop;
          affine &= BottomMapping :: Dphi_set( coords, x, factor, J );
          affine &= TopMapping :: Dphi_set( coords, x, factor, Jtop );

          FieldType norm = FieldType( 0 );
          for( unsigned int i = 0; i < dim-1; ++i )
          {
            Jtop[ i ] -= J[ i ];
            norm += Jtop[ i ].two_norm2();
            J[ i ].axpy( xn, Jtop[ i ] );
          }
          affine &= (norm < 1e-12);
        }
        BottomMapping :: phi_set( coords, x, -factor, J[ dim-1 ] );
        TopMapping :: phi_add( coords, x, factor, J[ dim-1 ] );
        return affine;
      }

      template< unsigned int numCorners >
      static bool Dphi_add ( const GlobalCoordType *const (&coords)[ numCorners ],
                             const LocalCoordType &x,
                             const FieldType &factor,
                             JacobianTransposedType &J )
      {
        const FieldType xn = x[ dim-1 ];
        bool affine = true;
        if( alwaysAffine )
        {
          const FieldType cxn = FieldType( 1 ) - xn;
          BottomMapping :: Dphi_add( coords, x, factor * cxn, J );
          TopMapping :: Dphi_add( coords, x, factor * xn, J );
        }
        else
        {
          JacobianTransposedType Jbottom, Jtop;
          affine &= BottomMapping :: Dphi_set( coords, x, FieldType( 1 ), Jbottom );
          affine &= TopMapping :: Dphi_set( coords, x, FieldType( 1 ), Jtop );

          FieldType norm = FieldType( 0 );
          for( unsigned int i = 0; i < dim-1; ++i )
          {
            Jtop[ i ] -= Jbottom[ i ];
            norm += Jtop[ i ].two_norm2();
            J[ i ].axpy( factor, Jbottom[ i ] );
            J[ i ].axpy( factor*xn, Jtop[ i ] );
          }
          affine &= (norm < 1e-12);
        }
        BottomMapping :: phi_add( coords, x, -factor, J[ dim-1 ] );
        TopMapping :: phi_add( coords, x, factor, J[ dim-1 ] );
        return affine;
      }

      // check if x/fac is in domain of phi
      static bool inDomain ( const LocalCoordType &x, FieldType factor )
      {
        const FieldType xn = x[ dim-1 ];
        const FieldType cxn = factor - xn;
        return (xn > -1e-12) && (cxn > -1e-12)
               && BottomMapping :: inDomain( x, factor );
      }

    };

    template< class BaseTopology, class Traits, unsigned int offset >
    class GenericMapping < Pyramid< BaseTopology >, Traits, offset >
    {
      typedef Pyramid< BaseTopology > Topology;

      typedef GenericMapping< BaseTopology, Traits, offset > BottomMapping;
      typedef GenericMapping< Point, Traits, offset + BaseTopology :: numCorners >
      TopMapping;

    public:
      enum {dim  = Topology::dimension};

      enum { dimW = Traits :: dimW };
      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianTransposedType JacobianTransposedType;

      enum { alwaysAffine = (BottomMapping :: alwaysAffine || Traits :: affine) };

      template< unsigned int numCorners >
      static const GlobalCoordType &
      origin ( const GlobalCoordType *const (&coords)[ numCorners ] )
      {
        return BottomMapping :: origin( coords );
      }

      template< unsigned int numCorners >
      static void phi_set ( const GlobalCoordType *const (&coords)[ numCorners ],
                            const LocalCoordType &x,
                            const FieldType &factor,
                            GlobalCoordType &p )
      {
        const FieldType xn = x[ dim-1 ];
        if( alwaysAffine )
        {
          const GlobalCoordType &top = TopMapping :: origin( coords );
          const GlobalCoordType &bottom = BottomMapping :: origin( coords );

          BottomMapping :: phi_set( coords, x, factor, p );
          for( unsigned int i = 0; i < dimW; ++i )
            p[ i ] += (factor * xn) * (top[ i ] - bottom[ i ]);
        }
        else
        {
          TopMapping :: phi_set( coords, x, factor * xn, p );
          const FieldType cxn = FieldType( 1 ) - xn;
          if( cxn > 1e-12 )
          {
            const FieldType icxn = FieldType( 1 ) / cxn;
            LocalCoordType xb;
            for( unsigned int i = 0; i < dim-1; ++i )
              xb[ i ] = icxn * x[ i ];

            BottomMapping :: phi_add( coords, xb, factor * cxn, p );
          }
        }
      }

      template< unsigned int numCorners >
      static void phi_add ( const GlobalCoordType *const (&coords)[ numCorners ],
                            const LocalCoordType &x,
                            const FieldType &factor,
                            GlobalCoordType &p )
      {
        const FieldType xn = x[ dim-1 ];
        if( alwaysAffine )
        {
          const GlobalCoordType &top = TopMapping :: origin( coords );
          const GlobalCoordType &bottom = BottomMapping :: origin( coords );

          BottomMapping :: phi_add( coords, x, factor, p );
          for( unsigned int i = 0; i < dimW; ++i )
            p[ i ] += (factor * xn) * (top[ i ] - bottom[ i ]);
        }
        else
        {
          TopMapping :: phi_add( coords, x, factor * xn, p );
          const FieldType cxn = FieldType( 1 ) - xn;
          if( cxn > 1e-12 )
          {
            const FieldType icxn = FieldType( 1 ) / cxn;
            LocalCoordType xb;
            for( unsigned int i = 0; i < dim-1; ++i )
              xb[ i ] = icxn * x[ i ];

            BottomMapping :: phi_add( coords, xb, factor * cxn, p );
          }
        }
      }

      template< unsigned int numCorners >
      static bool Dphi_set ( const GlobalCoordType *const (&coords)[ numCorners ],
                             const LocalCoordType &x,
                             const FieldType &factor,
                             JacobianTransposedType &J )
      {
        GlobalCoordType &q = J[ dim-1 ];
        bool affine;
        if( alwaysAffine )
        {
          const GlobalCoordType &top = TopMapping :: origin( coords );
          const GlobalCoordType &bottom = BottomMapping :: origin( coords );

          affine = BottomMapping :: Dphi_set( coords, x, factor, J );
          for( unsigned int i = 0; i < dimW; ++i )
            q[ i ] = factor * (top[ i ] - bottom[ i ]);
        }
        else
        {
          const FieldType xn = x[ dim-1 ];
          const FieldType cxn = FieldType( 1 ) - xn;
          const FieldType icxn = FieldType( 1 ) / cxn;
          LocalCoordType xb;
          for( unsigned int i = 0; i < dim-1; ++i )
            xb[ i ] = icxn * x[ i ];
          affine = BottomMapping :: Dphi_set( coords, xb, factor, J );

          TopMapping :: phi_set( coords, x, factor, q );
          BottomMapping :: phi_add( coords, xb, -factor, q );
          xb *= factor;
          for( unsigned int j = 0; j < dim-1; ++j )
          {
            for( unsigned int i = 0; i < dimW; ++i )
              q[ i ] += J[ j ][ i ] * xb[ j ];
          }
        }
        return affine;
      }

      template< unsigned int numCorners >
      static bool Dphi_add ( const GlobalCoordType *const (&coords)[ numCorners ],
                             const LocalCoordType &x,
                             const FieldType &factor,
                             JacobianTransposedType &J )
      {
        GlobalCoordType &q = J[ dim-1 ];
        bool affine;
        if( alwaysAffine )
        {
          const GlobalCoordType &top = TopMapping :: origin( coords );
          const GlobalCoordType &bottom = BottomMapping :: origin( coords );

          affine = BottomMapping :: Dphi_add( coords, x, factor, J );
          for( unsigned int i = 0; i < dimW; ++i )
            q[ i ] = factor * (top[ i ] - bottom[ i ]);
        }
        else
        {
          const FieldType xn = x[ dim-1 ];
          const FieldType cxn = FieldType( 1 ) - xn;
          const FieldType icxn = FieldType( 1 ) / cxn;
          LocalCoordType xb;
          for( unsigned int i = 0; i < dim-1; ++i )
            xb[ i ] = icxn * x[ i ];
          affine = BottomMapping :: Dphi_add( coords, xb, factor, J );

          TopMapping :: phi_add( coords, x, factor, q );
          BottomMapping :: phi_add( coords, xb, -factor, q );
          xb *= factor;
          for( unsigned int j = 0; j < dim-1; ++j )
          {
            for( unsigned int i = 0; i < dimW; ++i )
              q[ i ] += J[ j ][ i ] * xb[ j ];
          }
        }
        return affine;
      }

      static bool inDomain ( const LocalCoordType &x, FieldType factor )
      {
        const FieldType xn = x[ dim-1 ];
        const FieldType cxn = factor - xn;
        return (xn > -1e-12) && (cxn > -1e-12)
               && BottomMapping :: inDomain( x, factor * cxn );
      }

    };



    // CornerMapping
    // -------------


    template< class Topology, class CoordTraits >
    class CornerMapping
    {
      typedef CornerMapping< Topology, CoordTraits > ThisType;

    public:
      typedef MappingTraits< Topology :: dimension, CoordTraits > Traits;

      enum { dimG = Traits :: dimG };
      enum { dimW = Traits :: dimW };
      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianType JacobianType;
      typedef typename Traits :: JacobianTransposedType JacobianTransposedType;

      enum { numCorners = Topology :: numCorners };

      typedef GenericGeometry :: GenericMapping< Topology, Traits > GenericMapping;
      typedef GenericGeometry :: ReferenceElement< Topology, FieldType > ReferenceElement;

    protected:
      const GlobalCoordType *coords_[ numCorners ];

      mutable JacobianTransposedType jT_;
      mutable bool jTComputed;

    public:
      template< class CoordVector >
      explicit CornerMapping ( const CoordVector &coords )
        : jTComputed( false )
      {
        for( int i = 0; i < numCorners; ++i )
          coords_[ i ] = &(coords[ i ]);
      }

      bool affine () const
      {
        if( GenericMapping :: alwaysAffine )
          return true;
        jacobianT( baryCenter() );
        return jTComputed;
      }

      const GlobalCoordType &operator[] ( int i ) const
      {
        return *(coords_[ i ]);
      }

      int corners () const
      {
        return numCorners;
      }

      GlobalCoordType global( const LocalCoordType &x ) const
      {
        GlobalCoordType p;
        if( jTComputed )
        {
          MatrixHelper< CoordTraits > :: template ATx< dimG, dimW >( jT_, x, p );
          p += (*this)[ 0 ];
        }
        else
          GenericMapping :: phi_set( coords_, x, FieldType( 1 ), p );
        return p;
      }

      static bool checkInside ( const LocalCoordType &x )
      {
        return GenericMapping :: inDomain( x, FieldType( 1 ) );
      }

      const JacobianTransposedType &jacobianT ( const LocalCoordType &x ) const
      {
        if( !jTComputed )
        {
          jTComputed = GenericMapping :: Dphi_set( coords_, x, FieldType( 1 ), jT_ );
        }
        return jT_;
      }

    protected:
      static const LocalCoordType &baryCenter ()
      {
        return ReferenceElement :: template baryCenter< 0 >( 0 );
      }
    };

  }

}

#endif
