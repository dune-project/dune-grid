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



    // GenericCornerMapping
    // --------------------

    template< class Topology, class Traits, unsigned int offset = 0 >
    class GenericCornerMapping;

    template< class Traits, unsigned int offset >
    class GenericCornerMapping < Point, Traits, offset >
    {
      typedef Point Topology;

    public:
      static const unsigned int dim = Topology :: dimension;

      static const unsigned int dimW = Traits :: dimW;
      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianTransposedType JacobianTransposedType;

      static const bool alwaysAffine = true;

      template< class CoordStorage >
      static const GlobalCoordType &origin ( const CoordStorage &coords )
      {
        dune_static_assert( CoordStorage :: size, "Invalid offset." );
        return coords[ offset ];
      }

      template< class CoordStorage >
      static void phi_set ( const CoordStorage &coords,
                            const LocalCoordType &x,
                            const FieldType &factor,
                            GlobalCoordType &p )
      {
        const GlobalCoordType &y = origin( coords );
        for( unsigned int i = 0; i < dimW; ++i )
          p[ i ] = factor * y[ i ];
      }

      template< class CoordStorage >
      static void phi_add ( const CoordStorage &coords,
                            const LocalCoordType &x,
                            const FieldType &factor,
                            GlobalCoordType &p )
      {
        const GlobalCoordType &y = origin( coords );
        for( unsigned int i = 0; i < dimW; ++i )
          p[ i ] += factor * y[ i ];
      }

      template< class CoordStorage >
      static bool Dphi_set ( const CoordStorage &coords,
                             const LocalCoordType &x,
                             const FieldType &factor,
                             JacobianTransposedType &J )
      {
        return true;
      }

      template< class CoordStorage >
      static bool Dphi_add ( const CoordStorage &coords,
                             const LocalCoordType &x,
                             const FieldType &factor,
                             JacobianTransposedType &J )
      {
        return true;
      }
    };


    template< class BaseTopology, class Traits, unsigned int offset >
    class GenericCornerMapping< Prism< BaseTopology >, Traits, offset >
    {
      typedef Prism< BaseTopology > Topology;

      typedef GenericCornerMapping< BaseTopology, Traits, offset > BottomMapping;
      typedef GenericCornerMapping< BaseTopology, Traits, offset + BaseTopology :: numCorners >
      TopMapping;

    public:
      static const unsigned int dim = Topology :: dimension;

      static const unsigned int dimW = Traits :: dimW;
      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianTransposedType JacobianTransposedType;

      static const bool alwaysAffine = ((dim < 2) || Traits :: affine);

      template< class CoordStorage >
      static const GlobalCoordType &origin ( const CoordStorage &coords )
      {
        return BottomMapping :: origin( coords );
      }

      template< class CoordStorage >
      static void phi_set ( const CoordStorage &coords,
                            const LocalCoordType &x,
                            const FieldType &factor,
                            GlobalCoordType &p )
      {
        const FieldType xn = x[ dim-1 ];
        const FieldType cxn = FieldType( 1 ) - xn;
        BottomMapping :: phi_set( coords, x, factor * cxn, p );
        TopMapping :: phi_add( coords, x, factor * xn, p );
      }

      template< class CoordStorage >
      static void phi_add ( const CoordStorage &coords,
                            const LocalCoordType &x,
                            const FieldType &factor,
                            GlobalCoordType &p )
      {
        const FieldType xn = x[ dim-1 ];
        const FieldType cxn = FieldType( 1 ) - xn;
        BottomMapping :: phi_add( coords, x, factor * cxn, p );
        TopMapping :: phi_add( coords, x, factor * xn, p );
      }

      template< class CoordStorage >
      static bool Dphi_set ( const CoordStorage &coords,
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

      template< class CoordStorage >
      static bool Dphi_add ( const CoordStorage &coords,
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
    };


    template< class BaseTopology, class Traits, unsigned int offset >
    class GenericCornerMapping < Pyramid< BaseTopology >, Traits, offset >
    {
      typedef Pyramid< BaseTopology > Topology;

      typedef GenericCornerMapping< BaseTopology, Traits, offset > BottomMapping;
      typedef GenericCornerMapping< Point, Traits, offset + BaseTopology :: numCorners >
      TopMapping;

    public:
      static const unsigned int dim = Topology :: dimension;

      static const unsigned int dimW = Traits :: dimW;
      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianTransposedType JacobianTransposedType;

      static const bool alwaysAffine = (BottomMapping :: alwaysAffine || Traits :: affine);

      template< class CoordStorage >
      static const GlobalCoordType &origin ( const CoordStorage &coords )
      {
        return BottomMapping :: origin( coords );
      }

      template< class CoordStorage >
      static void phi_set ( const CoordStorage &coords,
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

      template< class CoordStorage >
      static void phi_add ( const CoordStorage &coords,
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

      template< class CoordStorage >
      static bool Dphi_set ( const CoordStorage &coords,
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

      template< class CoordStorage >
      static bool Dphi_add ( const CoordStorage &coords,
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
    };



    // SubMappingCoords
    // ----------------

    template< class Mapping, unsigned int codim >
    class SubMappingCoords
    {
      typedef typename Mapping :: GlobalCoordType GlobalCoordType;
      typedef typename Mapping :: ReferenceElement ReferenceElement;

      enum { dimension = ReferenceElement :: dimension };

      const Mapping &mapping_;
      const unsigned int i_;

    public:
      SubMappingCoords ( const Mapping &mapping, unsigned int i )
        : mapping_( mapping ), i_( i )
      {}

      const GlobalCoordType &operator[] ( unsigned int j ) const
      {
        const unsigned int k
          = ReferenceElement :: template subNumbering< codim, dimension - codim >( i_, j );
        return mapping_.corner( k );
      }
    };



    // CoordPointerStorage
    // -------------------

    template< class Topology, class Coordinate >
    class CoordPointerStorage
    {
      typedef CoordPointerStorage< Topology, Coordinate > This;

    public:
      static const unsigned int size = Topology :: numCorners;

      template< unsigned int codim, unsigned int i >
      struct SubTopology
      {
        typedef typename GenericGeometry :: SubTopology< Topology, codim, i > :: type type;
        typedef CoordPointerStorage< type, Coordinate > CornerStorage;
      };

    private:
      const Coordinate *coords_[ size ];

    public:
      template< class CoordVector >
      explicit CoordPointerStorage ( const CoordVector &coords )
      {
        for( unsigned int i = 0; i < size; ++i )
          coords_[ i ] = &(coords[ i ]);
      }

      const Coordinate &operator[] ( unsigned int i ) const
      {
        return *(coords_[ i ]);
      }
    };



    // CornerMapping
    // -------------

    template< class Topology, class CoordTraits, class CornerStorage >
    class CornerMapping
    {
      typedef CornerMapping< Topology, CoordTraits, CornerStorage > This;

    public:
      typedef MappingTraits< Topology :: dimension, CoordTraits > Traits;

      typedef CornerStorage CornerStorageType;

      static const unsigned int dimG = Traits :: dimG;
      static const unsigned int dimW = Traits :: dimW;

      typedef typename Traits :: FieldType FieldType;
      typedef typename Traits :: LocalCoordType LocalCoordType;
      typedef typename Traits :: GlobalCoordType GlobalCoordType;
      typedef typename Traits :: JacobianType JacobianType;
      typedef typename Traits :: JacobianTransposedType JacobianTransposedType;

      template< unsigned int codim, unsigned int i >
      struct SubTopology
      {
        typedef typename GenericGeometry :: SubTopology< Topology, codim, i > :: type type;

        typedef typename CornerStorage :: template SubTopology< codim, i > :: CornerStorage
        CornerStorageType;

        typedef CornerMapping< type, CoordTraits, CornerStorageType > Trace;

        typedef SubMappingCoords< This, codim > TraceCoordVector;
      };

    private:
      typedef GenericGeometry :: GenericCornerMapping< Topology, Traits > GenericMapping;

    public:
      static const bool alwaysAffine = GenericMapping :: alwaysAffine;

    protected:
      CornerStorageType coords_;

    public:
      template< class CoordVector >
      explicit CornerMapping ( const CoordVector &coords )
        : coords_( coords )
      {}

      const GlobalCoordType &corner ( int i ) const
      {
        return coords_[ i ];
      }

      void global ( const LocalCoordType &x, GlobalCoordType &ret ) const
      {
        GenericMapping :: phi_set( coords_, x, FieldType( 1 ), ret );
      }

      bool jacobianTransposed ( const LocalCoordType &x,
                                JacobianTransposedType &ret ) const
      {
        return GenericMapping :: Dphi_set( coords_, x, FieldType( 1 ), ret );
      }

      template< unsigned int codim, unsigned int i >
      typename SubTopology< codim, i > :: Trace trace () const
      {
        typedef typename SubTopology< codim, i > :: Trace Trace;
        typedef typename SubTopology< codim, i > :: TraceCoordVector CoordVector;
        return Trace( CoordVector( *this, i ) );
      }
    };

  }

}

#endif
