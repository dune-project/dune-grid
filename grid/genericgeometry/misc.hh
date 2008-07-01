// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_MISC_HH
#define DUNE_GENERICGEOMETRY_MISC_HH

namespace Dune
{

  namespace GenericGeometry
  {

    template< unsigned int n >
    struct Faculty
    {
      enum { value = n * Faculty< n-1 > :: value };
    };

    template<>
    struct Faculty< 0 >
    {
      enum { value = 1 };
    };



    template< bool condition, class TypeTrue, class TypeFalse >
    struct TypeIf;

    template< class TypeTrue, class TypeFalse >
    struct TypeIf< false, TypeTrue, TypeFalse >
    {
      typedef TypeFalse type;
    };

    template< class TypeTrue, class TypeFalse >
    struct TypeIf< true, TypeTrue, TypeFalse >
    {
      typedef TypeTrue type;
    };


    template< bool condition, template< bool > class True, template< bool > class False >
    struct ProtectedIf;

    template< template< bool > class True, template< bool > class False >
    struct ProtectedIf< true, True, False >
      : public True< true >
    {};

    template< template< bool > class True, template< bool > class False >
    struct ProtectedIf< false, True, False >
      : public False< false >
    {};



    // *******************************
    // Helpers
    // *******************************
    // * For the description of the traits class
    // * see mapping.hh
    template < class CoordTraits>
    struct QRDecompose
    {
      enum {dimW = CoordTraits :: dimW};
      enum {dimG = CoordTraits :: dimG};
      typedef typename CoordTraits :: field_type FieldType;
      typedef typename CoordTraits :: local_type LocalCoordType;
      typedef typename CoordTraits :: global_type GlobalCoordType;
      typedef typename CoordTraits :: derivative_type JacobianType;
      typedef typename CoordTraits :: localmapping_type SquareMappingType;
      typedef typename CoordTraits :: derivativeT_type JacobianTransposeType;
      typedef typename CoordTraits :: coord_vector CoordVector;

      static void
      multiplyTransposed ( JacobianTransposeType& A,
                           SquareMappingType &B,
                           JacobianType &C )
      {
        for( int i = 0; i < dimW; ++i )
        {
          LocalCoordType &row = C[ i ];
          for( int j = 0; j < dimG; ++j )
          {
            FieldType &field = row[ j ];

            field = FieldType( 0 );
            for( int k = 0; k < dimG; ++k )
              field += A[ k ][ i ] * B[ k ][ j ];
          }
        }
        return C;
      }

      static FieldType
      compute ( JacobianTransposeType& QT,
                SquareMappingType& R,
                JacobianType& inverseQ)
      {
        // Compute a LQ decomposition for a
        //   (non square) matrix d:
        //   d=LQ (d^T=Q^TR) where
        //   L is a square lower triangular matrix
        //   and QQ^T = I
        // In addition the integrationElement
        //   sqrt(det(dd^T) is computed using
        //   dd^T=LQ Q^TL^T = LL^T
        //   and det(dd^T)=det(L)det(L^T)=(sum(l_ii))^2
        //   we thus return |sum(l_ii)|
        FieldType integrationElement = FieldType( 1 );
        SquareMappingType RInvT (0) ;
        for( unsigned int i = 0; i < dimG; ++i )
        {
          for( unsigned int j = 0; j < i; ++j )
          {
            FieldType dot = QT[ i ] * QT[ j ];
            R[ j ][ i ] = dot;
            QT[ i ].axpy(  -dot, QT[ j ] );
          }

          FieldType norm = sqrt( QT[ i ] * QT[ i ] );
          integrationElement *= norm;

          LocalCoordType& RInvT_i = RInvT[ i ];
          RInvT_i[ i ] = (FieldType( 1. ) / norm);
          QT[ i ] *= RInvT_i[ i ];

          for( unsigned int j = i; j > 0; --j )
          {
            for( unsigned int k = j; k <= i; ++k )
              RInvT_i[ j-1 ] -= R[ j-1 ][ k ] * RInvT_i[ k ];
            RInvT_i[ j-1 ] *= RInvT[ j-1 ][ j-1 ];
          }
        }

        multiplyTransposed( QT, RInvT, inverseQ );

        return integrationElement;
      }
    };

  }
}

#endif
