// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_ALGEBRA_HH
#define DUNE_ALBERTA_ALGEBRA_HH

#include <dune/common/misc.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dune
{

  namespace Alberta
  {

    template< class K >
    inline static FieldVector< K, 3 >
    vectorProduct ( const FieldVector< K, 3 > &u, const FieldVector< K, 3 > &v )
    {
      FieldVector< K, 3 > w;
      w[ 0 ] = u[ 1 ] * v[ 2 ] - u[ 2 ] * v[ 1 ];
      w[ 1 ] = u[ 2 ] * v[ 0 ] - u[ 0 ] * v[ 2 ];
      w[ 2 ] = u[ 0 ] * v[ 1 ] - u[ 1 ] * v[ 0 ];
      return w;
    }


    template< class K, int m >
    inline static K determinant ( const FieldMatrix< K, 0, m > &matrix )
    {
      return K( 1 );
    }

    template< class K >
    inline static K determinant ( const FieldMatrix< K, 1, 1 > &matrix )
    {
      return matrix[ 0 ][ 0 ];
    }

    template< class K, int m >
    inline static K determinant ( const FieldMatrix< K, 1, m > &matrix )
    {
      K sum = SQR( matrix[ 0 ][ 0 ] );
      for( int i = 1; i < m; ++i )
        sum += SQR( matrix[ 0 ][ i ] );
      return sqrt( sum );
    }

    template< class K >
    inline static K determinant ( const FieldMatrix< K, 2, 2 > &matrix )
    {
      return matrix[ 0 ][ 0 ] * matrix[ 1 ][ 1 ] - matrix[ 0 ][ 1 ] * matrix[ 1 ][ 0 ];
    }

    template< class K >
    inline static K determinant ( const FieldMatrix< K, 2, 3 > &matrix )
    {
      return vectorProduct( matrix[ 0 ], matrix[ 1 ] ).two_norm();
    }

    template< class K, int m >
    inline static K determinant ( const FieldMatrix< K, 2, m > &matrix )
    {
      const K tmpA = SQR( matrix[ 0 ] );
      const K tmpB = SQR( matrix[ 1 ] );
      const K tmpC = matrix[ 0 ] * matrix[ 1 ];
      return sqrt( tmpA * tmpB - SQR( tmpC ) );
    }

    template< class K >
    inline static K determinant ( const FieldMatrix< K, 3, 3 > &matrix )
    {
      return matrix[ 0 ] * vectorProduct( matrix[ 1 ], matrix[ 2 ] );
    }


    template< class K, int m >
    inline static K invert ( const FieldMatrix< K, 0, m > &matrix,
                             FieldMatrix< K, m, 0 > &inverse )
    {
      return K( 1 );
    }

    template< class K >
    inline static K invert ( const FieldMatrix< K, 1, 1 > &matrix,
                             FieldMatrix< K, 1, 1 > &inverse )
    {
      inverse[ 0 ][ 0 ] = K( 1 ) / matrix[ 0 ][ 0 ];
      return matrix[ 0 ][ 0 ];
    }

    template< class K, int m >
    inline static K invert ( const FieldMatrix< K, 1, m > &matrix,
                             FieldMatrix< K, m, 1 > &inverse )
    {
      K detSqr = SQR( matrix[ 0 ] );
      K invDetSqr = K( 1 ) / detSqr;
      for( int i = 0; i < m; ++i )
        inverse[ i ][ 0 ] = invDetSqr * matrix[ 0 ][ i ];
      return sqrt( detSqr );
    }

    template< class K >
    inline static K invert ( const FieldMatrix< K, 2, 2 > &matrix,
                             FieldMatrix< K, 2, 2 > &inverse )
    {
      K det = determinant( matrix );
      K invDet = K( 1 ) / det;
      inverse[ 0 ][ 0 ] =   invDet * matrix[ 1 ][ 1 ];
      inverse[ 0 ][ 1 ] = - invDet * matrix[ 0 ][ 1 ];
      inverse[ 1 ][ 0 ] = - invDet * matrix[ 1 ][ 0 ];
      inverse[ 1 ][ 1 ] =   invDet * matrix[ 0 ][ 0 ];
      return det;
    }

    template< class K, int m >
    inline static K invert ( const FieldMatrix< K, 2, m > &matrix,
                             FieldMatrix< K, m, 2 > &inverse )
    {
      const K tmpA = SQR( matrix[ 0 ] );
      const K tmpB = SQR( matrix[ 1 ] );
      const K tmpC = matrix[ 0 ] * matrix[ 1 ];
      const K detSqr = tmpA * tmpB - SQR( tmpC );
      const K invDetSqr = K( 1 ) / detSqr;
      for( int i = 0; i < m; ++i )
      {
        inverse[ i ][ 0 ] = invDetSqr * (tmpB * matrix[ 0 ][ i ] - tmpC * matrix[ 1 ][ i ]);
        inverse[ i ][ 1 ] = invDetSqr * (tmpA * matrix[ 1 ][ i ] - tmpC * matrix[ 0 ][ i ]);
      }
      return sqrt( detSqr );
    }

    template< class K >
    inline static K invert ( const FieldMatrix< K, 3, 3 > &matrix,
                             FieldMatrix< K, 3, 3 > &inverse )
    {
      return FMatrixHelp::invertMatrix( matrix, inverse );
    }
  }

}

#endif // #ifndef DUNE_ALBERTA_ALGEBRA_HH
