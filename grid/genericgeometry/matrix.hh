// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_MATRIXHELPER_HH
#define DUNE_GENERICGEOMETRY_MATRIXHELPER_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/static_assert.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    template< class Traits >
    struct MatrixHelper
    {
      typedef typename Traits :: ctype FieldType;

      template< int m, int n >
      static void
      Ax ( const typename Traits :: template Matrix< m, n > :: type &A,
           const typename Traits :: template Vector< n > :: type &x,
           typename Traits :: template Vector< m > :: type &ret )
      {
        for( int i = 0; i < m; ++i )
        {
          ret[ i ] = FieldType( 0 );
          for( int j = 0; j < n; ++j )
            ret[ i ] += A[ i ][ j ] * x[ j ];
        }
      }

      template< int m, int n >
      static void
      ATx ( const typename Traits :: template Matrix< m, n > :: type &A,
            const typename Traits :: template Vector< m > :: type &x,
            typename Traits :: template Vector< n > :: type &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          ret[ i ] = FieldType( 0 );
          for( int j = 0; j < m; ++j )
            ret[ i ] += A[ j ][ i ] * x[ j ];
        }
      }

      template< int m, int n, int p >
      static void
      AB ( const typename Traits :: template Matrix< m, n > :: type &A,
           const typename Traits :: template Matrix< n, p > :: type &B,
           typename Traits :: template Matrix< m, p > :: type &ret )
      {
        for( int i = 0; i < m; ++i )
        {
          for( int j = 0; j < p; ++j )
          {
            ret[ i ][ j ] = FieldType( 0 );
            for( int k = 0; k < n; ++k )
              ret[ i ][ j ] += A[ i ][ k ] * B[ k ][ j ];
          }
        }
      }

      template< int m, int n, int p >
      static void
      ATBT ( const typename Traits :: template Matrix< m, n > :: type &A,
             const typename Traits :: template Matrix< p, m > :: type &B,
             typename Traits :: template Matrix< n, p > :: type &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          for( int j = 0; j < p; ++j )
          {
            ret[ i ][ j ] = FieldType( 0 );
            for( int k = 0; k < m; ++k )
              ret[ i ][ j ] += A[ k ][ i ] * B[ j ][ k ];
          }
        }
      }

      template< int m, int n >
      static void
      ATA_L ( const typename Traits :: template Matrix< m, n > :: type &A,
              typename Traits :: template Matrix< n, n > :: type &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          for( int j = 0; j <= i; ++j )
          {
            ret[ i ][ j ] = FieldType( 0 );
            for( int k = 0; k < m; ++k )
              ret[ i ][ j ] += A[ k ][ i ] * A[ k ][ j ];
          }
        }
      }

      template< int m, int n >
      static void
      ATA ( const typename Traits :: template Matrix< m, n > :: type &A,
            typename Traits :: template Matrix< n, n > :: type &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          for( int j = 0; j <= i; ++j )
          {
            ret[ i ][ j ] = FieldType( 0 );
            for( int k = 0; k < m; ++k )
              ret[ i ][ j ] += A[ k ][ i ] * A[ k ][ j ];
            ret[ j ][ i ] = ret[ i ][ j ];
          }

          ret[ i ][ i ] = FieldType( 0 );
          for( int k = 0; k < m; ++k )
            ret[ i ][ i ] += A[ k ][ i ] * A[ k ][ i ];
        }
      }

      template< int m, int n >
      static void
      AAT_L ( const typename Traits :: template Matrix< m, n > :: type &A,
              typename Traits :: template Matrix< m, m > :: type &ret )
      {
        /*
           if (m==2) {
           ret[0][0] = A[0]*A[0];
           ret[1][1] = A[1]*A[1];
           ret[1][0] = A[0]*A[1];
           }
           else
         */
        for( int i = 0; i < m; ++i )
        {
          for( int j = 0; j <= i; ++j )
          {
            FieldType &retij = ret[ i ][ j ];
            retij = A[ i ][ 0 ] * A[ j ][ 0 ];
            for( int k = 1; k < n; ++k )
              retij += A[ i ][ k ] * A[ j ][ k ];
          }
        }
      }

      template< int m, int n >
      static void
      AAT ( const typename Traits :: template Matrix< m, n > :: type &A,
            typename Traits :: template Matrix< m, m > :: type &ret )
      {
        for( int i = 0; i < m; ++i )
        {
          for( int j = 0; j < i; ++j )
          {
            ret[ i ][ j ] = FieldType( 0 );
            for( int k = 0; k < n; ++k )
              ret[ i ][ j ] += A[ i ][ k ] * A[ j ][ k ];
            ret[ j ][ i ] = ret[ i ][ j ];
          }
          ret[ i ][ i ] = FieldType( 0 );
          for( int k = 0; k < n; ++k )
            ret[ i ][ i ] += A[ i ][ k ] * A[ i ][ k ];
        }
      }

      template< int n >
      static void
      Lx ( const typename Traits :: template Matrix< n, n > :: type &L,
           const typename Traits :: template Vector< n > :: type &x,
           typename Traits :: template Vector< n > :: type &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          ret[ i ] = FieldType( 0 );
          for( int j = 0; j <= i; ++j )
            ret[ i ] += L[ i ][ j ] * x[ j ];
        }
      }

      template< int n >
      static void
      LTx ( const typename Traits :: template Matrix< n, n > :: type &L,
            const typename Traits :: template Vector< n > :: type &x,
            typename Traits :: template Vector< n > :: type &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          ret[ i ] = FieldType( 0 );
          for( int j = i; j < n; ++j )
            ret[ i ] += L[ j ][ i ] * x[ j ];
        }
      }

      template< int n >
      static void
      LTL ( const typename Traits :: template Matrix< n, n > :: type &L,
            typename Traits :: template Matrix< n, n > :: type &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          for( int j = 0; j < i; ++j )
          {
            ret[ i ][ j ] = FieldType( 0 );
            for( int k = i; k < n; ++k )
              ret[ i ][ j ] += L[ k ][ i ] * L[ k ][ j ];
            ret[ j ][ i ] = ret[ i ][ j ];
          }
          ret[ i ][ i ] = FieldType( 0 );
          for( int k = i; k < n; ++k )
            ret[ i ][ i ] += L[ k ][ i ] * L[ k ][ i ];
        }
      }

      template< int n >
      static void
      LLT ( const typename Traits :: template Matrix< n, n > :: type &L,
            typename Traits :: template Matrix< n, n > :: type &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          for( int j = 0; j < i; ++j )
          {
            ret[ i ][ j ] = FieldType( 0 );
            for( int k = 0; k <= j; ++k )
              ret[ i ][ j ] += L[ i ][ k ] * L[ j ][ k ];
            ret[ j ][ i ] = ret[ i ][ j ];
          }
          ret[ i ][ i ] = FieldType( 0 );
          for( int k = 0; k <= i; ++k )
            ret[ i ][ i ] += L[ i ][ k ] * L[ i ][ k ];
        }
      }

      template< int n >
      static void
      cholesky_L ( const typename Traits :: template Matrix< n, n > :: type &A,
                   typename Traits :: template Matrix< n, n > :: type &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          FieldType &rii = ret[ i ][ i ];

          FieldType x = A[ i ][ i ];
          for( int j = 0; j < i; ++j )
            x -= ret[ i ][ j ] * ret[ i ][ j ];
          assert( x > FieldType( 0 ) );
          rii = sqrt( x );

          FieldType invrii = FieldType( 1 ) / rii;
          for( int k = i+1; k < n; ++k )
          {
            FieldType x = A[ k ][ i ];
            for( int j = 0; j < i; ++j )
              x -= ret[ i ][ j ] * ret[ k ][ j ];
            ret[ k ][ i ] = invrii * x;
          }
        }
      }

      template< int n >
      static FieldType
      detL ( const typename Traits :: template Matrix< n, n > :: type &L )
      {
        FieldType det = FieldType( 1 );
        for( int i = 0; i < n; ++i )
          det *= L[ i ][ i ];
        return det;
      }

      template< int n >
      static FieldType
      invL ( typename Traits :: template Matrix< n, n > :: type &L )
      {
        FieldType det = FieldType( 1 );
        for( int i = 0; i < n; ++i )
        {
          FieldType &lii = L[ i ][ i ];
          det *= lii;
          lii = FieldType( 1 ) / lii;
          for( int j = 0; j < i; ++j )
          {
            FieldType &lij = L[ i ][ j ];
            FieldType x = lij * L[ j ][ j ];
            for( int k = j+1; k < i; ++k )
              x += L[ i ][ k ] * L[ k ][ j ];
            lij = (-lii) * x;
          }
        }
        return det;
      }

      // calculates x := L^{-1} x
      template< int n >
      static void
      invLx ( typename Traits :: template Matrix< n, n > :: type &L,
              typename Traits :: template Vector< n > :: type &x )
      {
        for( int i = 0; i < n; ++i )
        {
          for( int j = 0; j < i; ++j )
            x[ i ] -= L[ i ][ j ] * x[ j ];
          x[ i ] /= L[ i ][ i ];
        }
      }

      // calculates x := L^{-T} x
      template< int n >
      static void
      invLTx ( typename Traits :: template Matrix< n, n > :: type &L,
               typename Traits :: template Vector< n > :: type &x )
      {
        for( int i = n; i > 0; --i )
        {
          for( int j = i; j < n; ++j )
            x[ i-1 ] -= L[ j ][ i-1 ] * x[ j ];
          x[ i-1 ] /= L[ i-1 ][ i-1 ];
        }
      }

      template< int n >
      static FieldType
      spdDetA ( const typename Traits :: template Matrix< n, n > :: type &A )
      {
        // return A[0][0]*A[1][1]-A[1][0]*A[1][0];
        typename Traits :: template Matrix< n, n > :: type L;
        cholesky_L< n >( A, L );
        return detL< n >( L );
      }

      template< int n >
      static FieldType
      spdInvA ( typename Traits :: template Matrix< n, n > :: type &A )
      {
        typename Traits :: template Matrix< n, n > :: type L;
        cholesky_L< n >( A, L );
        const FieldType det = invL< n >( L );
        LTL< n >( L, A );
        return det;
      }

      // calculate x := A^{-1} x
      template< int n >
      static void
      spdInvAx ( typename Traits :: template Matrix< n, n > :: type &A,
                 typename Traits :: template Vector< n > :: type &x )
      {
        typename Traits :: template Matrix< n, n > :: type L;
        cholesky_L< n >( A, L );
        invLx< n >( L, x );
        invLTx< n >( L, x );
      }

      template< int m, int n >
      static FieldType
      detATA ( const typename Traits :: template Matrix< m, n > :: type &A )
      {
        if( m >= n )
        {
          typename Traits :: template Matrix< n, n > :: type ata;
          ATA_L< m, n >( A, ata );
          return spdDetA< n >( ata );
        }
        else
          return FieldType( 0 );
      }

      template< int m, int n >
      static FieldType
      detAAT ( const typename Traits :: template Matrix< m, n > :: type &A )
      {
        if( n >= m )
        {
          typename Traits :: template Matrix< m, m > :: type aat;
          AAT_L< m, n >( A, aat );
          return spdDetA< m >( aat );
        }
        else
          return FieldType( 0 );
      }

      // A^{-1}_L = (A^T A)^{-1} A^T
      // => A^{-1}_L A = I
      template< int m, int n >
      static FieldType
      leftInvA ( const typename Traits :: template Matrix< m, n > :: type &A,
                 typename Traits :: template Matrix< n, m > :: type &ret )
      {
        dune_static_assert( (m >= n), "Matrix has no left inverse." );
        typename Traits :: template Matrix< n, n > :: type ata;
        ATA_L< m, n >( A, ata );
        const FieldType det = spdInvA< n >( ata );
        ATBT< n, n, m >( ata, A, ret );
        return det;
      }

      template< int m, int n >
      static void
      leftInvAx ( const typename Traits :: template Matrix< m, n > :: type &A,
                  const typename Traits :: template Vector< m > :: type &x,
                  typename Traits :: template Vector< n > :: type &y )
      {
        dune_static_assert( (m >= n), "Matrix has no left inverse." );
        typename Traits :: template Matrix< n, n > :: type ata;
        ATx< m, n >( A, x, y );
        ATA_L< m, n >( A, ata );
        spdInvAx< n >( ata, y );
      }

      // A^{-1}_R = A^T (A A^T)^{-1}
      // => A A^{-1}_R = I
      template< int m, int n >
      static FieldType
      rightInvA ( const typename Traits :: template Matrix< m, n > :: type &A,
                  typename Traits :: template Matrix< n, m > :: type &ret )
      {
        dune_static_assert( (n >= m), "Matrix has no right inverse." );
        if( (n == 2) && (m == 2) )
        {
          const FieldType det = (A[ 0 ][ 0 ]*A[ 1 ][ 1 ] - A[ 1 ][ 0 ]*A[ 0 ][ 1 ]);
          const FieldType detInv = FieldType( 1 ) / det;
          ret[ 0 ][ 0 ] = A[ 1 ][ 1 ] * detInv;
          ret[ 1 ][ 1 ] = A[ 0 ][ 0 ] * detInv;
          ret[ 1 ][ 0 ] = -A[ 1 ][ 0 ] * detInv;
          ret[ 0 ][ 1 ] = -A[ 0 ][ 1 ] * detInv;
          return det;
        }
        else
        {
          typename Traits :: template Matrix< m , m > :: type aat;
          AAT_L< m, n >( A, aat );
          const FieldType det = spdInvA< m >( aat );
          ATBT< m, n, m >( A , aat , ret );
          return det;
        }
      }

      template< int m, int n >
      static void
      xTRightInvA ( const typename Traits :: template Matrix< m, n > :: type &A,
                    const typename Traits :: template Vector< n > :: type &x,
                    typename Traits :: template Vector< m > :: type &y )
      {
        dune_static_assert( (n >= m), "Matrix has no right inverse." );
        typename Traits :: template Matrix< m, m > :: type aat;
        Ax< m, n >( A, x, y );
        AAT_L< m, n >( A, aat );
        spdInvAx< m >( aat, y );
      }
    };



    // Auxilliary Functions
    // --------------------

    template< class F, int dim >
    F two_norm2 ( const Dune::FieldVector< F, dim > &f )
    {
      F norm( 0 );
      for( int d = 0; d < dim; ++d )
        norm += f[ d ]*f[ d ];
      return norm;
    }

    template< class F, int rows, int cols >
    F frobenius_norm2 ( const Dune::FieldMatrix< F, rows, cols > &m )
    {
      F norm( 0 );
      for( int i = 0; i < rows; ++i )
        norm += two_norm2( m[ i ] );
      return norm;
    }

  }

}

#endif
