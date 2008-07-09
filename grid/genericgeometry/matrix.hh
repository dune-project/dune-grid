// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GENERICGEOMETRY_MATRIXHELPER_HH
#define DUNE_GENERICGEOMETRY_MATRIXHELPER_HH

namespace Dune
{

  namespace GenericGeometry
  {

    template< class Traits >
    struct MatrixHelper
    {
      typedef typename Traits :: field_type FieldType;

      template< int m, int n >
      static void
      Ax ( const typename Traits :: template Matrix< m, n > :: Type &A,
           const typename Traits :: template Vector< n > :: Type &x,
           typename Traits :: template Vector< m > :: Type &ret )
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
      ATx ( const typename Traits :: template Matrix< m, n > :: Type &A,
            const typename Traits :: template Vector< m > :: Type &x,
            typename Traits :: template Vector< n > :: Type &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          ret[ i ] = FieldType( 0 );
          for( int j = 0; j < m; ++j )
            ret[ i ] += A[ j ][ i ] * x[ j ];
        }
      }

      template< int m, int n, int p>
      static void
      AB ( const typename Traits :: template Matrix< m, n > :: Type &A,
           const typename Traits :: template Matrix< n, p > :: Type &B,
           typename Traits :: template Matrix< m, p > :: Type &AB )
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

      template< int m, int n >
      static void
      ATA ( const typename Traits :: typename Matrix< m, n > :: Type &A,
            typename Traits :: template Matrix< n, n > :: Type &ret )
      {
        for( int i = 0; i < n; ++i )
        {
          for( int j = 0; j <= i; ++j )
          {
            ret[ i ][ j ] = FieldType( 0 );
            for( int k = 0; k < m; ++k )
              ret[ i ][ j ] += A[ i ][ k ] * A[ k ][ j ];
          }
          for( int j = 0; j < i; ++j )
            ret[ j ][ i ] = ret[ i ][ j ];
        }
      }

      template< int m, int n >
      static FieldType
      invAT ( const typename Traits :: template Matrix< m, n > :: Type &A,
              typename Traits :: template Matrix< m, n > :: Type &ret )
      {
        typename Traits :: template Matrix< n, n > ATA;
        MatrixHelper :: ATA< m, n >( A, ATA );
        const FieldType det = MatrixHelper :: invP< n > ( ATA );
        MatrixHelper :: AB< m, n, n >( A, ATA, ret );
        return det;
      }

      template< int n >
      static FieldType
      invPD ( typename Traits :: typename Matrix< n, n > :: Type &A )
      {
        FieldType det = FieldType( 1 );
        // in place LU decomposition
        for( int j = 0; j < n-1; ++j )
        {
          FieldType &ajj = A[ j ][ j ];
          FieldType invajj = FieldType( 1 ) / ajj;
          det *= ajj;
          for( int i = j+1; i < n; ++i )
          {
            FieldType &aij = A[ i ][ j ];
            aij *= invajj;
            for( int k = j+1; k < n; ++k )
              A[ i ][ k ] -= aij * A[ j ][ k ];
          }
        }
        // invert L in place
        for( int i = 1; i < n; ++i )
        {
          for( int j = 0; j < i; ++j )
          {
            FieldType &aij = A[ i ][ j ];
            for( int k = j+1; k < i; ++k )
              aij += A[ i ][ k ] * A[ k ][ j ];
            aij *= -1;
          }
        }
        // invert U in place
        for( int j = 1; j < n; ++j )
        {
          FieldType &ajj = A[ j ][ j ];
          FieldType ajj = FieldType( 1 ) / ajj;
          for( int i = 0; i < j; ++i )
          {
            FieldType &aij = A[ i ][ j ];
            for( int k = i+1; k < j; ++k )
              aij += A[ i ][ k ] * A[ k ][ j ];
            aij *= -ajj;
          }
        }
        // in place multiplication
        for( int j = n; j > 0; --j )
        {
          for( int i = n; i > 0; --i )
          {
            FieldType v = FieldType( 0 );
            for( int k = 0; k < std :: min( i, j ); ++k )
              v += A[ i-1 ][ k ] * A[ k ][ j-1 ];
            A[ i-1 ][ j-1 ] = v;
          }
          return det;
        }

      };

    }

  }

#endif
