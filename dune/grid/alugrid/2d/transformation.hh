// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_2D_TRANSFORMATION_HH
#define DUNE_ALUGRID_2D_TRANSFORMATION_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#if ENABLE_ALUGRID

namespace Dune
{

  template< int dimw >
  struct ALU2dTransformation
  {
    typedef alu2d_ctype ctype;

    static const int dimension = dimw;

    typedef FieldVector< ctype, dimension > WorldVector;
    typedef FieldMatrix< ctype, dimension, dimension > WorldMatrix;

    ALU2dTransformation ( const WorldMatrix &matrix, const WorldVector &shift )
      : matrix_( matrix ),
        shift_( shift )
    {}

    WorldVector evaluate ( const WorldVector &x ) const
    {
      WorldVector y = shift_;
      matrix_.umv( x, y );
      return y;
    }

    WorldVector evaluateInverse ( const WorldVector &y ) const
    {
      // Note: We assume the matrix to be orthogonal, here
      WorldVector ys = y - shift_;
      WorldVector x;
      matrix_.mtv( ys, x );
      return x;
    }

  private:
    WorldMatrix matrix_;
    WorldVector shift_;
  };

}

#endif // #if ENABLE_ALUGRID

#endif // #ifndef DUNE_ALUGRID_2D_TRANSFORMATION_HH
