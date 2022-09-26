// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_TRANSFORMATION_HH
#define DUNE_ALBERTA_TRANSFORMATION_HH

#include <dune/common/fvector.hh>

#include <dune/grid/albertagrid/misc.hh>

#if HAVE_ALBERTA

namespace Dune
{

  class AlbertaTransformation
  {
    typedef Alberta::GlobalSpace GlobalSpace;

  public:
    typedef Alberta::Real ctype;

    static const int dimension = Alberta::dimWorld;

    typedef FieldVector< ctype, dimension > WorldVector;

    explicit
    AlbertaTransformation ( const Alberta::AffineTransformation *trafo = NULL )
      : matrix_( (trafo != NULL ? trafo->M : GlobalSpace::identityMatrix()) ),
        shift_( (trafo != NULL ? trafo->t : GlobalSpace::nullVector()) )
    {}

    AlbertaTransformation ( const GlobalSpace::Matrix &matrix,
                            const GlobalSpace::Vector &shift )
      : matrix_( matrix ),
        shift_( shift )
    {}

    WorldVector evaluate ( const WorldVector &x ) const
    {
      WorldVector y;
      for( int i = 0; i < dimension; ++i )
      {
        const GlobalSpace::Vector &row = matrix_[ i ];
        y[ i ] = shift_[ i ];
        for( int j = 0; j < dimension; ++j )
          y[ i ] += row[ j ] * x[ j ];
      }
      return y;
    }

    WorldVector evaluateInverse ( const WorldVector &y ) const
    {
      // Note: ALBERTA requires the matrix to be orthogonal
      WorldVector x( ctype( 0 ) );
      for( int i = 0; i < dimension; ++i )
      {
        const GlobalSpace::Vector &row = matrix_[ i ];
        const ctype v = y[ i ] - shift_[ i ];
        for( int j = 0; j < dimension; ++j )
          x[ j ] += row[ j ] * v;
      }
      return x;
    }

  private:
    const GlobalSpace::Matrix &matrix_;
    const GlobalSpace::Vector &shift_;
  };

}

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_TRANSFORMATION_HH
