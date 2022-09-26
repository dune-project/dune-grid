// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_GEOMETRYCACHE_HH
#define DUNE_ALBERTA_GEOMETRYCACHE_HH

#include <dune/grid/albertagrid/misc.hh>
#include <dune/grid/albertagrid/algebra.hh>

#if HAVE_ALBERTA

namespace Dune
{

  namespace Alberta
  {

    // GeometryCache
    // -------------

    template< int dim >
    class GeometryCache
    {
      static const unsigned int flagIntegrationElement = (1 << 0);
      static const unsigned int flagJacobianTransposed = (1 << 1);
      static const unsigned int flagJacobianInverseTransposed = (1 << 2);

    public:
      typedef FieldMatrix< Real, dimWorld, dim > JacobianInverseTransposed;
      typedef FieldMatrix< Real, dim, dimWorld > JacobianTransposed;

      GeometryCache ()
        : flags_( 0 )
      {}

      const Real &integrationElement ( const ALBERTA EL_INFO &elInfo )
      {
        if( (flags_ & flagIntegrationElement) == 0 )
        {
          integrationElement_ = std::abs( determinant( jacobianTransposed( elInfo ) ) );
          assert( integrationElement_ > 1e-14 );
          flags_ |= flagIntegrationElement;
        }
        return integrationElement_;
      }

      const JacobianTransposed &jacobianTransposed ( const ALBERTA EL_INFO &elInfo )
      {
        if( (flags_ & flagJacobianTransposed) == 0 )
        {
          assert( (elInfo.fill_flag & FillFlags< dim >::coords) != 0 );
          const GlobalVector &x = elInfo.coord[ 0 ];
          for( int i = 0; i < dim; ++i )
          {
            const GlobalVector &y = elInfo.coord[ i+1 ];
            for( int j = 0; j < dimWorld; ++j )
              jacobianTransposed_[ i ][ j ] = y[ j ] - x[ j ];
          }
          flags_ |= flagJacobianTransposed;
        }
        return jacobianTransposed_;
      }

      const JacobianInverseTransposed &
      jacobianInverseTransposed ( const ALBERTA EL_INFO &elInfo )
      {
        if( (flags_ & flagJacobianInverseTransposed) == 0 )
        {
          integrationElement_ = std::abs( invert( jacobianTransposed( elInfo ), jacobianInverseTransposed_ ) );
          assert( integrationElement_ > 1e-14 );
          flags_ |= flagIntegrationElement | flagJacobianInverseTransposed;
        }
        return jacobianInverseTransposed_;
      }

    private:
      unsigned int flags_;
      Real integrationElement_;
      FieldMatrix< Real, dim, dimWorld > jacobianTransposed_;
      FieldMatrix< Real, dimWorld, dim > jacobianInverseTransposed_;
    };



    // GeometryCacheProxy
    // ------------------

    template< int dim >
    struct GeometryCacheProxy
    {
      typedef FieldMatrix< Real, dimWorld, dim > JacobianInverseTransposed;
      typedef FieldMatrix< Real, dim, dimWorld > JacobianTransposed;

      GeometryCacheProxy ( GeometryCache< dim > &geometryCache, const ALBERTA EL_INFO &elInfo )
        : geometryCache_( geometryCache ),
          elInfo_( elInfo )
      {}

      const Real &integrationElement ()
      {
        return geometryCache_.integrationElement( elInfo_ );
      }

      const JacobianTransposed &jacobianTransposed ()
      {
        return geometryCache_.jacobianTransposed( elInfo_ );
      }

      const JacobianInverseTransposed &jacobianInverseTransposed ()
      {
        return geometryCache_.jacobianInverseTransposed( elInfo_ );
      }

    private:
      GeometryCache< dim > &geometryCache_;
      const ALBERTA EL_INFO &elInfo_;
    };

  } // namespace Alberta

} // namespace Dune

#endif // #if HAVE_ALBERTA

#endif // #ifndef DUNE_ALBERTA_GEOMETRYCACHE_HH
