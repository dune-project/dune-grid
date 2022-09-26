// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/io/file/dgfparser/blocks/periodicfacetrans.hh>

namespace Dune
{

  namespace dgf
  {

    // PeriodicFaceTransformationBlock
    // -------------------------------

    PeriodicFaceTransformationBlock
    ::PeriodicFaceTransformationBlock ( std::istream &in, int dimworld )
      : BasicBlock( in, "PeriodicFaceTransformation" )
    {
      while( getnextline() )
      {
        AffineTransformation trafo( dimworld );
        for( int i = 0; i < dimworld; ++i )
        {
          if( i > 0 )
            match( ',' );

          for( int j = 0; j < dimworld; ++j )
          {
            if( !getnextentry( trafo.matrix( i, j ) ) )
            {
              DUNE_THROW( DGFException,
                          "Error in " << *this << ": "
                                      << "Not enough entries in matrix row " << i << "." );
            }
          }
        }

        match( '+' );
        for( int i = 0; i < dimworld; ++i )
        {
          if( !getnextentry( trafo.shift[ i ] ) )
          {
            DUNE_THROW( DGFException,
                        "Error in " << *this << ": "
                                    << "Not enough entries in shift." );
          }
        }

        transformations_.push_back( trafo );
      }
    }


    void PeriodicFaceTransformationBlock::match ( char what )
    {
      char c;
      if( !getnextentry( c ) || (c != what) )
      {
        DUNE_THROW( DGFException,
                    "Error in " << *this << ": " << what << "expected." );
      }
    }

  } // end namespace dgf

} // end namespace Dune
