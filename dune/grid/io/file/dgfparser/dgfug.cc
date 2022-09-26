// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/geometry/utility/typefromvertexcount.hh>

#include <dune/grid/io/file/dgfparser/dgfug.hh>

namespace Dune
{

  namespace dgf
  {

    // Implementation of UGGridParameterBlock
    // --------------------------------------

    UGGridParameterBlock::UGGridParameterBlock ( std::istream &input )
      : GridParameterBlock( input ),
        noClosure_( false ), // default value
        noCopy_( true )   // default value
    {
      // check closure
      if( findtoken( "closure" ) )
      {
        std::string closure;
        if( getnextentry( closure ) )
        {
          makeupcase( closure );
          if( closure == "NONE" )
            noClosure_ = true ;
          else if( closure != "GREEN" )
          {
            dwarn << "UGGridParameterBlock: Parameter 'closure' has invalid value: " << closure
                  << ", using default: 'GREEN'." << std::endl;
          }
        }
      }
      else
      {
        dwarn << "UGGridParameterBlock: Parameter 'closure' not specified"
              << ", using default: 'GREEN'." << std::endl;
      }

      if( findtoken( "copies" ) )
      {
        std::string copies;
        if( getnextentry( copies ) )
        {
          makeupcase( copies );
          if( copies == "YES" )
            noCopy_ = false;
          else if( copies != "NO" )
          {
            dwarn << "UGGridParameterBlock: Parameter 'copies' has invalid value: " << copies
                  << ", using default: 'NO'." << std::endl;
          }
        }
      }
      else
      {
        dwarn << "UGGridParameterBlock: Parameter 'copies' not specified"
              << ", using default: 'NO'." << std::endl;
      }
    }

  } // namespace dgf



  // Implementation of DGFGridFactory< UGGrid >
  // -------------------------------------------

#if HAVE_DUNE_UGGRID
  template< int dim >
  void DGFGridFactory< UGGrid< dim > >::generate ( std::istream &input )
  {
    dgf_.element = DuneGridFormatParser::General;

    if( !dgf_.readDuneGrid( input, dim, dim ) )
      DUNE_THROW( DGFException, "Error: Failed to build grid");

    dgf_.setOrientation( 0, 1 );

    // get grid parameter block
    dgf::UGGridParameterBlock gridParam( input );

    // create grid factory (passed grid is returned by createGrid method)
    for( int n = 0; n < dgf_.nofvtx; n++ )
    {
      FieldVector< double, dim > v;
      for( int j = 0; j < dim; j++ )
        v[ j ] = dgf_.vtx[ n ][ j ];
      factory_.insertVertex( v );
    }

    std::vector< unsigned int > el;
    for( int n = 0; n < dgf_.nofelements; n++ )
    {
      el.clear();
      for( size_t j = 0; j < dgf_.elements[ n ].size(); ++j )
        el.push_back( ( dgf_.elements[ n ][ j ] ) );
      factory_.insertElement( geometryTypeFromVertexCount( dim, el.size() ), el );
    }

    // create grid
    grid_ = factory_.createGrid().release();

    // set closure type to none if parameter say so
    if( gridParam.noClosure() )
      grid_->setClosureType( UGGrid< dim >::NONE );

    if ( !gridParam.noCopy() )
      grid_->setRefinementType( UGGrid< dim >::COPY );
  }


  template void DGFGridFactory< UGGrid< 2 > >::generate ( std::istream &input );
  template void DGFGridFactory< UGGrid< 3 > >::generate ( std::istream &input );
#endif // #if HAVE_DUNE_UGGRID

} // namespace Dune
