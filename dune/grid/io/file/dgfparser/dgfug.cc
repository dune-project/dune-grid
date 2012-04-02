// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

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
        noCopy_( true ),   // default value
        heapSize_( 0 )     // default value (see UGGrid constructor)
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

      if( findtoken( "heapsize" ) )
      {
        int heapSize;
        if( getnextentry( heapSize ) )
        {
          if( heapSize > 0 )
            heapSize_ = heapSize;
          else
          {
            dwarn << "UGGridParameterBlock: Parameter 'heapsize' is non-positive"
                  << ", using default: '500' MB." << std::endl;
          }
        }
      }
      else
      {
        dwarn << "UGGridParameterBlock: Parameter 'heapsize' not specified"
              << ", using default: '500' MB." << std::endl;
      }
    }

  } // namespace dgf



  // Implementation of DGFGridFactory< UGGrid >
  // -------------------------------------------

#if ENABLE_UG
  template< int dim >
  void DGFGridFactory< UGGrid< dim > >::generate ( std::istream &input )
  {
    dgf_.element = DuneGridFormatParser::General;

    if( !dgf_.readDuneGrid( input, dim, dim ) )
      DUNE_THROW( DGFException, "Error: Failed to build grid");

    dgf_.setOrientation( 0, 1 );

    // get grid parameter block
    dgf::UGGridParameterBlock gridParam( input );

    // create grid here to set heap size
    // create grid factory (passed grid is returned by createGrid method)
    if( gridParam.heapSize() > 0 )
      UGGrid< dim >::setDefaultHeapSize( gridParam.heapSize() );

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

      // simplices
      if( el.size() == std::size_t( dim+1 ) )
        factory_.insertElement( GeometryType( GeometryType::simplex, dim ), el );
      // cubes
      else if( el.size() == 1u << dim )
        factory_.insertElement( GeometryType( GeometryType::cube, dim ), el );
#ifdef EXPERIMENTAL_GRID_EXTENSIONS
      // pyramid
      else if( (dim == 3) && (el.size() == 5u) )
        factory_.insertElement( GeometryType( GeometryType::pyramid, dim ), el );
      // prisms
      else if( (dim == 3) && (el.size() == 6u) )
        factory_.insertElement( GeometryType( GeometryType::prism, dim ), el );
#endif
      else
        DUNE_THROW( DGFException, "Invalid number of element vertices: " << el.size() );
    }

    // create grid
    grid_ = factory_.createGrid();

    // set closure type to none if parameter say so
    if( gridParam.noClosure() )
      grid_->setClosureType( UGGrid< dim >::NONE );

    if ( !gridParam.noCopy() )
      grid_->setRefinementType( UGGrid< dim >::COPY );
  }


  template void DGFGridFactory< UGGrid< 2 > >::generate ( std::istream &input );
  template void DGFGridFactory< UGGrid< 3 > >::generate ( std::istream &input );
#endif // #if ENABLE_UG

} // namespace Dune
