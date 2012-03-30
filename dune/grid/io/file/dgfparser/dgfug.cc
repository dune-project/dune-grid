// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
//- local includes
#include <dune/grid/io/file/dgfparser/dgfug.hh>


namespace Dune
{

  namespace dgf
  {

    // Implementation of UGGridParameterBlock
    // --------------------------------------

    UGGridParameterBlock::UGGridParameterBlock( std::istream &input )
      : GridParameterBlock( input ),
        _noClosure( false ), // default value
        _noCopy( true ),   // default value
        _heapsize( 0 )     // default value (see UGGrid constructor)
    {
      // check closure
      if (findtoken( "closure") )
      {
        std::string clo;
        if( getnextentry(clo) )
        {
          makeupcase(clo);
          if(clo == "NONE" )
          {
            _noClosure = true ;
          }
        }
      }
      else
      {
        dwarn << "UGGridParameterBlock: Parameter 'closure' not specified, "
              << "defaulting to 'GREEN'." << std::endl;
      }

      if (findtoken( "copies") )
      {
        std::string copies;
        if( getnextentry(copies) )
        {
          makeupcase(copies);
          if(copies == "YES" )
          {
            _noCopy = false ;
          }
        }
      }
      else
      {
        dwarn << "UGGridParameterBlock: Parameter 'copies' not specified, "
              << "no copies will be generated." << std::endl;
      }

      bool foundHeapSize = false ;
      if (findtoken( "heapsize") )
      {
        int heap;
        if( getnextentry( heap ) )
        {
          if( heap > 0 )
          {
            _heapsize = heap;
            foundHeapSize = true ;
          }
        }
      }

      if( ! foundHeapSize )
      {
        dwarn << "UGGridParameterBlock: Parameter 'heapsize' not specified, "
              << "defaulting to '500' MB." << std::endl;
      }
    }

    bool UGGridParameterBlock::noClosure () const
    {
      return _noClosure;
    }

    bool UGGridParameterBlock::noCopy () const
    {
      return _noCopy;
    }

    size_t UGGridParameterBlock::heapSize() const
    {
      return _heapsize;
    }

  } // end namespace dgf



  // Implementation of DGFGridFactory< UGGrid >
  // -------------------------------------------

  template< int dim >
  inline void DGFGridFactory< UGGrid< dim > >::generate ( std::istream &input )
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

    // eval 2^dim
    size_t two_power_dim = (1 << dim );

    std::vector< unsigned int > el;
    for( int n = 0; n < dgf_.nofelements; n++ )
    {
      el.clear();
      for( size_t j = 0; j < dgf_.elements[ n ].size(); j++ )
      {
        el.push_back( ( dgf_.elements[ n ][ j ] ) );
      }

      // simplices
      if( (int) el.size() == dgf_.dimw+1 )
        factory_.insertElement( GeometryType( GeometryType::simplex, dim ), el );
      // cubes
      else if( el.size() == two_power_dim )
        factory_.insertElement( GeometryType( GeometryType::cube, dim), el );
#ifdef EXPERIMENTAL_GRID_EXTENSIONS
      // pyramid
      else if( el.size() == 5 )
        factory_.insertElement( GeometryType( GeometryType::pyramid, dim ), el );
      // prisms
      else if( el.size() == 6 )
        factory_.insertElement( GeometryType( GeometryType::prism, dim ), el );
#endif
      else
        DUNE_THROW( DGFException, "Wrong number of vertices for element" );
    }

    // create grid
    grid_ = factory_.createGrid();

    // set closure type to none if parameter say so
    if( gridParam.noClosure() )
    {
      grid_->setClosureType( UGGrid< dim >::NONE );
    }

    if ( !gridParam.noCopy() )
    {
      grid_->setRefinementType( UGGrid< dim >::COPY );
    }
  }

} // end namespace Dune
