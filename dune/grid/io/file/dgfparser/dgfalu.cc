// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**\cond */
namespace Dune
{

  namespace dgf
  {

    struct ALU2dGridParameterBlock
      : public GridParameterBlock
    {
      ALU2dGridParameterBlock( std::istream &in, const bool verbose )
        : GridParameterBlock( in ),
          tolerance_( 1e-8 )
      {
        if( findtoken( "tolerance" ) )
        {
          double x;
          if( getnextentry(x) )
            tolerance_ = x;
          else
          {
            if( verbose )
            {
              dwarn << "GridParameterBlock: found keyword `tolerance' but no value, "
                    << "defaulting to `" <<  tolerance_ <<"'!" << std::endl;
            }
          }
          if( tolerance_ <= 0 )
            DUNE_THROW( DGFException, "Nonpositive tolerance specified!" );
        }
        else
        {
          if( verbose )
          {
            dwarn << "GridParameterBlock: Parameter 'tolerance' not specified, "
                  << "defaulting to `" <<  tolerance_ <<"'!" << std::endl;
          }
        }
      }

      double tolerance () const { return tolerance_; }

    protected:
      double tolerance_;
    };

  }

  template < class G >
  inline bool DGFBaseFactory< G > ::
  generateALUGrid( const ALUGridElementType eltype,
                   std::istream &file, MPICommunicatorType communicator,
                   const std::string &filename )
  {
    typedef G DGFGridType ;

    const int dimworld = DGFGridType :: dimensionworld ;
    dgf_.element = ( eltype == simplex) ?
                   DuneGridFormatParser::Simplex :
                   DuneGridFormatParser::Cube ;
    dgf_.dimgrid = dimworld;
    dgf_.dimw = dimworld;

    const bool isDGF = dgf_.isDuneGridFormat( file );
    file.seekg( 0 );
    if( !isDGF )
      return false;

    int rank = 0;
#if ALU3DGRID_PARALLEL
    MPI_Comm_rank( communicator, &rank );
#endif

    dgf::GridParameterBlock parameter( file );

    typedef FieldVector< typename DGFGridType :: ctype, dimworld > CoordinateType ;

    if( rank == 0 )
    {
      if( !dgf_.readDuneGrid( file, dimworld, dimworld ) )
        DUNE_THROW( InvalidStateException, "DGF file not recognized on second call." );

      if( eltype == simplex )
      {
        dgf_.setOrientation( 2, 3 );
      }

      for( int n = 0; n < dgf_.nofvtx; ++n )
      {
        CoordinateType pos;
        for( int i = 0; i < dimworld; ++i )
          pos[ i ] = dgf_.vtx[ n ][ i ];
        factory_.insertVertex( pos );
      }

      GeometryType elementType( (eltype == simplex) ?
                                GeometryType::simplex :
                                GeometryType::cube, dimworld );

      const int nFaces = (eltype == simplex) ? dimworld+1 : 2*dimworld;
      for( int n = 0; n < dgf_.nofelements; ++n )
      {
        factory_.insertElement( elementType, dgf_.elements[ n ] );
        for( int face = 0; face <nFaces; ++face )
        {
          typedef DuneGridFormatParser::facemap_t::key_type Key;
          typedef DuneGridFormatParser::facemap_t::iterator Iterator;

          const Key key = ElementFaceUtil::generateFace( dimworld, dgf_.elements[ n ], face );
          const Iterator it = dgf_.facemap.find( key );
          if( it != dgf_.facemap.end() )
            factory_.insertBoundary( n, face, it->second.first );
        }
      }

      dgf::ProjectionBlock projectionBlock( file, dimworld );
      const DuneBoundaryProjection< dimworld > *projection
        = projectionBlock.defaultProjection< dimworld >();

      if( projection != 0 )
        factory_.insertBoundaryProjection( *projection );

      const size_t numBoundaryProjections = projectionBlock.numBoundaryProjections();
      for( size_t i = 0; i < numBoundaryProjections; ++i )
      {
        GeometryType type( (eltype == simplex) ?
                           GeometryType::simplex :
                           GeometryType::cube,
                           dimworld-1);

        const std::vector< unsigned int > &vertices = projectionBlock.boundaryFace( i );
        const DuneBoundaryProjection< dimworld > *projection
          = projectionBlock.boundaryProjection< dimworld >( i );
        factory_.insertBoundaryProjection( type, vertices, projection );
      }
    }

    typedef dgf::PeriodicFaceTransformationBlock::AffineTransformation Transformation;
    dgf::PeriodicFaceTransformationBlock trafoBlock( file, dimworld );
    const int size = trafoBlock.numTransformations();
    for( int k = 0; k < size; ++k )
    {
      const Transformation &trafo = trafoBlock.transformation( k );

      typename GridFactory::WorldMatrix matrix;
      for( int i = 0; i < dimworld; ++i )
        for( int j = 0; j < dimworld; ++j )
          matrix[ i ][ j ] = trafo.matrix( i, j );

      typename GridFactory::WorldVector shift;
      for( int i = 0; i < dimworld; ++i )
        shift[ i ] = trafo.shift[ i ];

      factory_.insertFaceTransformation( matrix, shift );
    }

    if ( ! parameter.dumpFileName().empty() )
      grid_ = factory_.createGrid( dgf_.facemap.empty(), false, parameter.dumpFileName() );
    else
      grid_ = factory_.createGrid( dgf_.facemap.empty(), true, filename );
    return true;
  }

  template <ALUGridElementType eltype, ALUGridRefinementType refinementtype, class Comm>
  inline bool DGFGridFactory< ALUGrid< 3, 3, eltype, refinementtype, Comm > >
  ::generate( std::istream &file, MPICommunicatorType communicator, const std::string &filename )
  {
    return BaseType :: generateALUGrid( eltype, file, communicator, filename );
  }


  // ALUGrid 2d version
  //-------------------

  template < class G >
  inline bool DGFBaseFactory< G > ::
  generateALU2dGrid( const ALUGridElementType eltype,
                     std::istream &file,
                     MPICommunicatorType communicator,
                     const std::string &filename )
  {
    const int dimgrid = 2;
    const int dimworld = G :: dimensionworld ;
    dgf_.element = (eltype == simplex) ?
                   DuneGridFormatParser::Simplex : DuneGridFormatParser::Cube ;
    dgf_.dimgrid = dimgrid;
    dgf_.dimw = dimworld;

    const bool isDGF = dgf_.isDuneGridFormat( file );
    file.seekg( 0 );
    if( !isDGF )
      return false;

    int rank = 0;
#if HAVE_MPI
    MPI_Comm_rank( communicator, &rank );
#endif

    // set verbosity of factory only for rank = 0
    factory_.setVerbosity( (rank == 0) );

    // only print warnings of ALU2dGridParameterBlock on rank = 0
    dgf::ALU2dGridParameterBlock parameter( file, (rank == 0) );

#if ALU2DGRID_PARALLEL
    // for a parallel ALUGrid implementation do the following only on rank 0
    if( rank == 0 )
#endif
    {
      factory_.setTolerance( parameter.tolerance() );

      if( !dgf_.readDuneGrid( file, dimgrid, dimworld ) )
        DUNE_THROW( InvalidStateException, "DGF file not recognized on second call." );

      for( int n = 0; n < dgf_.nofvtx; ++n )
      {
        FieldVector< double, dimworld > pos;
        for( int i = 0; i < dimworld; ++i )
          pos[ i ] = dgf_.vtx[ n ][ i ];
        factory_.insertVertex( pos );
      }

      GeometryType elementType( (eltype == simplex) ?
                                GeometryType::simplex :
                                GeometryType::cube, dimgrid );

      const int nFaces = (eltype == simplex) ? dimgrid+1 : 2*dimgrid;
      for( int n = 0; n < dgf_.nofelements; ++n )
      {
        factory_.insertElement( elementType, dgf_.elements[ n ] );
        for( int face = 0; face <nFaces; ++face )
        {
          typedef typename DuneGridFormatParser::facemap_t::key_type Key;
          typedef typename DuneGridFormatParser::facemap_t::iterator Iterator;

          const Key key = ElementFaceUtil::generateFace( dimgrid, dgf_.elements[ n ], face );
          const Iterator it = dgf_.facemap.find( key );
          if( it != dgf_.facemap.end() )
            factory_.insertBoundary( n, face, it->second.first );
        }
      }

      dgf::ProjectionBlock projectionBlock( file, dimworld );
      const DuneBoundaryProjection< dimworld > *projection
        = projectionBlock.defaultProjection< dimworld >();
      if( projection != 0 )
        factory_.insertBoundaryProjection( *projection );
      const size_t numBoundaryProjections = projectionBlock.numBoundaryProjections();
      for( size_t i = 0; i < numBoundaryProjections; ++i )
      {
        GeometryType type( (eltype == simplex) ?
                           GeometryType::simplex :
                           GeometryType::cube, dimgrid-1 );
        const std::vector< unsigned int > &vertices = projectionBlock.boundaryFace( i );
        const DuneBoundaryProjection< dimworld > *projection
          = projectionBlock.boundaryProjection< dimworld >( i );
        factory_.insertBoundaryProjection( type, vertices, projection );
      }
    }

    typedef dgf::PeriodicFaceTransformationBlock::AffineTransformation Transformation;
    dgf::PeriodicFaceTransformationBlock trafoBlock( file, dimworld );
    const int size = trafoBlock.numTransformations();
    for( int k = 0; k < size; ++k )
    {
      const Transformation &trafo = trafoBlock.transformation( k );

      typename GridFactory::WorldMatrix matrix;
      for( int i = 0; i < dimworld; ++i )
        for( int j = 0; j < dimworld; ++j )
          matrix[ i ][ j ] = trafo.matrix( i, j );

      typename GridFactory::WorldVector shift;
      for( int i = 0; i < dimworld; ++i )
        shift[ i ] = trafo.shift[ i ];

      factory_.insertFaceTransformation( matrix, shift );
    }

    if ( ! parameter.dumpFileName().empty() && (rank == 0) )
      grid_ = factory_.createGrid( dgf_.facemap.empty(), false, parameter.dumpFileName() );
    else
      grid_ = factory_.createGrid( dgf_.facemap.empty(), true, filename );
    return true;
  }

  template <int dimw, ALUGridElementType eltype,
      ALUGridRefinementType refinementtype, class Comm >
  inline bool DGFGridFactory< ALUGrid< 2, dimw, eltype, refinementtype, Comm > >
  ::generate( std::istream &file, MPICommunicatorType communicator, const std::string &filename )
  {
    return BaseType :: generateALU2dGrid( eltype, file, communicator, filename );
  }

} // end namespace Dune
  /** \endcond */
