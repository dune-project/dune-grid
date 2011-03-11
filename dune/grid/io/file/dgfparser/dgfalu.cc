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
      ALU2dGridParameterBlock( std::istream &in )
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
            dwarn << "GridParameterBlock: found keyword `tolerance' but no value, "
                  << "defaulting to `" <<  tolerance_ <<"'!" << std::endl;
          }
          if( tolerance_ <= 0 )
            DUNE_THROW( DGFException, "Nonpositive tolerance specified!" );
        }
        else
        {
          dwarn << "GridParameterBlock: Parameter 'tolerance' not specified, "
                << "defaulting to `" <<  tolerance_ <<"'!" << std::endl;
        }
      }

      double tolerance () const { return tolerance_; }

    protected:
      double tolerance_;
    };

  }



  inline bool DGFGridFactory< ALUSimplexGrid< 3, 3 > >
  ::generate( std::istream &file, MPICommunicatorType communicator, const std::string &filename )
  {
    const int dimworld = 3;
    dgf_.element = DuneGridFormatParser::Simplex;
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

    if( rank == 0 )
    {
      if( !dgf_.readDuneGrid( file, dimworld, dimworld ) )
        DUNE_THROW( InvalidStateException, "DGF file not recognized on second call." );

      dgf_.setOrientation( 2, 3 );

      for( int n = 0; n < dgf_.nofvtx; ++n )
      {
        FieldVector< ALUSimplexGrid< 3, 3 > :: ctype, dimworld > pos;
        for( int i = 0; i < dimworld; ++i )
          pos[ i ] = dgf_.vtx[ n ][ i ];
        factory_.insertVertex( pos );
      }

      GeometryType elementType( GeometryType::simplex, dimworld );
      for( int n = 0; n < dgf_.nofelements; ++n )
      {
        factory_.insertElement( elementType, dgf_.elements[ n ] );
        for( int face = 0; face <= dimworld; ++face )
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
        GeometryType type( GeometryType::simplex, dimworld-1 );
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

      GridFactory::WorldMatrix matrix;
      for( int i = 0; i < dimworld; ++i )
        for( int j = 0; j < dimworld; ++j )
          matrix[ i ][ j ] = trafo.matrix( i, j );

      GridFactory::WorldVector shift;
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


  inline bool DGFGridFactory< ALUCubeGrid< 3, 3 > >
  ::generate( std::istream &file, MPICommunicatorType communicator, const std::string &filename )
  {
    const int dimworld = 3;
    dgf_.element = DuneGridFormatParser::Cube;
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

    if( rank == 0 )
    {
      if( !dgf_.readDuneGrid( file, dimworld, dimworld ) )
        DUNE_THROW( InvalidStateException, "DGF file not recognized on second call." );

      for( int n = 0; n < dgf_.nofvtx; ++n )
      {
        FieldVector< ALUCubeGrid< 3, 3 >::ctype , dimworld > pos;
        for( int i = 0; i < dimworld; ++i )
          pos[ i ] = dgf_.vtx[ n ][ i ];
        factory_.insertVertex( pos );
      }

      GeometryType elementType( GeometryType::cube, dimworld );
      for( int n = 0; n < dgf_.nofelements; ++n )
      {
        factory_.insertElement( elementType, dgf_.elements[ n ] );
        for( int face = 0; face < 2*dimworld; ++face )
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
        GeometryType type( GeometryType::cube, dimworld-1 );
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

      GridFactory::WorldMatrix matrix;
      for( int i = 0; i < dimworld; ++i )
        for( int j = 0; j < dimworld; ++j )
          matrix[ i ][ j ] = trafo.matrix( i, j );

      GridFactory::WorldVector shift;
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


  template <int dimw>
  inline bool DGFGridFactory< ALUSimplexGrid< 2, dimw > >
  ::generate( std::istream &file, MPICommunicatorType communicator, const std::string &filename )
  {
    const int dimgrid = 2;
    const int dimworld = dimw;
    dgf_.element = DuneGridFormatParser::Simplex;
    dgf_.dimgrid = dimgrid;
    dgf_.dimw = dimworld;

    const bool isDGF = dgf_.isDuneGridFormat( file );
    file.seekg( 0 );
    if( !isDGF )
      return false;

    int rank = 0;
#if ALU2DGRID_PARALLEL
    MPI_Comm_rank( communicator, &rank );
#endif

    dgf::ALU2dGridParameterBlock parameter( file );

    if( rank == 0 )
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

      GeometryType elementType( GeometryType::simplex, dimgrid );
      for( int n = 0; n < dgf_.nofelements; ++n )
      {
        factory_.insertElement( elementType, dgf_.elements[ n ] );
        for( int face = 0; face <= dimgrid; ++face )
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
        GeometryType type( GeometryType::simplex, dimgrid-1 );
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


  template <int dimw>
  inline bool DGFGridFactory< ALUCubeGrid< 2, dimw > >
  ::generate( std::istream &file, MPICommunicatorType communicator, const std::string &filename )
  {
    const int dimgrid = 2;
    const int dimworld = dimw;
    dgf_.element = DuneGridFormatParser::Cube;
    dgf_.dimgrid = dimgrid;
    dgf_.dimw = dimworld;

    const bool isDGF = dgf_.isDuneGridFormat( file );
    file.seekg( 0 );
    if( !isDGF )
      return false;

    int rank = 0;
#if ALU2DGRID_PARALLEL
    MPI_Comm_rank( communicator, &rank );
#endif

    dgf::ALU2dGridParameterBlock parameter( file );

    if( rank == 0 )
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

      GeometryType elementType( GeometryType::cube, dimgrid );
      for( int n = 0; n < dgf_.nofelements; ++n )
      {
        factory_.insertElement( elementType, dgf_.elements[ n ] );
        for( int face = 0; face < 2*dimgrid; ++face )
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
        GeometryType type( GeometryType::simplex, dimgrid-1 );
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


  template <int dimw>
  inline bool DGFGridFactory< ALUConformGrid< 2, dimw > >
  ::generate( std::istream &file, MPICommunicatorType communicator, const std::string &filename )
  {
    const int dimgrid = 2;
    const int dimworld = dimw;
    dgf_.element = DuneGridFormatParser::Simplex;
    dgf_.dimgrid = dimgrid;
    dgf_.dimw = dimworld;

    const bool isDGF = dgf_.isDuneGridFormat( file );
    file.seekg( 0 );
    if( !isDGF )
      return false;

    int rank = 0;
#if ALU2DGRID_PARALLEL
    MPI_Comm_rank( communicator, &rank );
#endif

    dgf::ALU2dGridParameterBlock parameter( file );

    if( rank == 0 )
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

      GeometryType elementType( GeometryType::simplex, dimgrid );
      for( int n = 0; n < dgf_.nofelements; ++n )
      {
        factory_.insertElement( elementType, dgf_.elements[ n ] );
        for( int face = 0; face <= dimgrid; ++face )
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
        GeometryType type( GeometryType::cube, dimgrid-1 );
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

} // end namespace Dune
  /** \endcond */
