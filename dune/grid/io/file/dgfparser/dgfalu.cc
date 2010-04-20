// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**\cond */
namespace Dune
{

  inline void DGFGridFactory< ALUSimplexGrid< 3, 3 > >
  ::generate( const std::string &filename,
              MPICommunicatorType communicator )
  {
    const int dimworld = 3;
    dgf_.element = DuneGridFormatParser::Simplex;
    dgf_.dimgrid = dimworld;
    dgf_.dimw = dimworld;
    int rank = 0;
#if ALU3DGRID_PARALLEL
    MPI_Comm_rank( communicator, &rank );
#endif

    std::ifstream file( filename.c_str() );
    if (! file)
    {
      DUNE_THROW(DGFException,
                 "Macrofile " << filename << " not found");
    }

    if( dgf_.readDuneGrid( file, dimworld, dimworld ) )
    {
      if( dgf_.dimw != dimworld )
        DUNE_THROW( DGFException, "Macrofile " << filename << " is for "
                                               << "dimension " << dgf_.dimw
                                               << " and cannot be used to initialize an "
                                               << "ALUGrid of dimension " << dimension << ".");
      if (dimworld == 3)
        dgf_.setOrientation( 2, 3 );

      dgf::GridParameterBlock parameter( file );
      if( rank == 0 )
      {
        for( int n = 0; n < dgf_.nofvtx; ++n )
        {
          FieldVector< double, dimworld > pos;
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
            typedef MacroGrid::facemap_t::key_type Key;
            typedef MacroGrid::facemap_t::iterator Iterator;

            const Key key = ElementFaceUtil::generateFace( dimworld, dgf_.elements[ n ], face );
            const Iterator it = dgf_.facemap.find( key );
            if( it != dgf_.facemap.end() )
              factory_.insertBoundary( n, face, it->second );
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

      if ( ! parameter.dumpFileName().empty() )
        grid_ = factory_.createGrid( dgf_.facemap.empty(), false, parameter.dumpFileName() );
      else
        grid_ = factory_.createGrid( dgf_.facemap.empty(), true, filename );
    }
    else
    {
      grid_ =  callDirectly ("ALUSimplexGrid< 3 , 3 >", rank, filename.c_str(), communicator);
    }
  }
  inline void DGFGridFactory< ALUCubeGrid< 3, 3 > >
  ::generate( const std::string &filename,
              MPICommunicatorType communicator )
  {
    const int dimworld = 3;
    dgf_.element = DuneGridFormatParser::Cube;
    dgf_.dimgrid = dimworld;
    dgf_.dimw = dimworld;
    int rank = 0;
#if ALU3DGRID_PARALLEL
    MPI_Comm_rank( communicator, &rank );
#endif

    std::ifstream file( filename.c_str() );
    if (! file)
    {
      DUNE_THROW(DGFException,
                 "Macrofile " << filename << " not found");
    }

    if( dgf_.readDuneGrid( file, dimworld, dimworld ) )
    {
      if( dgf_.dimw != dimworld )
        DUNE_THROW( DGFException, "Macrofile " << filename << " is for "
                                               << "dimension " << dgf_.dimw
                                               << " and cannot be used to initialize an "
                                               << "ALUGrid of dimension " << dimension << ".");

      dgf::GridParameterBlock parameter( file );
      if( rank == 0 )
      {
        for( int n = 0; n < dgf_.nofvtx; ++n )
        {
          FieldVector< double, dimworld > pos;
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
            typedef MacroGrid::facemap_t::key_type Key;
            typedef MacroGrid::facemap_t::iterator Iterator;

            const Key key = ElementFaceUtil::generateFace( dimworld, dgf_.elements[ n ], face );
            const Iterator it = dgf_.facemap.find( key );
            if( it != dgf_.facemap.end() )
              factory_.insertBoundary( n, face, it->second );
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

      if ( ! parameter.dumpFileName().empty() )
        grid_ = factory_.createGrid( dgf_.facemap.empty(), false, parameter.dumpFileName() );
      else
        grid_ = factory_.createGrid( dgf_.facemap.empty(), true, filename );
    }
    else
    {
      grid_ =  callDirectly ("ALUCubeGrid< 3 , 3 >", rank, filename.c_str(), communicator);
    }
  }

  inline void DGFGridFactory< ALUSimplexGrid< 2, 2 > >
  ::generate( const std::string &filename,
              MPICommunicatorType communicator )
  {
    const int dimworld = 2;
    dgf_.element = DuneGridFormatParser::Simplex;
    dgf_.dimgrid = dimworld;
    dgf_.dimw = dimworld;
    int rank = 0;
#if ALU2DGRID_PARALLEL
    MPI_Comm_rank( communicator, &rank );
#endif

    std::ifstream file( filename.c_str() );
    if (! file)
    {
      DUNE_THROW(DGFException,
                 "Macrofile " << filename << " not found");
    }

    if( dgf_.readDuneGrid( file, dimworld, dimworld ) )
    {
      if( dgf_.dimw != dimworld )
        DUNE_THROW( DGFException, "Macrofile " << filename << " is for "
                                               << "dimension " << dgf_.dimw
                                               << " and cannot be used to initialize an "
                                               << "ALUGrid of dimension " << dimension << ".");

      dgf::GridParameterBlock parameter( file );

      if( rank == 0 )
      {
        for( int n = 0; n < dgf_.nofvtx; ++n )
        {
          FieldVector< double, dimworld > pos;
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
            typedef MacroGrid::facemap_t::key_type Key;
            typedef MacroGrid::facemap_t::iterator Iterator;

            const Key key = ElementFaceUtil::generateFace( dimworld, dgf_.elements[ n ], face );
            const Iterator it = dgf_.facemap.find( key );
            if( it != dgf_.facemap.end() )
              factory_.insertBoundary( n, face, it->second );
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

      if ( ! parameter.dumpFileName().empty() )
        grid_ = factory_.createGrid( dgf_.facemap.empty(), false, parameter.dumpFileName() );
      else
        grid_ = factory_.createGrid( dgf_.facemap.empty(), true, filename );
    }
    else
    {
      if( fileExists( filename.c_str() ) )
        grid_ = new Grid( filename );
      else
        DUNE_THROW( GridError, "Unable to create a 2d ALUGrid from '"
                    << filename << "'." );
    }
  }
  inline void DGFGridFactory< ALUConformGrid< 2, 2 > >
  ::generate( const std::string &filename,
              MPICommunicatorType communicator )
  {
    const int dimworld = 2;
    dgf_.element = DuneGridFormatParser::Simplex;
    dgf_.dimgrid = dimworld;
    dgf_.dimw = dimworld;
    int rank = 0;
#if ALU2DGRID_PARALLEL
    MPI_Comm_rank( communicator, &rank );
#endif

    std::ifstream file( filename.c_str() );
    if (! file)
    {
      DUNE_THROW(DGFException,
                 "Macrofile " << filename << " not found");
    }

    if( dgf_.readDuneGrid( file, dimworld, dimworld ) )
    {
      if( dgf_.dimw != dimworld )
        DUNE_THROW( DGFException, "Macrofile " << filename << " is for "
                                               << "dimension " << dgf_.dimw
                                               << " and cannot be used to initialize an "
                                               << "ALUGrid of dimension " << dimension << ".");

      dgf::GridParameterBlock parameter( file );

      if( rank == 0 )
      {
        for( int n = 0; n < dgf_.nofvtx; ++n )
        {
          FieldVector< double, dimworld > pos;
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
            typedef MacroGrid::facemap_t::key_type Key;
            typedef MacroGrid::facemap_t::iterator Iterator;

            const Key key = ElementFaceUtil::generateFace( dimworld, dgf_.elements[ n ], face );
            const Iterator it = dgf_.facemap.find( key );
            if( it != dgf_.facemap.end() )
              factory_.insertBoundary( n, face, it->second );
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

      if ( ! parameter.dumpFileName().empty() )
        grid_ = factory_.createGrid( dgf_.facemap.empty(), false, parameter.dumpFileName() );
      else
        grid_ = factory_.createGrid( dgf_.facemap.empty(), true, filename );
    }
    else
    {
      if( fileExists( filename.c_str() ) )
        grid_ = new Grid( filename );
      else
        DUNE_THROW( GridError, "Unable to create a 2d ALUGrid from '"
                    << filename << "'." );
    }
  }
} // end namespace Dune
  /** \endcond */
