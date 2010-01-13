// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**\cond */
namespace Dune
{

  template<>
  inline ALUSimplexGrid< 3, 3 > *MacroGrid :: Impl< ALUSimplexGrid< 3, 3 > >
  :: generate ( MacroGrid &macroGrid,
                const char *filename,
                MPICommunicatorType communicator )
  {
    macroGrid.element = Simplex;

    int rank = 0;
#if ALU3DGRID_PARALLEL
    MPI_Comm_rank( communicator, &rank );
#endif

    std::ifstream file( filename );
    if( macroGrid.readDuneGrid( file, 3, 3 ) )
    {
      if( macroGrid.dimw != 3 )
        DUNE_THROW( DGFException, "Macrofile " << filename << " is for "
                                               << "dimension " << macroGrid.dimw
                                               << " and cannot be used to initialize an "
                                               << "ALUGrid of dimension 3." );
      macroGrid.setOrientation( 2, 3 );

      dgf::GridParameterBlock parameter( file );
      GridFactory< ALUSimplexGrid< 3, 3 > > factory( parameter.dumpFileName(), communicator );
      if( rank == 0 )
      {
        for( int n = 0; n < macroGrid.nofvtx; ++n )
        {
          FieldVector< double, 3 > pos;
          for( unsigned int i = 0; i < 3; ++i )
            pos[ i ] = macroGrid.vtx[ n ][ i ];
          factory.insertVertex( pos );
        }

        GeometryType elementType( GeometryType::simplex, 3 );
        for( int n = 0; n < macroGrid.nofelements; ++n )
        {
          factory.insertElement( elementType, macroGrid.elements[ n ] );
          for( int face = 0; face <= 3; ++face )
          {
            typedef MacroGrid::facemap_t::key_type Key;
            typedef MacroGrid::facemap_t::iterator Iterator;

            const Key key = ElementFaceUtil::generateFace( 3, macroGrid.elements[ n ], face );
            const Iterator it = macroGrid.facemap.find( key );
            if( it != macroGrid.facemap.end() )
              factory.insertBoundary( n, face, it->second );
          }
        }

        const int dimworld = 3;
        dgf::ProjectionBlock projectionBlock( file, dimworld );
        const DuneBoundaryProjection< dimworld > *projection
          = projectionBlock.defaultProjection< dimworld >();
        if( projection != 0 )
          factory.insertBoundaryProjection( *projection );
        const size_t numBoundaryProjections = projectionBlock.numBoundaryProjections();
        for( size_t i = 0; i < numBoundaryProjections; ++i )
        {
          GeometryType type( GeometryType::simplex, dimworld-1 );
          const std::vector< unsigned int > &vertices = projectionBlock.boundaryFace( i );
          const DuneBoundaryProjection< dimworld > *projection
            = projectionBlock.boundaryProjection< dimworld >( i );
          factory.insertBoundaryProjection( type, vertices, projection );
        }
      }

      return factory.createGrid( macroGrid.facemap.empty(), filename );
    }
    else
    {
      return MacroGrid :: Impl< ALUCubeGrid< 3, 3 > > ::
             callDirectly < ALUSimplexGrid< 3, 3 > >
               ("ALUSimplexGrid< 3 , 3 >", rank, filename, communicator);
    }
  }



  template<>
  inline ALUCubeGrid< 3, 3 > *MacroGrid :: Impl< ALUCubeGrid< 3, 3 > >
  :: generate ( MacroGrid &macroGrid,
                const char *filename,
                MPICommunicatorType communicator )
  {
    macroGrid.element = Cube;

    int rank = 0;
#if ALU3DGRID_PARALLEL
    MPI_Comm_rank( communicator, &rank );
#endif

    std :: ifstream file( filename );
    if( macroGrid.readDuneGrid( file, 3, 3 ) )
    {
      if( macroGrid.dimw != 3 )
        DUNE_THROW( DGFException, "Macrofile " << filename << " is for "
                                               << "dimension " << macroGrid.dimw
                                               << " and cannot be used to initialize an "
                                               << "ALUGrid of dimension 3." );

      dgf::GridParameterBlock parameter( file );
      GridFactory< ALUCubeGrid< 3, 3 > > factory( parameter.dumpFileName(), communicator );
      if( rank == 0 )
      {
        for( int n = 0; n < macroGrid.nofvtx; ++n )
        {
          FieldVector< double, 3 > pos;
          for( unsigned int i = 0; i < 3; ++i )
            pos[ i ] = macroGrid.vtx[ n ][ i ];
          factory.insertVertex( pos );
        }

        GeometryType elementType( GeometryType::cube, 3 );
        for( int n = 0; n < macroGrid.nofelements; ++n )
        {
          factory.insertElement( elementType, macroGrid.elements[ n ] );
          for( int face = 0; face < 6; ++face )
          {
            typedef MacroGrid::facemap_t::key_type Key;
            typedef MacroGrid::facemap_t::iterator Iterator;

            const Key key = ElementFaceUtil::generateFace( 3, macroGrid.elements[ n ], face );
            const Iterator it = macroGrid.facemap.find( key );
            if( it != macroGrid.facemap.end() )
              factory.insertBoundary( n, face, it->second );
          }
        }

        const int dimworld = 3;
        dgf::ProjectionBlock projectionBlock( file, dimworld );
        const DuneBoundaryProjection< dimworld > *projection
          = projectionBlock.defaultProjection< dimworld >();
        if( projection != 0 )
          factory.insertBoundaryProjection( *projection );
        const size_t numBoundaryProjections = projectionBlock.numBoundaryProjections();
        for( size_t i = 0; i < numBoundaryProjections; ++i )
        {
          GeometryType type( GeometryType::cube, dimworld-1 );
          const std::vector< unsigned int > &vertices = projectionBlock.boundaryFace( i );
          const DuneBoundaryProjection< dimworld > *projection
            = projectionBlock.boundaryProjection< dimworld >( i );
          factory.insertBoundaryProjection( type, vertices, projection );
        }
      }

      return factory.createGrid( macroGrid.facemap.empty(), filename );
    }
    else
    {
      return callDirectly< ALUCubeGrid< 3, 3 > >
               ("ALUCubeGrid< 3 , 3 >", rank, filename, communicator);
    }
  }

  template <>
  template <class GridImp>
  inline GridImp*  MacroGrid :: Impl< ALUCubeGrid< 3, 3 > >
  :: callDirectly( const char* gridname,
                   const int rank,
                   const char* filename,
                   MPICommunicatorType communicator )
  {
#if ALU3DGRID_PARALLEL
    std :: stringstream tmps;
    tmps << filename << "." << rank;
    const std :: string &tmp = tmps.str();

    if( fileExists( tmp.c_str() ) )
      return new GridImp( tmp.c_str(), communicator );
#endif
    if( fileExists( filename ) )
    {
      if( rank == 0 )
        return new GridImp( filename );
      else
        return new GridImp( );
    }
    DUNE_THROW( GridError, "Unable to create " << gridname << " from '"
                                               << filename << "'." );
  }

  template <>
  template <class GridImp>
  inline GridImp* MacroGrid :: Impl< ALUConformGrid< 2, 2 > >
  :: generate ( const char* gridname,
                MacroGrid &macroGrid,
                const char* filename,
                MPICommunicatorType communicator )
  {
    macroGrid.element = Simplex;

    int rank = 0;
#if ALU2DGRID_PARALLEL
    MPI_Comm_rank( communicator, &rank );
#endif

    std :: ifstream file( filename );
    if( macroGrid.readDuneGrid( file, 2, 2 ) )
    {
      if( macroGrid.dimw != 2 )
        DUNE_THROW( DGFException, "Macrofile " << filename << " is for "
                                               << "dimension " << macroGrid.dimw
                                               << " and cannot be used to initialize an "
                                               << "ALUGrid of dimension 2." );

      dgf::GridParameterBlock parameter( file );

      // for parallel runs only rank 0 should create dumpFile
      // all other processors create temporary file
      typedef GridFactory< GridImp > GridFactoryType;
      std::auto_ptr< GridFactoryType > factPtr
        ( ( rank == 0 ) ?
        new GridFactoryType( parameter.dumpFileName() ) :
        new GridFactoryType( )
        );

      // get reference
      GridFactoryType& factory = * factPtr;

      for( int n = 0; n < macroGrid.nofvtx; ++n )
      {
        FieldVector< double, 2 > pos;
        for( unsigned int i = 0; i < 2; ++i )
          pos[ i ] = macroGrid.vtx[ n ][ i ];
        factory.insertVertex( pos );
      }

      GeometryType elementType( GeometryType::simplex, 2 );
      for( int n = 0; n < macroGrid.nofelements; ++n )
      {
        factory.insertElement( elementType, macroGrid.elements[ n ] );
        for( int face = 0; face < 3; ++face )
        {
          typedef MacroGrid::facemap_t::key_type Key;
          typedef MacroGrid::facemap_t::iterator Iterator;

          const Key key = ElementFaceUtil::generateFace( 2, macroGrid.elements[ n ], face );
          const Iterator it = macroGrid.facemap.find( key );
          if( it != macroGrid.facemap.end() )
            factory.insertBoundary( n, face, it->second );
        }
      }
      return factory.createGrid( macroGrid.facemap.empty(), filename );
    }
    else
    {
      if( fileExists( filename ) )
        return new GridImp( filename );
    }
    DUNE_THROW( GridError, "Unable to create " << gridname << " from '"
                                               << filename << "'." );
  }

} // end namespace Dune
  /** \endcond */
