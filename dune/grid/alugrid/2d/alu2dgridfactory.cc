// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// Of course, we would like to put this into the lib, but unfortunately HAVE_MPI
// is unpredictable.
//#include <config.h>

#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>

#include <dune/grid/alugrid/2d/alu2dgridfactory.hh>

namespace Dune
{
  template< template< int, int > class ALUGrid >
  ALU2dGridFactory<ALUGrid>
  :: ALU2dGridFactory ( bool removeGeneratedFile )
    : globalProjection_ ( 0 ),
      numFacesInserted_ ( 0 )
  {}

  template< template< int, int > class ALUGrid >
  ALU2dGridFactory<ALUGrid>
  :: ALU2dGridFactory ( const std::string &filename )
    : globalProjection_ ( 0 ),
      numFacesInserted_ ( 0 )
  {}


  template< template< int, int > class ALUGrid >
  ALU2dGridFactory<ALUGrid> :: ~ALU2dGridFactory ()
  {}


  template< template< int, int > class ALUGrid >
  void ALU2dGridFactory<ALUGrid> :: insertVertex ( const VertexType &pos )
  {
    vertices_.push_back( pos );
  }


  template< template< int, int > class ALUGrid >
  void ALU2dGridFactory<ALUGrid>
  :: insertElement ( const GeometryType &geometry,
                     const std::vector< unsigned int > &vertices )
  {
    if( !geometry.isTriangle() )
      DUNE_THROW( GridError, "Only triangles can be inserted into "
                  "ALUGrid< 2, 2 >." );

    if( vertices.size() != numCorners )
      DUNE_THROW( GridError, "Wrong number of vertices." );

    elements_.push_back( vertices );
  }


  template< template< int, int > class ALUGrid >
  void ALU2dGridFactory<ALUGrid>
  :: insertBoundary ( const GeometryType &geometry,
                      const std::vector< unsigned int > &vertices,
                      const int id )
  {
    // lines can be either cube or simplex
    assert( geometry.isSimplex() || geometry.isCube() );
    if( geometry.dim() != dimension-1 )
    {
      DUNE_THROW( GridError, "Only 2-dimensional boundaries can be inserted "
                  "into a 2-dimensional ALUGrid." );
    }
    if( vertices.size() != numFaceCorners )
      DUNE_THROW( GridError, "Wrong number of vertices." );

    std::pair< FaceType, int > boundaryId;
    for( unsigned int i = 0; i < numFaceCorners; ++i )
      boundaryId.first[ i ] = vertices[ i ];
    boundaryId.second = id;
    boundaryIds_.push_back( boundaryId );
    ++numFacesInserted_;
  }


  template< template< int, int > class ALUGrid >
  void ALU2dGridFactory<ALUGrid>
  ::insertBoundary ( const int element, const int face, const int id )
  {
    if( (element < 0) || (element >= (int)elements_.size()) )
      DUNE_THROW( RangeError, "ALU2dGridFactory::insertBoundary: invalid element index given." );

    std::pair< FaceType, int > boundaryId;
    generateFace( elements_[ element ], face, boundaryId.first );
    boundaryId.second = id;
    boundaryIds_.push_back( boundaryId );
    ++numFacesInserted_;
  }

  template< template< int, int > class ALUGrid >
  void ALU2dGridFactory< ALUGrid > ::
  insertBoundaryProjection( const DuneBoundaryProjectionType& bndProjection )
  {
    if( globalProjection_ )
      DUNE_THROW(InvalidStateException,"You can only insert one globalProjection");

    globalProjection_ = &bndProjection;
  }

  template< template< int, int > class ALUGrid >
  void ALU2dGridFactory< ALUGrid > ::
  insertBoundaryProjection ( const GeometryType &type,
                             const std::vector< unsigned int > &vertices,
                             const DuneBoundaryProjectionType *projection )
  {
    assert( type.isSimplex() || type.isCube() );
    if( (int)type.dim() != dimension-1 )
      DUNE_THROW( GridError, "Inserting boundary face of wrong dimension: " << type.dim() );

    FaceType faceId;
    copyAndSort( vertices, faceId );

    if( vertices.size() != numFaceCorners )
      DUNE_THROW( GridError, "Wrong number of face vertices passed: " << vertices.size() << "." );

    if( boundaryProjections_.find( faceId ) != boundaryProjections_.end() )
      DUNE_THROW( GridError, "Only one boundary projection can be attached to a face." );
    boundaryProjections_[ faceId ] = projection;
  }

  template< template< int, int > class ALUGrid >
  void ALU2dGridFactory< ALUGrid > ::
  insertBoundarySegment ( const std::vector< unsigned int >& vertices )
  {
    DUNE_THROW( NotImplemented, "insertBoundarySegment with a single argument" );
  }

  template< template< int, int > class ALUGrid >
  void ALU2dGridFactory< ALUGrid > ::
  insertBoundarySegment ( const std::vector< unsigned int >& vertices,
                          const shared_ptr<BoundarySegment<2,2> >& boundarySegment )
  {
    FaceType faceId;
    copyAndSort( vertices, faceId );

    if( vertices.size() != numFaceCorners )
      DUNE_THROW( GridError, "Wrong number of face vertices passed: " << vertices.size() << "." );

    if( boundaryProjections_.find( faceId ) != boundaryProjections_.end() )
      DUNE_THROW( GridError, "Only one boundary projection can be attached to a face." );

    const size_t numVx = vertices.size();
    GeometryType type( (numVx == 3) ? GeometryType::simplex : GeometryType::cube, dimension-1 );

    std::vector< VertexType > coords( numVx );
    for( size_t i = 0; i < numVx; ++i )
    {
      // get global coordinate and copy it
      const VertexType &x = vertices_[ vertices[ i ] ];
      for( unsigned int j = 0; j < dimensionworld; ++j )
        coords[ i ][ j ] = x[ j ];
    }
    BoundarySegmentWrapperType* prj
      = new BoundarySegmentWrapperType( type, coords, boundarySegment );
    boundaryProjections_[ faceId ] = prj;
#ifndef NDEBUG
    // consistency check
    for( size_t i = 0; i < numVx; ++i )
    {
      VertexType global = (*prj)( coords [ i ] );
      if( (global - coords[ i ]).two_norm() > 1e-6 )
      {
        DUNE_THROW(GridError,"BoundarySegment does not map face vertices to face vertices.");
      }
    }
#endif
  }

  template< template< int, int > class ALUGrid >
  ALUGrid< 2, 2 > *ALU2dGridFactory< ALUGrid >::createGrid ()
  {
    return createGrid( true, true, "" );
  }

  template< template< int, int > class ALUGrid >
  ALUGrid< 2, 2 > *ALU2dGridFactory< ALUGrid >
  ::createGrid ( const bool addMissingBoundaries, const std::string dgfName )
  {
    return createGrid( addMissingBoundaries, true, dgfName );
  }

  template< template< int, int > class ALUGrid >
  ALUGrid< 2, 2 > *ALU2dGridFactory< ALUGrid >
  ::createGrid ( const bool addMissingBoundaries, bool temporary, const std::string name )
  {
    if( addMissingBoundaries )
      recreateBoundaryIds();

#ifndef ALUGRID_NOTEMPFILE_2D
    std::string filename ( temporary ?
                           temporaryFileName( name ) :
                           name );
#else
    std::string filename ( name );
#endif

    std::ofstream outfile;
    std::stringstream temp;

    std::ostream* outptr = 0;
#ifdef ALUGRID_NOTEMPFILE_2D
    if( temporary )
      outptr = & temp;
    else
#endif
    {
      outfile.open( filename.c_str() , std::ios::out );
      outptr = & outfile;
    }
    std::ostream& out = *outptr ;

    out.setf( std::ios_base::scientific, std::ios_base::floatfield );
    out.precision( 16 );
    out << "!Triangles" << std::endl;

    const unsigned int numVertices = vertices_.size();
    // now start writing grid
    out << numVertices << std :: endl;
    typedef typename std :: vector< VertexType > :: iterator VertexIteratorType;
    const VertexIteratorType endV = vertices_.end();
    for( VertexIteratorType it = vertices_.begin(); it != endV; ++it )
    {
      const VertexType &vertex = *it;
      out << vertex[ 0 ];
      for( unsigned int i = 1; i < dimensionworld; ++i )
        out << " " << vertex[ i ];
      out << std :: endl;
    }

    out << elements_.size() << std :: endl;
    typedef typename ElementVector::iterator ElementIteratorType;
    const ElementIteratorType endE = elements_.end();
    for( ElementIteratorType it = elements_.begin(); it != endE; ++it )
    {
      array< unsigned int, numCorners > element;
      for( unsigned int i = 0; i < numCorners; ++i )
      {
        element[ i ] = (*it)[ i ];
      }

      out << element[ 0 ];
      for( unsigned int i = 1; i < numCorners; ++i )
        out << "  " << element[ i ];
      out << std :: endl;
    }

    const size_t boundarySegments = boundaryIds_.size();
    out << boundarySegments << std :: endl;
    typedef typename std :: vector< std :: pair< FaceType, int > > :: iterator
    BoundaryIdIteratorType;

    BoundaryProjectionVector* bndProjections = 0;
    const size_t bndProjectionSize = boundaryProjections_.size();
    if( bndProjectionSize > 0 )
    {
      // the memory is freed by the grid on destruction
      bndProjections = new BoundaryProjectionVector( boundarySegments,
                                                     (DuneBoundaryProjectionType*) 0 );
    }

    const BoundaryIdIteratorType endB = boundaryIds_.end();
    int segmentIndex = 0;
    for( BoundaryIdIteratorType it = boundaryIds_.begin(); it != endB; ++it, ++segmentIndex )
    {
      const std :: pair< FaceType, int > &boundaryId = *it;
      out << "-" << boundaryId.second ;
      for( unsigned int i = 0; i < numFaceCorners; ++i )
        out << "  " << boundaryId.first[ i ];
      out << std :: endl;

      if( bndProjectionSize > 0 )
      {
        // generate boundary segment pointer
        FaceType faceId (boundaryId.first);
        std::sort( faceId.begin(), faceId.end() );

        const DuneBoundaryProjectionType* projection = boundaryProjections_[ faceId ];

        // if no projection given we use global projection, otherwise identity
        if( ! projection && globalProjection_ )
        {
          typedef BoundaryProjectionWrapper< dimensionworld > ProjectionWrapperType;
          // we need to wrap the global projection because of
          // deleting in desctructor of ALUGrid
          projection = new ProjectionWrapperType( *globalProjection_ );
        }

        // copy pointer
        (*bndProjections)[ segmentIndex ] = projection;
      }
    }

    outfile.close();

    vertices_.clear();
    elements_.clear();
    boundaryIds_.clear();
    boundaryProjections_.clear();

#ifdef ALUGRID_NOTEMPFILE_2D
    std::istream& inFile = temp;
#endif

    // if we have a vector reset global projection
    // because empty positions are filled with global projection anyway
    if( bndProjections ) globalProjection_ = 0;

    // ALUGrid is taking ownership of the bndProjections pointer
    Grid *grid =
#ifdef ALUGRID_NOTEMPFILE_2D
      ( temporary ) ? new Grid( filename, inFile, globalProjection_ , bndProjections ) :
#endif
      new Grid( filename, globalProjection_ , bndProjections );

#ifndef ALUGRID_NOTEMPFILE_2D
    if( temporary )
      std::remove( filename.c_str() );
#endif

    // remove pointer
    globalProjection_ = 0;

    return grid;
  }


  template< template< int, int > class ALUGrid >
  inline std::string
  ALU2dGridFactory<ALUGrid>::temporaryFileName ( const std::string& dgfName )
  {
    std::string filename ( dgfName );
    if( filename.empty() )
    {
      filename = "ALU2dGrid.XXXXXX";
    }
    else
    {
      filename += ".XXXXXX";
    }
    char filetemp[ FILENAME_MAX ];
    std :: strcpy( filetemp, filename.c_str() );
    const int fd = mkstemp( filetemp );
    if( fd < 0 )
      DUNE_THROW( IOError, "Unable to create temporary file." );
    close( fd );
    return std :: string( filetemp );
  }


  template< template< int, int > class ALUGrid >
  inline void ALU2dGridFactory<ALUGrid>
  ::generateFace ( const ElementType &element, const int f, FaceType &face )
  {
    if (f==0) { face[0]=element[0]; face[1]=element[1]; }
    if (f==1) { face[0]=element[0]; face[1]=element[2]; }
    if (f==2) { face[0]=element[1]; face[1]=element[2]; }
  }

  template< template< int, int > class ALUGrid >
  inline void ALU2dGridFactory<ALUGrid>
  ::recreateBoundaryIds ( const int defaultId )
  {
    typedef std::pair< unsigned int, int > SubEntity;
    typedef std::map< FaceType, SubEntity, FaceLess > FaceMap;
    typedef typename FaceMap::iterator FaceIterator;
    FaceMap faceMap;

    const GeometryType faceGeo( dimension-1 );

    const unsigned int numElements = elements_.size();
    for( unsigned int n = 0; n < numElements; ++n )
    {
      for( unsigned int face = 0; face < numFaces; ++face )
      {
        FaceType key;
        generateFace( elements_[ n ], face, key );
        std::sort( key.begin(), key.end() );

        const FaceIterator pos = faceMap.find( key );
        if( pos != faceMap.end() )
          faceMap.erase( key );
        else
          faceMap.insert( std::make_pair( key, SubEntity( n, face ) ) );
      }
    }

    // swap current boundary ids with an empty vector
    BoundaryIdVector boundaryIds;
    boundaryIds_.swap( boundaryIds );
    assert( boundaryIds_.size() == 0 );

    // add all current boundary ids again (with their reordered keys)
    typedef typename BoundaryIdVector::iterator BoundaryIterator;
    const BoundaryIterator bndEnd = boundaryIds.end();
    for( BoundaryIterator bndIt = boundaryIds.begin(); bndIt != bndEnd; ++bndIt )
    {
      FaceType key = bndIt->first;
      std::sort( key.begin(), key.end() );
      const FaceIterator pos = faceMap.find( key );
      if( pos == faceMap.end() )
        continue;
      insertBoundary( pos->second.first, pos->second.second, bndIt->second );
      faceMap.erase( pos );
    }

    // add all new boundaries (with defaultId)
    const FaceIterator faceEnd = faceMap.end();
    for( FaceIterator faceIt = faceMap.begin(); faceIt != faceEnd; ++faceIt )
      insertBoundary( faceIt->second.first, faceIt->second.second, defaultId );
  }

}
