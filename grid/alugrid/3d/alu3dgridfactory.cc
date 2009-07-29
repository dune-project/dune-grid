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

#include <dune/grid/alugrid/3d/alu3dgridfactory.hh>

namespace Dune
{

  template< template< int, int > class ALUGrid >
  ALU3dGridFactory< ALUGrid >
  :: ALU3dGridFactory ( const MPICommunicatorType &communicator,
                        bool removeGeneratedFile )
    : filename_( temporaryFileName() ),
      communicator_( communicator ),
      removeGeneratedFile_( removeGeneratedFile )
  {
#if ALU3DGRID_PARALLEL
    MPI_Comm_rank( communicator, &rank_ );
#endif
  }


  template< template< int, int > class ALUGrid >
  ALU3dGridFactory< ALUGrid >
  :: ALU3dGridFactory ( const std::string &filename,
                        const MPICommunicatorType &communicator )
    : filename_( filename + ".ALU3dGrid" ),
      communicator_( communicator ),
      removeGeneratedFile_( false )
  {
#if ALU3DGRID_PARALLEL
    MPI_Comm_rank( communicator, &rank_ );
#endif
  }


  template< template< int, int > class ALUGrid >
  ALU3dGridFactory< ALUGrid > :: ~ALU3dGridFactory ()
  {}


  template< template< int, int > class ALUGrid >
  void ALU3dGridFactory< ALUGrid > :: insertVertex ( const VertexType &pos )
  {
#if ALU3DGRID_PARALLEL
    if( rank_ != 0 )
      DUNE_THROW( GridError, "ALU3dGridFactory allows insertion only for rank 0." );
#endif
    vertices_.push_back( pos );
  }


  template< template< int, int > class ALUGrid >
  void ALU3dGridFactory< ALUGrid >
  :: insertElement ( const GeometryType &geometry,
                     const std::vector< unsigned int > &vertices )
  {
#if ALU3DGRID_PARALLEL
    if( rank_ != 0 )
      DUNE_THROW( GridError, "ALU3dGridFactory allows insertion only for rank 0." );
#endif
    assertGeometryType( geometry );
    if( geometry.dim() != dimension )
      DUNE_THROW( GridError, "Only 3-dimensional elements can be inserted "
                  "into a 3-dimensional ALUGrid." );
    if( vertices.size() != numCorners )
      DUNE_THROW( GridError, "Wrong number of vertices." );

    elements_.push_back( vertices );
  }


  template< template< int, int > class ALUGrid >
  void ALU3dGridFactory< ALUGrid >
  :: insertBoundary ( const GeometryType &geometry,
                      const std::vector< unsigned int > &vertices,
                      const int id )
  {
#if ALU3DGRID_PARALLEL
    if( rank_ != 0 )
      DUNE_THROW( GridError, "ALU3dGridFactory allows insertion only for rank 0." );
#endif
    assertGeometryType( geometry );
    if( geometry.dim() != dimension-1 )
    {
      DUNE_THROW( GridError, "Only 2-dimensional boundaries can be inserted "
                  "into a 3-dimensional ALUGrid." );
    }
    if( vertices.size() != numFaceCorners )
      DUNE_THROW( GridError, "Wrong number of vertices." );

    std::pair< FaceType, int > boundaryId;
    for( unsigned int i = 0; i < numFaceCorners; ++i )
    {
      const unsigned int j = FaceTopologyMappingType::dune2aluVertex( i );
      boundaryId.first[ j ] = vertices[ i ];
    }
    boundaryId.second = id;
    boundaryIds_.push_back( boundaryId );
  }


  template< template< int, int > class ALUGrid >
  void ALU3dGridFactory< ALUGrid >
  ::insertBoundary ( const int element, const int face, const int id )
  {
#if ALU3DGRID_PARALLEL
    if( rank_ != 0 )
      DUNE_THROW( GridError, "ALU3dGridFactory allows insertion only for rank 0." );
#endif

    if( (element < 0) || (element >= (int)elements_.size()) )
      DUNE_THROW( RangeError, "ALU3dGridFactory::insertBoundary: invalid element index given." );

    std::pair< FaceType, int > boundaryId;
    generateFace( elements_[ element ], face, boundaryId.first );
    boundaryId.second = id;
    boundaryIds_.push_back( boundaryId );
  }


  template< template< int, int > class ALUGrid >
  ALUGrid< 3, 3 > *ALU3dGridFactory< ALUGrid >::createGrid ()
  {
    return createGrid( true );
  }


  template< template< int, int > class ALUGrid >
  ALUGrid< 3, 3 > *ALU3dGridFactory< ALUGrid >
  ::createGrid ( const bool addMissingBoundaries )
  {
#if ALU3DGRID_PARALLEL
    if( rank_ != 0 )
      return new Grid( communicator_ );
#endif

    correctElementOrientation();
    if( addMissingBoundaries )
      recreateBoundaryIds();

    std :: ofstream out( filename_.c_str() );
    if( elementType == tetra )
      out << "!Tetrahedra";
    else
      out << "!Hexahedra";

    const unsigned int numVertices = vertices_.size();
    // print information about vertices and elements
    // to header to have an easy check
    out << "  ( noVertices = " << numVertices;
    out << " | noElements = " << elements_.size() << " )" << std :: endl;

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
        const unsigned int j = ElementTopologyMappingType::dune2aluVertex( i );
        element[ j ] = (*it)[ i ];
      }

      out << element[ 0 ];
      for( unsigned int i = 1; i < numCorners; ++i )
        out << "  " << element[ i ];
      out << std :: endl;
    }

    out << boundaryIds_.size() << std :: endl;
    typedef typename std :: vector< std :: pair< FaceType, int > > :: iterator
    BoundaryIdIteratorType;
    const BoundaryIdIteratorType endB = boundaryIds_.end();
    for( BoundaryIdIteratorType it = boundaryIds_.begin(); it != endB; ++it )
    {
      const std :: pair< FaceType, int > &boundaryId = *it;
      out << "-" << boundaryId.second << "  " << numFaceCorners;
      for( unsigned int i = 0; i < numFaceCorners; ++i )
        out << "  " << boundaryId.first[ i ];
      out << std :: endl;
    }

    for( unsigned int i = 0; i < numVertices; ++i )
      out << i << "  -1" << std :: endl;
    out.close();

    vertices_.clear();
    elements_.clear();
    boundaryIds_.clear();

#if ALU3DGRID_PARALLEL
    Grid *grid = new Grid( filename_, communicator_ );
#else
    Grid *grid = new Grid( filename_ );
#endif
    if( removeGeneratedFile_ )
      remove( filename_.c_str() );
    return grid;
  }


  template< template< int, int > class ALUGrid >
  inline std::string
  ALU3dGridFactory< ALUGrid >::temporaryFileName ()
  {
    char filetemp[ FILENAME_MAX ];
    std :: strcpy( filetemp, "ALU3dGrid.XXXXXX" );
    const int fd = mkstemp( filetemp );
    if( fd < 0 )
      DUNE_THROW( IOError, "Unable to create temporary file." );
    close( fd );
    return std :: string( filetemp );
  }


  template< template< int, int > class ALUGrid >
  inline void ALU3dGridFactory< ALUGrid >
  ::generateFace ( const ElementType &element, const int f, FaceType &face )
  {
    const int falu = ElementTopologyMapping< elementType >::generic2aluFace( f );
    for( unsigned int i = 0; i < numFaceCorners; ++i )
    {
      const int j = ElementTopologyMapping< elementType >::faceVertex( falu, i );
      const int k = ElementTopologyMapping< elementType >::alu2genericVertex( j );
      face[ i ] = element[ k ];
    }
  }


  template< template< int, int > class ALUGrid >
  inline void
  ALU3dGridFactory< ALUGrid >::correctElementOrientation ()
  {
    const typename ElementVector::iterator elementEnd = elements_.end();
    for( typename ElementVector::iterator elementIt = elements_.begin();
         elementIt != elementEnd; ++elementIt )
    {
      ElementType &element = *elementIt;

      const VertexType &p0 = vertices_[ element[ 0 ] ];
      VertexType p1, p2, p3;

      if( elementType == tetra )
      {
        p1 = vertices_[ element[ 1 ] ];
        p2 = vertices_[ element[ 2 ] ];
        p3 = vertices_[ element[ 3 ] ];
      }
      else
      {
        p1 = vertices_[ element[ 1 ] ];
        p2 = vertices_[ element[ 2 ] ];
        p3 = vertices_[ element[ 4 ] ];
      }
      p1 -= p0; p2 -= p0; p3 -= p0;

      VertexType n;
      n[ 0 ] = p1[ 1 ] * p2[ 2 ] - p2[ 1 ] * p1[ 2 ];
      n[ 1 ] = p1[ 2 ] * p2[ 0 ] - p2[ 2 ] * p1[ 0 ];
      n[ 2 ] = p1[ 0 ] * p2[ 1 ] - p2[ 0 ] * p1[ 1 ];

      if( n * p3 > 0 )
        continue;

      if( elementType == hexa )
      {
        for( int i = 0; i < 4; ++i )
          exchange( element[ i ], element[ i+4 ] );
      }
      else
        exchange( element[ 2 ], element[ 3 ] );
    } // end of loop over all elements
  }


  template< template< int, int > class ALUGrid >
  inline void ALU3dGridFactory< ALUGrid >
  ::recreateBoundaryIds ( const int defaultId )
  {
    typedef std::pair< unsigned int, int > SubEntity;
    typedef std::map< FaceType, SubEntity, FaceLess > FaceMap;
    typedef typename FaceMap::iterator FaceIterator;
    FaceMap faceMap;

    const GeometryType faceGeo( elementType == tetra
                                ? GeometryType::simplex : GeometryType::cube,
                                dimension-1 );

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
