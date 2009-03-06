// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// Of course, we would like to put this into the lib, but unfortunately HAVE_MPI
// is unpredictable.
//#include <config.h>

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
      DUNE_THROW( GridError,
                  "ALU3dGridFactory allows insertion only for rank 0." );
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
      DUNE_THROW( GridError,
                  "ALU3dGridFactory allows insertion only for rank 0." );
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
    {
      DUNE_THROW( GridError,
                  "ALU3dGridFactory allows insertion only for rank 0." );
    }
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
  ::insertBoundary ( const GeometryType &geometry,
                     const DGFEntityKey< unsigned int > &key,
                     const int id )
  {
#if ALU3DGRID_PARALLEL
    if( rank_ != 0 )
    {
      DUNE_THROW( GridError,
                  "ALU3dGridFactory allows insertion only for rank 0." );
    }
#endif
    assertGeometryType( geometry );
    if( geometry.dim() != dimension-1 )
    {
      DUNE_THROW( GridError, "Only 2-dimensional boundaries can be inserted "
                  "into a 3-dimensional ALUGrid." );
    }
    if( (unsigned int)key.size() != numFaceCorners )
      DUNE_THROW( GridError, "Wrong number of vertices." );

    std::pair< FaceType, int > boundaryId;
    for( unsigned int i = 0; i < numFaceCorners; ++i )
      boundaryId.first[ i ] = key.origKey( i );
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
      return new GridType( communicator_ );
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
    GridType *grid = new GridType( filename_, communicator_ );
#else
    GridType *grid = new GridType( filename_ );
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
  ::correctElementOrientation ()
  {
    // orientation is only corrected for tetrahedra
    if( elementType == hexa )
      return;

    const typename ElementVector::iterator elementEnd = elements_.end();
    for( typename ElementVector::iterator elementIt = elements_.begin();
         elementIt != elementEnd; ++elementIt )
    {
      ElementType &element = *elementIt;

      const VertexType &p1 = vertices_[ element[ 1 ] ];

      VertexType p2 = vertices_[ element[ 2 ] ]; p2 -= p1;
      VertexType p3 = vertices_[ element[ 3 ] ]; p3 -= p1;
      VertexType p0 = vertices_[ element[ 0 ] ]; p0 -= p1;

      VertexType n;
      n[ 0 ] = p2[ 1 ] * p3[ 2 ] - p3[ 1 ] * p2[ 2 ];
      n[ 1 ] = p2[ 2 ] * p3[ 0 ] - p3[ 2 ] * p2[ 0 ];
      n[ 2 ] = p2[ 0 ] * p3[ 1 ] - p3[ 0 ] * p2[ 1 ];

      if( n * p0 < 0 )
      {
        int tmp = element[ 2 ];
        element[ 2 ] = element[ 3 ];
        element[ 3 ] = tmp;
      }
    }
  }


  template< template< int, int > class ALUGrid >
  inline void ALU3dGridFactory< ALUGrid >
  ::recreateBoundaryIds ( const int defaultId )
  {
    typedef std::set< DGFEntityKey< unsigned int > > FaceSet;
    typedef typename FaceSet::iterator FaceIterator;
    FaceSet faceSet;

    const GeometryType faceGeo( elementType == tetra
                                ? GeometryType::simplex : GeometryType::cube,
                                dimension-1 );

    // Add all (potential) boundary faces to faceSet (see DGF parser)
    const typename ElementVector::iterator elementEnd = elements_.end();
    for( typename ElementVector::iterator elementIt = elements_.begin();
         elementIt != elementEnd; ++elementIt )
    {
      const int numFaces = ElementFaceUtil::nofFaces( dimension, *elementIt );
      for( int face = 0; face < numFaces; ++face )
      {
        DGFEntityKey< unsigned int > key
          = ElementFaceUtil::generateFace( dimension, *elementIt, face );

        const FaceIterator pos = faceSet.find( key );
        if( pos != faceSet.end() )
          faceSet.erase( key );
        else
          faceSet.insert( key );
      }
    }

    // swap current boundary ids with an empty vector
    BoundaryIdVector boundaryIds;
    boundaryIds_.swap( boundaryIds );
    assert( boundaryIds_.size() == 0 );

    // add all current boundary ids again (with their reordered keys)
    const typename BoundaryIdVector::iterator bndEnd = boundaryIds.end();
    for( typename BoundaryIdVector::iterator bndIt = boundaryIds.begin();
         bndIt != bndEnd; ++bndIt )
    {
      std::vector< unsigned int > face( numFaceCorners );
      for( unsigned int i = 0; i < numFaceCorners; ++i )
        face[ i ] = bndIt->first[ i ];
      DGFEntityKey< unsigned int > key( face );

      const FaceIterator pos = faceSet.find( key );
      if( pos == faceSet.end() )
        continue;

      insertBoundary( faceGeo, *pos, bndIt->second );
      faceSet.erase( pos );
    }

    // add all new boundaries (with defaultId)
    const FaceIterator faceEnd = faceSet.end();
    for( FaceIterator faceIt = faceSet.begin(); faceIt != faceEnd; ++faceIt )
      insertBoundary( faceGeo, *faceIt, defaultId );
  }

}
