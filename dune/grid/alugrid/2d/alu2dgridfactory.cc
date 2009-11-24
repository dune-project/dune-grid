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
    : filename_( temporaryFileName() ),
      removeGeneratedFile_( removeGeneratedFile )
  {}


  template< template< int, int > class ALUGrid >
  ALU2dGridFactory<ALUGrid>
  :: ALU2dGridFactory ( const std::string &filename )
    : filename_( filename.empty() ? temporaryFileName() : filename ),
      removeGeneratedFile_( filename.empty() )
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
    assertGeometryType( geometry );
    if( geometry.dim() != dimension )
      DUNE_THROW( GridError, "Only 2-dimensional elements can be inserted "
                  "into a 2-dimensional ALUGrid." );
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
    assertGeometryType( geometry );
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
  }


  template< template< int, int > class ALUGrid >
  ALUGrid< 2, 2 > *ALU2dGridFactory<ALUGrid>::createGrid ()
  {
    return createGrid( true );
  }


  template< template< int, int > class ALUGrid >
  ALUGrid< 2, 2 > *ALU2dGridFactory<ALUGrid>::createGrid ( const bool addMissingBoundaries )
  {
    if( addMissingBoundaries )
      recreateBoundaryIds();

    std::ofstream out( filename_.c_str() );
    out.setf( std::ios_base::scientific, std::ios_base::floatfield );
    out.precision( 16 );
    out << "!Triangles";

    const unsigned int numVertices = vertices_.size();
    // print information about vertices and elements
    // to header to have an easy check
    //   out << "  ( noVertices = " << numVertices;
    //   out << " | noElements = " << elements_.size() << " )" << std :: endl;

    out << std::endl;

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

    out << boundaryIds_.size() << std :: endl;
    typedef typename std :: vector< std :: pair< FaceType, int > > :: iterator
    BoundaryIdIteratorType;
    const BoundaryIdIteratorType endB = boundaryIds_.end();
    for( BoundaryIdIteratorType it = boundaryIds_.begin(); it != endB; ++it )
    {
      const std :: pair< FaceType, int > &boundaryId = *it;
      out << "-" << boundaryId.second ;
      for( unsigned int i = 0; i < numFaceCorners; ++i )
        out << "  " << boundaryId.first[ i ];
      out << std :: endl;
    }

    // for( unsigned int i = 0; i < numVertices; ++i )
    //  out << i << "  -1" << std :: endl;
    out.close();

    vertices_.clear();
    elements_.clear();
    boundaryIds_.clear();

    Grid *grid = new Grid( filename_ );
    if( removeGeneratedFile_ )
      remove( filename_.c_str() );
    return grid;
  }


  template< template< int, int > class ALUGrid >
  inline std::string
  ALU2dGridFactory<ALUGrid>::temporaryFileName ()
  {
    char filetemp[ FILENAME_MAX ];
    std :: strcpy( filetemp, "ALU2dGrid.XXXXXX" );
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
