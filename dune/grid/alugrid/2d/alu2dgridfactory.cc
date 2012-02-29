// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>

#include <dune/grid/alugrid/2d/alu2dgridfactory.hh>
#include <dune/grid/alugrid/common/declaration.hh>
#include <dune/grid/alugrid/2d/alugrid.hh>

namespace Dune
{

  template< class GridImp >
  void ALU2dGridFactory< GridImp >::insertVertex ( const VertexType &pos )
  {
    vertices_.push_back( pos );
  }


  template< class GridImp >
  void ALU2dGridFactory< GridImp >
  ::insertElement ( const GeometryType &geometry,
                    const std::vector< unsigned int > &vertices )
  {
    switch( elementType )
    {
    case ALU2DSPACE triangle :
      if( !geometry.isTriangle() )
        DUNE_THROW( GridError, "Only triangles can be inserted into "
                    "ALUGrid< 2, " << dimensionworld << ", triangle >." );
      break;
    case ALU2DSPACE quadrilateral :
      if( !geometry.isCube() )
        DUNE_THROW( GridError, "Only cubes can be inserted into "
                    "ALUGrid< 2, " << dimensionworld << ", quadrilateral >." );
      break;
    default :
      assert( geometry.isSimplex() || geometry.isCube() );
    }
    if( (geometry.isSimplex() && (vertices.size() != 3))
        || (geometry.isCube() && (vertices.size() != 4)) )
      DUNE_THROW( GridError, "Wrong number of vertices." );

    elements_.push_back( vertices );
  }


  template< class GridImp >
  void ALU2dGridFactory< GridImp >
  ::insertBoundary ( const GeometryType &geometry,
                     const std::vector< unsigned int > &vertices,
                     const int id )
  {
    // lines can be either cube or simplex
    assert( geometry.isSimplex() || geometry.isCube() );
    if( geometry.dim() != dimension-1 )
    {
      DUNE_THROW( GridError, "Only 1-dimensional boundaries can be inserted "
                  "into a 2-dimensional ALUGrid." );
    }
    if( vertices.size() != numFaceCorners )
      DUNE_THROW( GridError, "Wrong number of vertices." );

    if( id <= 0 )
      DUNE_THROW( GridError, "Boundary ids must be positive." );

    std::pair< FaceType, int > boundaryId;
    for( unsigned int i = 0; i < numFaceCorners; ++i )
      boundaryId.first[ i ] = vertices[ i ];
    boundaryId.second = -id; // store negative id
    boundaryIds_.push_back( boundaryId );
  }


  template< class GridImp >
  void ALU2dGridFactory< GridImp >
  ::insertBoundary ( const int element, const int face, const int id )
  {
    if( (element < 0) || (element >= (int)elements_.size()) )
      DUNE_THROW( RangeError, "ALU2dGridFactory::insertBoundary: invalid element index given." );

    std::pair< FaceType, int > boundaryId;
    generateFace( elements_[ element ], face, boundaryId.first );
    boundaryId.second = -id; // store negative id
    boundaryIds_.push_back( boundaryId );
  }


  template< class GridImp >
  void ALU2dGridFactory< GridImp >
  ::insertBoundaryProjection ( const DuneBoundaryProjectionType &bndProjection )
  {
    if( globalProjection_ )
      DUNE_THROW(InvalidStateException,"You can only insert one globalProjection");

    globalProjection_ = &bndProjection;
  }


  template< class GridImp >
  void ALU2dGridFactory< GridImp >
  ::insertBoundaryProjection ( const GeometryType &type,
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


  template< class GridImp >
  void ALU2dGridFactory< GridImp >
  ::insertBoundarySegment ( const std::vector< unsigned int >& vertices )
  {
    FaceType faceId;
    copyAndSort( vertices, faceId );

    if( vertices.size() != numFaceCorners )
      DUNE_THROW( GridError, "Wrong number of face vertices passed: " << vertices.size() << "." );

    if( boundaryProjections_.find( faceId ) != boundaryProjections_.end() )
      DUNE_THROW( GridError, "Only one boundary projection can be attached to a face." );

    boundaryProjections_[ faceId ] = 0;

    std::pair< FaceType, int > boundaryId;
    for( unsigned int i = 0; i < numFaceCorners; ++i )
      boundaryId.first[ i ] = vertices[ i ];
    boundaryId.second = -1; // some default
    boundaryIds_.push_back( boundaryId );
  }


  template< class GridImp >
  void ALU2dGridFactory< GridImp >
  ::insertBoundarySegment ( const std::vector< unsigned int > &vertices,
                            const shared_ptr< BoundarySegment< 2, dimensionworld > > &boundarySegment )
  {
    FaceType faceId;
    copyAndSort( vertices, faceId );

    if( vertices.size() != numFaceCorners )
      DUNE_THROW( GridError, "Wrong number of face vertices passed: " << vertices.size() << "." );

    if( boundaryProjections_.find( faceId ) != boundaryProjections_.end() )
      DUNE_THROW( GridError, "Only one boundary projection can be attached to a face." );

    const size_t numVx = vertices.size();
    GeometryType type;
    if( numVx == 3 )
      type.makeSimplex( dimension-1 );
    else
      type.makeCube( dimension-1 );

    std::vector< VertexType > coords( numVx );
    for( size_t i = 0; i < numVx; ++i )
    {
      // if this assertions is thrown vertices were not inserted at first
      assert( vertices_.size() > vertices[ i ] );

      // get global coordinate and copy it
      const VertexType &x = vertices_[ vertices[ i ] ];
      for( int j = 0; j < dimensionworld; ++j )
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
      if( (global - coords[ i ]).two_norm() >= epsilon_ )
      {
        DUNE_THROW(GridError,"BoundarySegment does not map face vertices to face vertices.");
      }
    }
#endif

    std::pair< FaceType, int > boundaryId;
    for( unsigned int i = 0; i < numFaceCorners; ++i )
      boundaryId.first[ i ] = vertices[ i ];
    boundaryId.second = -1; // some default
    boundaryIds_.push_back( boundaryId );
  }


  template< class GridImp >
  void ALU2dGridFactory< GridImp >
  :: insertFaceTransformation ( const WorldMatrix &matrix, const WorldVector &shift )
  {
    faceTransformations_.push_back( Transformation( matrix, shift ) );
  }


  template< class GridImp >
  GridImp  *ALU2dGridFactory< GridImp >
  ::createGrid ( const bool addMissingBoundaries, bool temporary, const std::string name )
  {
    numFacesInserted_ = boundaryIds_.size();
    if( addMissingBoundaries || !faceTransformations_.empty() )
      recreateBoundaryIds();

    std::string filename ( temporary ?
                           temporaryFileName( name ) :
                           name );

    std::ofstream outfile;
    std::stringstream temp;

    std::ostream* outptr = 0;
    if( temporary )
      outptr = & temp;
    else
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
      for( int i = 1; i < dimensionworld; ++i )
        out << " " << vertex[ i ];
      out << std :: endl;
    }

    out << elements_.size() << std :: endl;
    typedef typename ElementVector::iterator ElementIteratorType;
    const ElementIteratorType endE = elements_.end();
    for( ElementIteratorType it = elements_.begin(); it != endE; ++it )
    {
      out << (*it)[ 0 ];
      out << "  " << (*it)[ 1 ];
      if ( it->size() == 4 )
        out << "  " << (*it)[ 3 ];
      out << "  " << (*it)[ 2 ];
      out << std :: endl;
    }

    const size_t boundarySegments = boundaryIds_.size();
    out << boundarySegments << std::endl;
    typedef typename std::vector< std::pair< FaceType, int > >::iterator BoundaryIdIterator;

    BoundaryProjectionVector* bndProjections = 0;
    const size_t bndProjectionSize = boundaryProjections_.size();
    if( bndProjectionSize > 0 )
    {
      // the memory is freed by the grid on destruction
      bndProjections = new BoundaryProjectionVector( boundarySegments,
                                                     (DuneBoundaryProjectionType*)0 );
    }

    const BoundaryIdIterator endB = boundaryIds_.end();
    int segmentIndex = 0;
    for( BoundaryIdIterator it = boundaryIds_.begin(); it != endB; ++it, ++segmentIndex )
    {
      const std::pair< FaceType, int > &boundaryId = *it;
      out << boundaryId.second;
      for( unsigned int i = 0; i < numFaceCorners; ++i )
        out << "  " << boundaryId.first[ i ];
      if( boundaryId.second == periodicBndId )
      {
        FaceType key = boundaryId.first;
        std::sort( key.begin(), key.end() );

        typename PeriodicNeighborMap::const_iterator pos = periodicNeighborMap_.find( key );
        if( pos != periodicNeighborMap_.end() )
          out << "  " << pos->second;
        else
          DUNE_THROW( InvalidStateException, "Periodic Neighbor not found." );
      }
      out << std::endl;

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

    std::istream& inFile = temp;

    // if we have a vector reset global projection
    // because empty positions are filled with global projection anyway
    if( bndProjections ) globalProjection_ = 0;

    // ALUGrid is taking ownership of the bndProjections pointer
    Grid *grid = createGridObj( temporary, filename, inFile, bndProjections );

    if( temporary )
      std::remove( filename.c_str() );

    // remove pointers
    globalProjection_ = 0;
    // is removed by the grid
    bndProjections    = 0;

    return grid;
  }


  template< class GridImp >
  std::string
  ALU2dGridFactory< GridImp >::temporaryFileName ( const std::string &dgfName )
  {
    const std::string filename( (dgfName.empty() ? std::string( "ALU2dGrid" ) : dgfName) + ".XXXXXX" );
    char filetemp[ FILENAME_MAX ];
    std::strcpy( filetemp, filename.c_str() );
    const int fd = mkstemp( filetemp );
    if( fd < 0 )
      DUNE_THROW( IOError, "Unable to create temporary file." );
    close( fd );
    return std::string( filetemp );
  }


  template< class GridImp >
  inline void ALU2dGridFactory< GridImp >
  ::generateFace ( const ElementType &element, const int f, FaceType &face )
  {
    if (element.size() == 3)
      switch (f)
      {
      case 0 :
        face[0]=element[0]; face[1]=element[1];
        break;
      case 1 :
        face[0]=element[0]; face[1]=element[2];
        break;
      case 2 :
        face[0]=element[1]; face[1]=element[2];
        break;
      }
    else
      switch (f)
      {
      case 0 :
        face[0]=element[0]; face[1]=element[2];
        break;
      case 1 :
        face[0]=element[1]; face[1]=element[3];
        break;
      case 2 :
        face[0]=element[0]; face[1]=element[1];
        break;
      case 3 :
        face[0]=element[2]; face[1]=element[3];
        break;
      }
  }


  template< class GridImp >
  typename ALU2dGridFactory< GridImp >::FaceMap::const_iterator
  ALU2dGridFactory< GridImp >
  ::findPeriodicNeighbor ( const FaceMap &faceMap, const FaceType &key ) const
  {
    typedef typename FaceTransformationVector::const_iterator TrafoIterator;
    typedef typename FaceMap::const_iterator FaceMapIterator;

    const WorldVector v[ 2 ] = { vertices_[ key[ 0 ] ], vertices_[ key[ 1 ] ] };

    // Note This should be done using a kd-tree
    const TrafoIterator trend = faceTransformations_.end();
    for( TrafoIterator trit = faceTransformations_.begin(); trit != trend; ++trit )
    {
      WorldVector w[ 2 ];
      w[ 0 ] = trit->evaluate( v[ 0 ] );
      w[ 1 ] = trit->evaluate( v[ 1 ] );

      const FaceMapIterator fend = faceMap.end();
      for( FaceMapIterator fit = faceMap.begin(); fit != fend; ++fit )
      {
        const WorldVector vv[ 2 ] = { vertices_[ fit->first[ 0 ] ], vertices_[ fit->first[ 1 ] ] };

        if( ((vv[ 0 ] - w[ 0 ]).two_norm() < epsilon_) && ((vv[ 1 ] - w[ 1 ]).two_norm() < epsilon_) )
          return fit;
        if( ((vv[ 0 ] - w[ 1 ]).two_norm() < epsilon_) && ((vv[ 1 ] - w[ 0 ]).two_norm() < epsilon_) )
          return fit;

        WorldVector ww[ 2 ];
        ww[ 0 ] = trit->evaluate( vv[ 0 ] );
        ww[ 1 ] = trit->evaluate( vv[ 1 ] );

        if( ((v[ 0 ] - ww[ 0 ]).two_norm() < epsilon_) && ((v[ 1 ] - ww[ 1 ]).two_norm() < epsilon_) )
          return fit;
        if( ((v[ 0 ] - ww[ 1 ]).two_norm() < epsilon_) && ((v[ 1 ] - ww[ 0 ]).two_norm() < epsilon_) )
          return fit;
      }
    }
    return faceMap.end();
  }


  template< class GridImp >
  void ALU2dGridFactory< GridImp >
  ::reinsertBoundary ( const FaceMap &faceMap, const typename FaceMap::const_iterator &pos, const int id )
  {
    typedef typename PeriodicNeighborMap::const_iterator PeriodicNbIterator;
    typedef typename FaceMap::const_iterator FaceMapIterator;
    typedef std::pair< FaceType, int > BoundaryId;

    const size_t index = boundaryIds_.size();
    insertBoundary( pos->second.first, pos->second.second, -id );
    BoundaryId &boundaryId = boundaryIds_[ index ];

    const FaceType &key = pos->first;
    const PeriodicNbIterator nb = periodicNeighborMap_.find( key );
    if( nb != periodicNeighborMap_.end() )
    {
      boundaryId.second = periodicBndId;
      const BoundaryId &nbBoundaryId = boundaryIds_[ nb->second ];
      FaceType nbKey = nbBoundaryId.first;
      std::sort( nbKey.begin(), nbKey.end() );
      periodicNeighborMap_[ nbKey ] = index;
    }
    else
    {
      const FaceMapIterator nbPos = findPeriodicNeighbor( faceMap, key );
      if( nbPos != faceMap.end() )
      {
        boundaryId.second = periodicBndId;
        const FaceType &nbKey = nbPos->first;
        if( periodicNeighborMap_.find( nbKey ) != periodicNeighborMap_.end() )
          DUNE_THROW( GridError, "One boundary segment can only identified with one other." );
        periodicNeighborMap_[ nbKey ] = index;
      }
    }
  }


  template< class GridImp >
  void ALU2dGridFactory< GridImp >::recreateBoundaryIds ( const int defaultId )
  {
    typedef typename FaceMap::iterator FaceIterator;
    typedef typename PeriodicNeighborMap::const_iterator PeriodicNbIterator;

    if( defaultId <= 0 )
      DUNE_THROW( GridError, "Boundary ids must be positive." );

    FaceMap faceMap;

    const unsigned int numElements = elements_.size();
    for( unsigned int n = 0; n < numElements; ++n )
    {
      for( unsigned int face = 0; face < elements_[ n ].size(); ++face )
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
        DUNE_THROW( GridError, "Inserted boundary segment is not part of the boundary." );

      reinsertBoundary( faceMap, pos, bndIt->second );
      faceMap.erase( pos );
    }

    // add all new boundaries (with defaultId)
    const FaceIterator faceEnd = faceMap.end();
    for( FaceIterator faceIt = faceMap.begin(); faceIt != faceEnd; ++faceIt )
      reinsertBoundary( faceMap, faceIt, defaultId );
  }

}



// Template Instantiation
// ----------------------

template class Dune::ALU2dGridFactory< Dune::ALUConformGrid< 2, 2 > >;
template class Dune::ALU2dGridFactory< Dune::ALUSimplexGrid< 2, 2 > >;

template class Dune::ALU2dGridFactory< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming,    Dune::No_Comm > >;
template class Dune::ALU2dGridFactory< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::nonconforming, Dune::No_Comm > >;
#if HAVE_MPI
template class Dune::ALU2dGridFactory< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming,    MPI_Comm > >;
template class Dune::ALU2dGridFactory< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::nonconforming, MPI_Comm > >;
#endif // #if HAVE_MPI

#ifdef ALUGRID_SURFACE_2D
template class Dune::ALU2dGridFactory< Dune::ALUCubeGrid< 2, 2 > >;
template class Dune::ALU2dGridFactory< Dune::ALUGrid< 2, 2, Dune::cube, Dune::nonconforming, Dune::No_Comm > >;

template class Dune::ALU2dGridFactory< Dune::ALUConformGrid< 2, 3 > >;
template class Dune::ALU2dGridFactory< Dune::ALUSimplexGrid< 2, 3 > >;
template class Dune::ALU2dGridFactory< Dune::ALUCubeGrid< 2, 3 > >;

// new versions
template class Dune::ALU2dGridFactory< Dune::ALUGrid< 2, 3, Dune::simplex, Dune::conforming,    Dune::No_Comm > >;
template class Dune::ALU2dGridFactory< Dune::ALUGrid< 2, 3, Dune::simplex, Dune::nonconforming, Dune::No_Comm > >;
template class Dune::ALU2dGridFactory< Dune::ALUGrid< 2, 3, Dune::cube,    Dune::nonconforming, Dune::No_Comm > >;

#if HAVE_MPI
template class Dune::ALU2dGridFactory< Dune::ALUGrid< 2, 2, Dune::cube,    Dune::nonconforming, MPI_Comm > >;

template class Dune::ALU2dGridFactory< Dune::ALUGrid< 2, 3, Dune::simplex, Dune::conforming,    MPI_Comm > >;
template class Dune::ALU2dGridFactory< Dune::ALUGrid< 2, 3, Dune::simplex, Dune::nonconforming, MPI_Comm > >;
template class Dune::ALU2dGridFactory< Dune::ALUGrid< 2, 3, Dune::cube,    Dune::nonconforming, MPI_Comm > >;
#endif // #if HAVE_MPI

#endif // #ifdef ALUGRID_SURFACE_2D
