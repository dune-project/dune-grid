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
#include <dune/grid/alugrid/2d/alugrid.hh>

namespace Dune
{

  template< template< int, int > class ALUGrid, int dimw >
  void ALU2dGridFactory< ALUGrid, dimw >::insertVertex ( const VertexType &pos )
  {
    vertices_.push_back( pos );
  }


  template< template< int, int > class ALUGrid, int dimw >
  void ALU2dGridFactory< ALUGrid, dimw >
  ::insertElement ( const GeometryType &geometry,
                    const std::vector< unsigned int > &vertices )
  {
    switch( elementType )
    {
    case ALU2DSPACE triangle :
      if( !geometry.isTriangle() )
        DUNE_THROW( GridError, "Only triangles can be inserted into "
                    "ALUGrid< 2, " << dimw << ", triangle >." );
      break;
    case ALU2DSPACE quadrilateral :
      if( !geometry.isCube() )
        DUNE_THROW( GridError, "Only cubes can be inserted into "
                    "ALUGrid< 2, " << dimw << ", quadrilateral >." );
      break;
    default :
      assert( geometry.isSimplex() || geometry.isCube() );
    }
    if( (geometry.isSimplex() && (vertices.size() != 3))
        || (geometry.isCube() && (vertices.size() != 4)) )
      DUNE_THROW( GridError, "Wrong number of vertices." );

    elements_.push_back( vertices );
  }


  template< template< int, int > class ALUGrid, int dimw >
  void ALU2dGridFactory< ALUGrid, dimw >
  ::insertBoundary ( const GeometryType &geometry,
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

    if( id <= 0 )
      DUNE_THROW( GridError, "Boundary ids must be positive." );

    std::pair< FaceType, int > boundaryId;
    for( unsigned int i = 0; i < numFaceCorners; ++i )
      boundaryId.first[ i ] = vertices[ i ];
    boundaryId.second = -id; // store negative id
    boundaryIds_.push_back( boundaryId );
  }


  template< template< int, int > class ALUGrid, int dimw >
  void ALU2dGridFactory< ALUGrid, dimw >
  ::insertBoundary ( const int element, const int face, const int id )
  {
    if( (element < 0) || (element >= (int)elements_.size()) )
      DUNE_THROW( RangeError, "ALU2dGridFactory::insertBoundary: invalid element index given." );

    std::pair< FaceType, int > boundaryId;
    generateFace( elements_[ element ], face, boundaryId.first );
    boundaryId.second = -id; // store negative id
    boundaryIds_.push_back( boundaryId );
  }


  template< template< int, int > class ALUGrid, int dimw >
  void ALU2dGridFactory< ALUGrid, dimw >
  ::insertBoundaryProjection ( const DuneBoundaryProjectionType &bndProjection )
  {
    if( globalProjection_ )
      DUNE_THROW(InvalidStateException,"You can only insert one globalProjection");

    globalProjection_ = &bndProjection;
  }


  template< template< int, int > class ALUGrid, int dimw >
  void ALU2dGridFactory< ALUGrid, dimw >
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


  template< template< int, int > class ALUGrid, int dimw >
  void ALU2dGridFactory< ALUGrid, dimw >
  ::insertBoundarySegment ( const std::vector< unsigned int >& vertices )
  {
    FaceType faceId;
    copyAndSort( vertices, faceId );

    if( vertices.size() != numFaceCorners )
      DUNE_THROW( GridError, "Wrong number of face vertices passed: " << vertices.size() << "." );

    if( boundaryProjections_.find( faceId ) != boundaryProjections_.end() )
      DUNE_THROW( GridError, "Only one boundary projection can be attached to a face." );

    boundaryProjections_[ faceId ] = 0;
  }


  template< template< int, int > class ALUGrid, int dimw >
  void ALU2dGridFactory< ALUGrid, dimw >
  ::insertBoundarySegment ( const std::vector< unsigned int > &vertices,
                            const shared_ptr< BoundarySegment< 2, dimw > > &boundarySegment )
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
      if( (global - coords[ i ]).two_norm() > 1e-6 )
      {
        DUNE_THROW(GridError,"BoundarySegment does not map face vertices to face vertices.");
      }
    }
#endif
  }


  template< template< int, int > class ALUGrid, int dimw >
  void ALU2dGridFactory< ALUGrid, dimw >
  :: insertFaceTransformation ( const WorldMatrix &matrix, const WorldVector &shift )
  {
    faceTransformations_.push_back( Transformation( matrix, shift ) );
  }


  template< template< int, int > class ALUGrid, int dimw >
  ALUGrid< 2, dimw > *ALU2dGridFactory< ALUGrid, dimw >
  ::createGrid ( const bool addMissingBoundaries, bool temporary, const std::string name )
  {
    numFacesInserted_ = boundaryIds_.size();
    if( addMissingBoundaries || !faceTransformations_.empty() )
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
      const std::pair< FaceType, int > &boundaryId = *it;
      out << boundaryId.second;
      for( unsigned int i = 0; i < numFaceCorners; ++i )
        out << "  " << boundaryId.first[ i ];
      if( boundaryId.second == periodicBndId )
      {
        typedef typename PeriodicNeighborMap::const_iterator Iterator;
        Iterator pos = periodicNeighborMap_.find( boundaryId.first );
        if( pos != periodicNeighborMap_.end() )
          out << "  " << pos->second;
        //else
        //  DUNE_THROW( InvalidStateException, "Periodic Neighbor not found." );
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

#ifdef ALUGRID_NOTEMPFILE_2D
    std::istream& inFile = temp;
#endif

    // if we have a vector reset global projection
    // because empty positions are filled with global projection anyway
    if( bndProjections ) globalProjection_ = 0;

    // ALUGrid is taking ownership of the bndProjections pointer
    Grid *grid =
#ifdef ALUGRID_NOTEMPFILE_2D
      ( temporary ) ? new Grid( filename, inFile, globalProjection_ , bndProjections, grdVerbose_ ) :
#endif
      new Grid( filename, globalProjection_ , bndProjections, grdVerbose_ );

#ifndef ALUGRID_NOTEMPFILE_2D
    if( temporary )
      std::remove( filename.c_str() );
#endif

    // remove pointers
    globalProjection_ = 0;
    // is removed by the grid
    bndProjections    = 0;

    return grid;
  }


  template< template< int, int > class ALUGrid, int dimw >
  std::string
  ALU2dGridFactory<ALUGrid,dimw>::temporaryFileName ( const std::string& dgfName )
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


  template< template< int, int > class ALUGrid, int dimw >
  inline void ALU2dGridFactory<ALUGrid,dimw>
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


  template< template< int, int > class ALUGrid, int dimw >
  typename ALU2dGridFactory< ALUGrid, dimw >::FaceMap::iterator
  ALU2dGridFactory< ALUGrid, dimw >
  ::findPeriodicNeighbor ( const FaceType &key, FaceMap &faceMap ) const
  {
    typedef typename FaceTransformationVector::const_iterator TrafoIterator;
    typedef typename FaceMap::iterator FaceIterator;

    const WorldVector v[ 2 ] = { vertices_[ key[ 0 ] ], vertices_[ key[ 1 ] ] };

    // Note This should be done using a kd-tree
    const TrafoIterator trend = faceTransformations_.end();
    for( TrafoIterator trit = faceTransformations_.begin(); trit != trend; ++trit )
    {
      WorldVector w[ 2 ];
      w[ 0 ] = trit->evaluate( v[ 0 ] );
      w[ 1 ] = trit->evaluate( v[ 1 ] );

      const FaceIterator fend = faceMap.end();
      for( FaceIterator fit = faceMap.begin(); fit != fend; ++fit )
      {
        const WorldVector vv[ 2 ] = { vertices_[ fit->first[ 0 ] ], vertices_[ fit->first[ 1 ] ] };

        if( ((vv[ 0 ] - w[ 0 ]).two_norm() < 1e-8) && ((vv[ 1 ] - w[ 1 ]).two_norm() < 1e-8) )
          return fit;
        if( ((vv[ 0 ] - w[ 1 ]).two_norm() < 1e-8) && ((vv[ 1 ] - w[ 0 ]).two_norm() < 1e-8) )
          return fit;

        WorldVector ww[ 2 ];
        ww[ 0 ] = trit->evaluate( vv[ 0 ] );
        ww[ 1 ] = trit->evaluate( vv[ 1 ] );

        if( ((v[ 0 ] - ww[ 0 ]).two_norm() < 1e-8) && ((v[ 1 ] - ww[ 1 ]).two_norm() < 1e-8) )
          return fit;
        if( ((v[ 0 ] - ww[ 1 ]).two_norm() < 1e-8) && ((v[ 1 ] - ww[ 0 ]).two_norm() < 1e-8) )
          return fit;
      }
    }
    return faceMap.end();
  }


  template< template< int, int > class ALUGrid, int dimw >
  void ALU2dGridFactory< ALUGrid, dimw >::recreateBoundaryIds ( const int defaultId )
  {
    typedef typename FaceMap::iterator FaceIterator;
    typedef typename PeriodicNeighborMap::const_iterator PeriodicNbIterator;

    if( defaultId <= 0 )
      DUNE_THROW( GridError, "Boundary ids must be positive." );

    FaceMap faceMap;

    const unsigned int numElements = elements_.size();
    for( unsigned int n = 0; n < numElements; ++n )
    {
      for( unsigned int face = 0; face < elements_[n].size(); ++face )
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
      int id = -bndIt->second;

      const FaceIterator pos = faceMap.find( key );
      if( pos == faceMap.end() )
        continue;

      const PeriodicNbIterator nb = periodicNeighborMap_.find( key );
      if( nb != periodicNeighborMap_.end() )
      {
        id = -periodicBndId;
        periodicNeighborMap_[ boundaryIds_[ nb->second ].first ] = boundaryIds_.size();
      }
      else
      {
        const FaceIterator nbpos = findPeriodicNeighbor( key, faceMap );
        if( nbpos != faceMap.end() )
        {
          id = -periodicBndId;
          if( periodicNeighborMap_.find( nbpos->first ) != periodicNeighborMap_.end() )
            DUNE_THROW( GridError, "One boundary segment can only identified with one other." );
          periodicNeighborMap_[ nbpos->first ] = boundaryIds_.size();
        }
      }

      insertBoundary( pos->second.first, pos->second.second, id );
      faceMap.erase( pos );
    }

    // add all new boundaries (with defaultId)
    const FaceIterator faceEnd = faceMap.end();
    for( FaceIterator faceIt = faceMap.begin(); faceIt != faceEnd; ++faceIt )
    {
      int id = defaultId;

      const PeriodicNbIterator nb = periodicNeighborMap_.find( faceIt->first );
      if( nb != periodicNeighborMap_.end() )
      {
        id = -periodicBndId;
        periodicNeighborMap_[ boundaryIds_[ nb->second ].first ] = boundaryIds_.size();
      }
      else
      {
        const FaceIterator nbpos = findPeriodicNeighbor( faceIt->first, faceMap );
        if( nbpos != faceMap.end() )
        {
          id = -periodicBndId;
          if( periodicNeighborMap_.find( nbpos->first ) != periodicNeighborMap_.end() )
            DUNE_THROW( GridError, "One boundary segment can only identified with one other." );
          periodicNeighborMap_[ nbpos->first ] = boundaryIds_.size();
        }
      }

      insertBoundary( faceIt->second.first, faceIt->second.second, id );
    }
  }

}



// Template Instantiation
// ----------------------

template class Dune::ALU2dGridFactory< Dune::ALUConformGrid, 2 >;
template class Dune::ALU2dGridFactory< Dune::ALUSimplexGrid, 2 >;
#ifdef ALUGRID_SURFACE_2D
template class Dune::ALU2dGridFactory< Dune::ALUCubeGrid, 2 >;

template class Dune::ALU2dGridFactory< Dune::ALUConformGrid, 3 >;
template class Dune::ALU2dGridFactory< Dune::ALUSimplexGrid, 3 >;
template class Dune::ALU2dGridFactory< Dune::ALUCubeGrid, 3 >;
#endif // #ifdef ALUGRID_SURFACE_2D
