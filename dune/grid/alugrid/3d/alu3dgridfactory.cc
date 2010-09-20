// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU3DGRID_FACTORY_CC
#define DUNE_ALU3DGRID_FACTORY_CC

#if COMPILE_ALUGRID_INLINE == 0
#include <config.h>
#endif

#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>

#include <dune/grid/alugrid/3d/alu3dgridfactory.hh>

#ifdef ENABLE_ALUGRID

#if COMPILE_ALUGRID_INLINE
#define alu_inline inline
#else
#define alu_inline
#endif

namespace Dune
{

  template< class ALUGrid >
  alu_inline
  ALU3dGridFactory< ALUGrid > :: ~ALU3dGridFactory ()
  {}


  template< class ALUGrid >
  alu_inline
  void ALU3dGridFactory< ALUGrid > :: insertVertex ( const VertexType &pos )
  {
    if( ! allowGridGeneration_ )
      DUNE_THROW( GridError, "ALU3dGridFactory allows insertion only for rank 0." );

    vertices_.push_back( pos );
  }


  template< class ALUGrid >
  alu_inline
  void ALU3dGridFactory< ALUGrid >
  :: insertElement ( const GeometryType &geometry,
                     const std::vector< unsigned int > &vertices )
  {
    if( ! allowGridGeneration_ )
      DUNE_THROW( GridError, "ALU3dGridFactory allows insertion only for rank 0." );

    assertGeometryType( geometry );
    if( geometry.dim() != dimension )
      DUNE_THROW( GridError, "Only 3-dimensional elements can be inserted "
                  "into a 3-dimensional ALUGrid." );
    if( vertices.size() != numCorners )
      DUNE_THROW( GridError, "Wrong number of vertices." );

    elements_.push_back( vertices );
  }


  template< class ALUGrid >
  alu_inline
  void ALU3dGridFactory< ALUGrid >
  :: insertBoundary ( const GeometryType &geometry,
                      const std::vector< unsigned int > &vertices,
                      const int id )
  {
    if( ! allowGridGeneration_ )
      DUNE_THROW( GridError, "ALU3dGridFactory allows insertion only for rank 0." );

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


  template< class ALUGrid >
  alu_inline
  void ALU3dGridFactory< ALUGrid >
  ::insertBoundary ( const int element, const int face, const int id )
  {
    if( ! allowGridGeneration_ )
      DUNE_THROW( GridError, "ALU3dGridFactory allows insertion only for rank 0." );

    if( (element < 0) || (element >= (int)elements_.size()) )
      DUNE_THROW( RangeError, "ALU3dGridFactory::insertBoundary: invalid element index given." );

    std::pair< FaceType, int > boundaryId;
    generateFace( elements_[ element ], face, boundaryId.first );
    boundaryId.second = id;
    boundaryIds_.push_back( boundaryId );
  }

  template< class ALUGrid >
  alu_inline
  void ALU3dGridFactory< ALUGrid > ::
  insertBoundaryProjection( const DuneBoundaryProjectionType& bndProjection )
  {
    if( globalProjection_ )
      DUNE_THROW(InvalidStateException,"You can only insert one globalProjection");

    globalProjection_ = &bndProjection;
  }


  template< class ALUGrid >
  alu_inline
  void ALU3dGridFactory< ALUGrid > ::
  insertBoundaryProjection ( const GeometryType &type,
                             const std::vector< unsigned int > &vertices,
                             const DuneBoundaryProjectionType *projection )
  {
    if( (int)type.dim() != dimension-1 )
      DUNE_THROW( GridError, "Inserting boundary face of wrong dimension: " << type.dim() );
    assert( type.isCube() || type.isSimplex() );

    FaceType faceId;
    copyAndSort( vertices, faceId );

    if( vertices.size() != numFaceCorners )
      DUNE_THROW( GridError, "Wrong number of face vertices passed: " << vertices.size() << "." );

    if( boundaryProjections_.find( faceId ) != boundaryProjections_.end() )
      DUNE_THROW( GridError, "Only one boundary projection can be attached to a face." );
    boundaryProjections_[ faceId ] = projection;
  }

#if 0
  template< class ALUGrid >
  void ALU3dGridFactory< ALUGrid > ::
  insertBoundarySegment ( const std::vector< unsigned int >& vertices )
  {
    DUNE_THROW( NotImplemented, "insertBoundarySegment with a single argument" );
  }

  template< class ALUGrid >
  void ALU3dGridFactory< ALUGrid > ::
  insertBoundarySegment ( const std::vector< unsigned int >& vertices,
                          const shared_ptr<BoundarySegment<3,3> >& boundarySegment )
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
        DUNE_THROW(GridError,"BoundarySegment does not map face vertices to face vertices.");
    }
#endif
  }
#endif


  template< class ALUGrid >
  alu_inline
  void ALU3dGridFactory< ALUGrid >
  ::insertFaceTransformation ( const WorldMatrix &matrix, const WorldVector &shift )
  {
    faceTransformations_.push_back( Transformation( matrix, shift ) );
  }


  template< class ALUGrid >
  alu_inline
  ALUGrid *ALU3dGridFactory< ALUGrid >::createGrid ()
  {
    return createGrid( true, true, "" );
  }

  template< class ALUGrid >
  alu_inline
  ALUGrid *ALU3dGridFactory< ALUGrid >
  ::createGrid ( const bool addMissingBoundaries, const std::string dgfName )
  {
    return createGrid( addMissingBoundaries, true, dgfName );
  }

  template< class ALUGrid >
  alu_inline
  ALUGrid *ALU3dGridFactory< ALUGrid >
  ::createGrid ( const bool addMissingBoundaries, bool temporary, const std::string name )
  {
    typedef typename std :: vector< std :: pair< FaceType, int > > :: iterator BoundaryIdIteratorType;
    BoundaryProjectionVector* bndProjections = 0;

    if( allowGridGeneration_ )
    {
      correctElementOrientation();
      numFacesInserted_ = boundaryIds_.size();
      if( addMissingBoundaries || !faceTransformations_.empty() )
        recreateBoundaryIds();

      // if dump file should be written
      if( ! temporary )
      {
        std::string filename ( name );

        std::ofstream out( filename.c_str() );
        out.setf( std::ios_base::scientific, std::ios_base::floatfield );
        out.precision( 16 );
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

        out << (boundaryIds_.size() + periodicBoundaries_.size()) << std :: endl;
        const BoundaryIdIteratorType endB = boundaryIds_.end();
        for( BoundaryIdIteratorType it = boundaryIds_.begin(); it != endB; ++it )
        {
          const std::pair< FaceType, int > &boundaryId = *it;
          out << "-" << boundaryId.second << "  " << numFaceCorners;

          for( unsigned int i = 0; i < numFaceCorners; ++i )
            out << "  " << boundaryId.first[ i ];
          out << std::endl;
        }
        const typename PeriodicBoundaryVector::iterator endP = periodicBoundaries_.end();
        for( typename PeriodicBoundaryVector::iterator it = periodicBoundaries_.begin(); it != endP; ++it )
        {
          const std::pair< FaceType, FaceType > &facePair = *it;
          out << "-20  " << (2*numFaceCorners);
          for( unsigned int i = 0; i < numFaceCorners; ++i )
            out << "  " << facePair.first[ numFaceCorners == 3 ? (3 - i) % 3 : i ];
          for( unsigned int i = 0; i < numFaceCorners; ++i )
            out << "  " << facePair.second[ numFaceCorners == 3 ? (3 - i) % 3 : i ];
          out << std::endl;
        }

        for( unsigned int i = 0; i < numVertices; ++i )
          out << i << "  -1" << std :: endl;
        out.close();
      }

      const size_t boundarySegments = boundaryIds_.size();

      const size_t bndProjectionSize = boundaryProjections_.size();
      if( bndProjectionSize > 0 )
      {
        // the memory is freed by the grid on destruction
        bndProjections = new BoundaryProjectionVector( boundarySegments,
                                                       (DuneBoundaryProjectionType*) 0 );
        const BoundaryIdIteratorType endB = boundaryIds_.end();
        int segmentIndex = 0;
        for( BoundaryIdIteratorType it = boundaryIds_.begin(); it != endB; ++it, ++segmentIndex )
        {
          // generate boundary segment pointer
          FaceType faceId ( (*it).first);
          std::sort( faceId.begin(), faceId.end() );

          const DuneBoundaryProjectionType* projection = boundaryProjections_[ faceId ];

          // if no projection given we use global projection, otherwise identity
          if( ! projection && globalProjection_ )
          {
            typedef BoundaryProjectionWrapper< dimensionworld > ProjectionWrapperType;
            // we need to wrap the global projection because of
            // delete in desctructor of ALUGrid
            projection = new ProjectionWrapperType( *globalProjection_ );
            assert( projection );
          }

          // copy pointer
          (*bndProjections)[ segmentIndex ] = projection;
        }

      } // end ! temporary
    } // end rank == 0

    // free memory
    boundaryProjections_.clear();

    // if we have a vector reset global projection
    // because empty positions are filled with global projection anyway
    if( bndProjections ) globalProjection_ = 0;

    // ALUGrid is taking ownership of bndProjections
    // and is going to delete this pointer
    Grid* grid = createGridObj( bndProjections , name );
    assert( grid );

    // remove pointers
    globalProjection_ = 0;
    // is removed by grid instance
    bndProjections    = 0;

    // insert grid using ALUGrid macro grid builder
    if( allowGridGeneration_ )
    {
      ALU3DSPACE MacroGridBuilder mgb ( grid->getBuilder()
#ifdef ALUGRID_VERTEX_PROJECTION
                                        , grid->vertexProjection()
#endif
                                        );

      // now start inserting grid
      const int vxSize = vertices_.size();
      for( int vxIdx = 0; vxIdx < vxSize ; ++vxIdx )
      {
        const VertexType &vertex = vertices_[ vxIdx ];
        // insert vertex
        mgb.InsertUniqueVertex( vertex[ 0 ], vertex[ 1 ], vertex[ 2 ],  vxIdx );
      }

      typedef typename ElementVector::iterator ElementIteratorType;
      const ElementIteratorType endE = elements_.end();
      for( ElementIteratorType it = elements_.begin(); it != endE; ++it )
      {
        if( elementType == hexa )
        {
          int element[ 8 ];
          for( unsigned int i = 0; i < 8; ++i )
          {
            const unsigned int j = ElementTopologyMappingType::dune2aluVertex( i );
            element[ j ] = (*it)[ i ];
          }
          mgb.InsertUniqueHexa( element );
        }
        else if( elementType == tetra )
        {
          int element[ 4 ];
          for( unsigned int i = 0; i < 4; ++i )
          {
            const unsigned int j = ElementTopologyMappingType::dune2aluVertex( i );
            element[ j ] = (*it)[ i ];
          }
          mgb.InsertUniqueTetra( element );
        }
        else
          DUNE_THROW( GridError, "Invalid element type");
      }

      const BoundaryIdIteratorType endB = boundaryIds_.end();
      for( BoundaryIdIteratorType it = boundaryIds_.begin(); it != endB; ++it )
      {
        const std :: pair< FaceType, int > &boundaryId = *it;
        ALU3DSPACE Gitter::hbndseg::bnd_t bndType = (ALU3DSPACE Gitter::hbndseg::bnd_t ) boundaryId.second;

        if( elementType == hexa )
        {
          int bndface[ 4 ];
          for( unsigned int i = 0; i < numFaceCorners; ++i )
          {
            bndface[ i ] = boundaryId.first[ i ];
          }
          mgb.InsertUniqueHbnd4( bndface, bndType );
        }
        else if( elementType == tetra )
        {
          int bndface[ 3 ];
          for( unsigned int i = 0; i < numFaceCorners; ++i )
          {
            bndface[ i ] = boundaryId.first[ i ];
          }
          mgb.InsertUniqueHbnd3( bndface, bndType );
        }
        else
          DUNE_THROW( GridError, "Invalid element type");
      }
      const typename PeriodicBoundaryVector::iterator endP = periodicBoundaries_.end();
      for( typename PeriodicBoundaryVector::iterator it = periodicBoundaries_.begin(); it != endP; ++it )
      {
        const std::pair< FaceType, FaceType > &facePair = *it;
        if( elementType == hexa )
        {
          int perel[ 8 ];
          for( unsigned int i = 0; i < numFaceCorners; ++i )
          {
            perel[ i+0 ] = facePair.first[ i ];
            perel[ i+4 ] = facePair.second[ i ];
          }
          mgb.InsertUniquePeriodic4( perel );
        }
        else if( elementType == tetra )
        {
          int perel[ 6 ];
          for( unsigned int i = 0; i < 3; ++i )
          {
            perel[ i+0 ] = facePair.first[ (3 - i) % 3 ];
            perel[ i+3 ] = facePair.second[ (3 - i) % 3 ];
          }
          mgb.InsertUniquePeriodic3( perel );
        }
        else
          DUNE_THROW( GridError, "Invalid element type" );
      }

    }

    // clear vertices
    vertices_ = VertexVector();
    // clear elements
    elements_ = ElementVector();
    // free memory
    boundaryIds_.clear();

#ifdef ALUGRID_EXPORT_MACROGRID_CHANGES
    if( realGrid_ )
    {
      // make changes in macro grid known in every partition
      grid->duneNotifyMacroGridChanges();
    }
#else
    if( grid->comm().size() > 1 )
      DUNE_THROW(NotImplemented,"ALUGrid factory not working in parallel right now!");
#endif // #ifdef ALUGRID_EXPORT_MACROGRID_CHANGES

    // reset wasRefined flags
    grid->postAdapt();
    // update additional information on grid
    grid->calcExtras();

    return grid;
  }


  template< class ALUGrid >
  alu_inline
  void ALU3dGridFactory< ALUGrid >
  ::generateFace ( const ElementType &element, const int f, FaceType &face )
  {
    typedef ElementTopologyMapping< elementType > ElementTopologyMappingType;
    const int falu = ElementTopologyMappingType :: generic2aluFace( f );
    for( unsigned int i = 0; i < numFaceCorners; ++i )
    {
      const int j = ElementTopologyMappingType :: faceVertex( falu, i );
      const int k = ElementTopologyMappingType :: alu2genericVertex( j );
      face[ i ] = element[ k ];
    }
  }


  template< class ALUGrid >
  alu_inline
  void
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
          std::swap( element[ i ], element[ i+4 ] );
      }
      else
        std::swap( element[ 2 ], element[ 3 ] );
    } // end of loop over all elements
  }


  template< class ALUGrid >
  alu_inline
  bool ALU3dGridFactory< ALUGrid >
  ::identifyFaces ( const Transformation &transformation, const FaceType &key1, const FaceType &key2 )
  {
    WorldVector w = transformation.evaluate( vertices_[ key1[ 0 ] ] );
    int org = -1;
    for( unsigned int i = 0; i < numFaceCorners; ++i )
    {
      if( (w - vertices_[ key2[ i ] ]).two_norm() < 1e-6 )
        org = i;
    }
    if( org < 0 )
      return false;

    FaceType key0;
    key0[ 0 ] = key2[ org ];
    for( unsigned int i = 1; i < numFaceCorners; ++i )
    {
      w = transformation.evaluate( vertices_[ key1[ i ] ] );
      const int j = ((org+numFaceCorners)-i) % numFaceCorners;
      if( (w - vertices_[ key2[ j ] ]).two_norm() >= 1e-6 )
        return false;
      key0[ i ] = key2[ j ];
    }

    periodicBoundaries_.push_back( std::make_pair( key0, key1 ) );

    return true;
  }


  template< class ALUGrid >
  alu_inline
  void ALU3dGridFactory< ALUGrid >
  ::searchPeriodicNeighbor ( FaceMap &faceMap, const typename FaceMap::iterator &pos )
  {
    typedef typename FaceTransformationVector::const_iterator TrafoIterator;
    typedef typename FaceMap::iterator FaceMapIterator;

    if( !faceTransformations_.empty() )
    {
      FaceType key1;
      generateFace( pos->second, key1 );

      const FaceMapIterator fend = faceMap.end();
      for( FaceMapIterator fit = faceMap.begin(); fit != fend; ++fit )
      {
        FaceType key2;
        generateFace( fit->second, key2 );

        const TrafoIterator trend = faceTransformations_.end();
        for( TrafoIterator trit = faceTransformations_.begin(); trit != trend; ++trit )
        {
          if( identifyFaces( *trit, key1, key2 ) || identifyFaces( *trit, key2, key1 ) )
          {
            faceMap.erase( fit );
            faceMap.erase( pos );
            return;
          }
        }
      }
    }
  }


  template< class ALUGrid >
  alu_inline
  void ALU3dGridFactory< ALUGrid >
  ::reinsertBoundary ( const FaceMap &faceMap, const typename FaceMap::const_iterator &pos, const int id )
  {
    insertBoundary( pos->second.first, pos->second.second, id );
  }


  template< class ALUGrid >
  alu_inline
  void ALU3dGridFactory< ALUGrid >
  ::recreateBoundaryIds ( const int defaultId )
  {
    typedef typename FaceMap::iterator FaceIterator;
    FaceMap faceMap;

#if 0
    const GeometryType faceGeo( elementType == tetra
                                ? GeometryType::simplex : GeometryType::cube,
                                dimension-1 );
#endif

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
        {
          faceMap.insert( std::make_pair( key, SubEntity( n, face ) ) );
          searchPeriodicNeighbor( faceMap, faceMap.find( key ) );
        }
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
      //DUNE_THROW( GridError, "Inserted boundary segment is not part of the boundary." );
      reinsertBoundary( faceMap, pos, bndIt->second );
      faceMap.erase( pos );
    }

    // add all new boundaries (with defaultId)
    const FaceIterator faceEnd = faceMap.end();
    for( FaceIterator faceIt = faceMap.begin(); faceIt != faceEnd; ++faceIt )
      reinsertBoundary( faceMap, faceIt, defaultId );
  }

#if COMPILE_ALUGRID_LIB
  template class ALU3dGridFactory< ALUCubeGrid< 3, 3 > >;
  template class ALU3dGridFactory< ALUSimplexGrid< 3, 3 > >;
#endif
}
#endif // end ENABLE_ALUGRID
#undef alu_inline

#endif // end DUNE_ALU3DGRID_FACTORY_CC
