// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// Of course, we would like to put this into the lib, but unfortunately HAVE_MPI
// is unpredictable.
//#include <config.h>

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>

#include <dune/grid/alugrid/3d/alu3dgridfactory.hh>

namespace Dune
{

  // GridFactory for ALUCubeGrid< 3, 3 >
  // -----------------------------------

  GridFactory< ALUSimplexGrid< 3, 3 > >
  :: GridFactory ( const MPICommunicatorType &communicator )
    : communicator_( communicator )
  {
#if ALU3DGRID_PARALLEL
    MPI_Comm_rank( communicator, &rank_ );
#endif
  }


  GridFactory< ALUSimplexGrid< 3, 3 > > :: ~GridFactory ()
  {}


  void GridFactory< ALUSimplexGrid< 3, 3 > >
  :: insertVertex ( const VertexType &pos )
  {
#if ALU3DGRID_PARALLEL
    if( rank_ != 0 )
      DUNE_THROW( GridError, "GridFactory< ALUSimplexGrid< 3, 3 > > allows "
                  "insertion only for rank 0." );
#endif
    vertices_.push_back( pos );
  }


  void GridFactory< ALUSimplexGrid< 3, 3 > >
  :: insertElement ( const GeometryType &type,
                     const std :: vector< unsigned int > &vertices )
  {
#if ALU3DGRID_PARALLEL
    if( rank_ != 0 )
      DUNE_THROW( GridError, "GridFactory< ALUSimplexGrid< 3, 3 > > allows "
                  "insertion only for rank 0." );
#endif
    if( !type.isSimplex() || (type.dim() != dimension) )
      DUNE_THROW( GridError, "ALUSimplexGrid supports only simplices." );
    if( vertices.size() != numCorners )
      DUNE_THROW( GridError, "Wrong number of vertices." );

    ElementType element;
    for( unsigned int i = 0; i < numCorners; ++i )
    {
      const unsigned int j = ElementTopologyMappingType :: dune2aluVertex( i );
      element[ j ] = vertices[ i ];
    }
    elements_.push_back( element );
  }


  void GridFactory< ALUSimplexGrid< 3, 3 > >
  :: insertBoundaryId ( const GeometryType &faceType,
                        const std :: vector< unsigned int > &faceVertices,
                        int id )
  {
#if ALU3DGRID_PARALLEL
    if( rank_ != 0 )
      DUNE_THROW( GridError, "GridFactory< ALUSimplexGrid< 3, 3 > > allows "
                  "insertion only for rank 0." );
#endif
    if( !faceType.isSimplex() || (faceType.dim() != dimension - 1) )
      DUNE_THROW( GridError, "ALUSimplexGrid supports only simplices." );
    if( faceVertices.size() != numFaceCorners )
      DUNE_THROW( GridError, "Wrong number of vertices." );

    std :: pair< FaceType, int > boundaryId;
    for( unsigned int i = 0; i < numFaceCorners; ++i )
    {
      const unsigned int j = FaceTopologyMappingType :: dune2aluVertex( i );
      boundaryId.first[ j ] = faceVertices[ i ];
    }
    boundaryId.second = id;
    boundaryIds_.push_back( boundaryId );
  }


  ALUSimplexGrid< 3, 3 > *GridFactory< ALUSimplexGrid< 3, 3 > > :: createGrid ()
  {
#if ALU3DGRID_PARALLEL
    if( rank_ != 0 )
      return new GridType( communicator_ );
#endif

    char filename[ FILENAME_MAX ];
    std :: strcpy( filename, "ALU3dSimplexGrid.XXXXXX" );
    mkstemp( filename );

    std :: ofstream out( filename );
    out << "!Tetrahedra" << std :: endl;

    const unsigned int numVertices = vertices_.size();

    out << numVertices << std :: endl;
    typedef std :: vector< VertexType > :: iterator VertexIteratorType;
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
    typedef std :: vector< ElementType > :: iterator
    ElementIteratorType;
    const ElementIteratorType endE = elements_.end();
    for( ElementIteratorType it = elements_.begin(); it != endE; ++it )
    {
      const ElementType &element = *it;
      out << element[ 0 ];
      for( unsigned int i = 1; i < numCorners; ++i )
        out << "  " << element[ i ];
      out << std :: endl;
    }

    out << boundaryIds_.size() << std :: endl;
    typedef std :: vector< std :: pair< FaceType, int > > :: iterator
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
    GridType *grid = new GridType( filename, communicator_ );
#else
    GridType *grid = new GridType( filename );
#endif

#ifdef NDEBUG
    remove( filename );
#endif
    return grid;
  }



  // GridFactory for ALUCubeGrid< 3, 3 >
  // -----------------------------------

  GridFactory< ALUCubeGrid< 3, 3 > >
  :: GridFactory ( const MPICommunicatorType &communicator )
    : communicator_( communicator )
  {
#if ALU3DGRID_PARALLEL
    MPI_Comm_rank( communicator, &rank_ );
#endif
  }


  GridFactory< ALUCubeGrid< 3, 3 > > :: ~GridFactory ()
  {}


  void GridFactory< ALUCubeGrid< 3, 3 > >
  :: insertVertex ( const VertexType &pos )
  {
#if ALU3DGRID_PARALLEL
    if( rank_ != 0 )
      DUNE_THROW( GridError, "GridFactory< ALUCubeGrid< 3, 3 > > allows "
                  "insertion only for rank 0." );
#endif
    vertices_.push_back( pos );
  }


  void GridFactory< ALUCubeGrid< 3, 3 > >
  :: insertElement ( const GeometryType &type,
                     const std :: vector< unsigned int > &vertices )
  {
#if ALU3DGRID_PARALLEL
    if( rank_ != 0 )
      DUNE_THROW( GridError, "GridFactory< ALUCubeGrid< 3, 3 > > allows "
                  "insertion only for rank 0." );
#endif
    if( !type.isCube() || (type.dim() != dimension) )
      DUNE_THROW( GridError, "ALUCubeGrid supports only cubes." );
    if( vertices.size() != numCorners )
      DUNE_THROW( GridError, "Wrong number of vertices." );

    ElementType element;
    for( unsigned int i = 0; i < numCorners; ++i )
    {
      const unsigned int j = ElementTopologyMappingType :: dune2aluVertex( i );
      element[ j ] = vertices[ i ];
    }
    elements_.push_back( element );
  }


  void GridFactory< ALUCubeGrid< 3, 3 > >
  :: insertBoundaryId ( const GeometryType &faceType,
                        const std :: vector< unsigned int > &faceVertices,
                        int id )
  {
#if ALU3DGRID_PARALLEL
    if( rank_ != 0 )
      DUNE_THROW( GridError, "GridFactory< ALUCubeGrid< 3, 3 > > allows "
                  "insertion only for rank 0." );
#endif
    if( !faceType.isCube() || (faceType.dim() != dimension - 1) )
      DUNE_THROW( GridError, "ALUCubeGrid supports only cubes." );
    if( faceVertices.size() != numFaceCorners )
      DUNE_THROW( GridError, "Wrong number of vertices." );

    std :: pair< FaceType, int > boundaryId;
    for( unsigned int i = 0; i < numFaceCorners; ++i )
    {
      const unsigned int j = FaceTopologyMappingType :: dune2aluVertex( i );
      boundaryId.first[ j ] = faceVertices[ i ];
    }
    boundaryId.second = id;
    boundaryIds_.push_back( boundaryId );
  }


  ALUCubeGrid< 3, 3 > *GridFactory< ALUCubeGrid< 3, 3 > > :: createGrid ()
  {
#if ALU3DGRID_PARALLEL
    if( rank_ != 0 )
      return new GridType( communicator_ );
#endif

    char filename[ FILENAME_MAX ];
    std :: strcpy( filename, "ALU3dCubeGrid.XXXXXX" );
    mkstemp( filename );

    std :: ofstream out( filename );
    out << "!Hexahedra" << std :: endl;

    const unsigned int numVertices = vertices_.size();

    out << numVertices << std :: endl;
    typedef std :: vector< VertexType > :: iterator VertexIteratorType;
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
    typedef std :: vector< ElementType > :: iterator ElementIteratorType;
    const ElementIteratorType endE = elements_.end();
    for( ElementIteratorType it = elements_.begin(); it != endE; ++it )
    {
      const ElementType &element = *it;
      out << element[ 0 ];
      for( unsigned int i = 1; i < numCorners; ++i )
        out << "  " << element[ i ];
      out << std :: endl;
    }

    out << boundaryIds_.size() << std :: endl;
    typedef std :: vector< std :: pair< FaceType, int > > :: iterator
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
    GridType *grid = new GridType( filename, communicator_ );
#else
    GridType *grid = new GridType( filename );
#endif

#ifdef NDEBUG
    remove( filename );
#endif
    return grid;
  }

}
