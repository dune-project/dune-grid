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

  template< template< int, int > class ALUGrid >
  ALU3dGridFactory< ALUGrid >
  :: ALU3dGridFactory ( const MPICommunicatorType &communicator )
    : filename_("") , communicator_( communicator )
  {
#if ALU3DGRID_PARALLEL
    MPI_Comm_rank( communicator, &rank_ );
#endif
  }


  template< template< int, int > class ALUGrid >
  ALU3dGridFactory< ALUGrid >
  :: ALU3dGridFactory ( const std::string filename,
                        const MPICommunicatorType &communicator )
    : filename_(filename) , communicator_( communicator )
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
                     const std :: vector< unsigned int > &vertices )
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

    ElementType element;
    for( unsigned int i = 0; i < numCorners; ++i )
    {
      const unsigned int j = ElementTopologyMappingType :: dune2aluVertex( i );
      element[ j ] = vertices[ i ];
    }
    elements_.push_back( element );
  }


  template< template< int, int > class ALUGrid >
  void ALU3dGridFactory< ALUGrid >
  :: insertBoundary ( const GeometryType &geometry,
                      const std :: vector< unsigned int > &vertices,
                      int id )
  {
#if ALU3DGRID_PARALLEL
    if( rank_ != 0 )
      DUNE_THROW( GridError,
                  "ALU3dGridFactory allows insertion only for rank 0." );
#endif
    assertGeometryType( geometry );
    if( geometry.dim() != dimension-1 )
      DUNE_THROW( GridError, "Only 2-dimensional boundaries can be inserted "
                  "into a 3-dimensional ALUGrid." );
    if( vertices.size() != numFaceCorners )
      DUNE_THROW( GridError, "Wrong number of vertices." );

    std :: pair< FaceType, int > boundaryId;
    for( unsigned int i = 0; i < numFaceCorners; ++i )
    {
      const unsigned int j = FaceTopologyMappingType :: dune2aluVertex( i );
      boundaryId.first[ j ] = vertices[ i ];
    }
    boundaryId.second = id;
    boundaryIds_.push_back( boundaryId );
  }


  template< template< int, int > class ALUGrid >
  ALUGrid< 3, 3 > *ALU3dGridFactory< ALUGrid > :: createGrid ()
  {
#if ALU3DGRID_PARALLEL
    if( rank_ != 0 )
      return new GridType( communicator_ );
#endif

    std::string filename(filename_);
    if( filename == "" )
    {
      char filetemp[ FILENAME_MAX ];
      std :: strcpy( filetemp, "ALU3dGrid.XXXXXX" );
      mkstemp( filetemp );
      filename = filetemp;
    }
    else
    {
      // append filename by alugrid
      filename += ".ALU3dGrid";
    }

    std :: ofstream out( filename.c_str() );
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
    typedef typename std :: vector< ElementType > :: iterator
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
