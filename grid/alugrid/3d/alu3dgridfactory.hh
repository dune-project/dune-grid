// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU3DGRID_FACTORY_HH
#define DUNE_ALU3DGRID_FACTORY_HH

#include <dune/common/fixedarray.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/alugrid.hh>
#include <dune/grid/io/file/dgfparser/entitykey.hh>

namespace Dune
{

  /** \brief Factory class for 3d ALUGrids */
  template< template< int, int > class ALUGrid >
  class ALU3dGridFactory
    : public GridFactoryInterface< ALUGrid< 3, 3 > >
  {
    typedef ALU3dGridFactory< ALUGrid > ThisType;
    typedef GridFactoryInterface< ALUGrid< 3, 3 > > BaseType;

  public:
    typedef ALUGrid< 3, 3 > GridType;

    typedef MPIHelper :: MPICommunicator MPICommunicatorType;

  private:
    typedef typename GridType :: ctype ctype;

    static const ALU3dGridElementType elementType = GridType :: elementType;

    dune_static_assert( (elementType == tetra || elementType == hexa),
                        "ALU3dGridFactory supports only grids containing "
                        "tetrahedrons or hexahedrons exclusively." );

    static const unsigned int dimension = GridType :: dimension;
    static const unsigned int dimensionworld = GridType :: dimensionworld;

    static const unsigned int numCorners
      = EntityCount< elementType > :: numVertices;
    static const unsigned int numFaceCorners
      = EntityCount< elementType > :: numVerticesPerFace;

    typedef ElementTopologyMapping< elementType > ElementTopologyMappingType;
    typedef FaceTopologyMapping< elementType > FaceTopologyMappingType;

    typedef FieldVector< ctype, dimensionworld > VertexType;
    typedef std::vector< unsigned int > ElementType;
    typedef array< unsigned int, numFaceCorners > FaceType;

  private:
    typedef std::vector< VertexType > VertexVector;
    typedef std::vector< ElementType > ElementVector;
    typedef std::vector< std::pair< FaceType, int > > BoundaryIdVector;

    const std::string filename_;
    MPICommunicatorType communicator_;
    bool removeGeneratedFile_;
#if ALU3DGRID_PARALLEL
    int rank_;
#endif
    VertexVector vertices_;
    ElementVector elements_;
    BoundaryIdVector boundaryIds_;

  public:
    /** \brief default constructor */
    explicit ALU3dGridFactory ( const MPICommunicatorType &communicator
                                  = MPIHelper :: getCommunicator(),
                                bool removeGeneratedFile = true );

    /** \brief constructor taking filename for temporary outfile */
    explicit ALU3dGridFactory ( const std::string &filename,
                                const MPICommunicatorType &communicator
                                  = MPIHelper :: getCommunicator() );

    /** \brief Destructor */
    virtual ~ALU3dGridFactory ();

    /** \brief insert a vertex into the coarse grid
     *
     *  \param[in]  pos  position of the vertex
     */
    virtual void insertVertex ( const VertexType &pos );

    /** \brief insert an element into the coarse grid
     *
     *  \note The order of the vertices must coincide with the vertex order in
     *        the corresponding DUNE reference element.
     *
     *  \param[in]  geometry  GeometryType of the new element
     *  \param[in]  vertices  vertices of the new element
     */
    virtual void
    insertElement ( const GeometryType &geometry,
                    const std::vector< unsigned int > &vertices );

    /** \brief insert a boundary element into the coarse grid
     *
     *  \note The order of the vertices must coincide with the vertex order in
     *        the corresponding DUNE reference element.
     *
     *  \param[in]  geometry    GeometryType of the boundary element
     *  \param[in]  vertices    vertices of the boundary element
     *  \param[in]  id          boundary identifier of the boundary element,
     *                          the default value is 0 (invalid boundary id)
     */
    virtual void
    insertBoundary ( const GeometryType &geometry,
                     const std::vector< unsigned int > &faceVertices,
                     const int id );

    void insertBoundary ( const GeometryType &geometry,
                          const DGFEntityKey< unsigned int > &key,
                          const int id );

    /** \brief finalize the grid creation and hand over the grid
     *
     *  The called takes responsibility for deleing the grid.
     */
    GridType *createGrid ();

    GridType *createGrid ( const bool addMissingBoundaries );

  private:
    void assertGeometryType( const GeometryType &geometry );
    static std::string temporaryFileName ();
    void correctElementOrientation ();
    void recreateBoundaryIds ( const int defaultId = 1 );
  };


  template< template< int, int > class ALUGrid >
  inline void ALU3dGridFactory< ALUGrid >
  :: assertGeometryType( const GeometryType &geometry )
  {
    if( elementType == tetra )
    {
      if( !geometry.isSimplex() )
        DUNE_THROW( GridError, "Only simplex geometries can be inserted into "
                    "ALUSimplexGrid< 3, 3 >." );
    }
    else
    {
      if( !geometry.isCube() )
        DUNE_THROW( GridError, "Only cube geometries can be inserted into "
                    "ALUCubeGrid< 3, 3 >." );
    }
  }


  /** \brief Specialization of the generic GridFactory for ALUSimplexGrid<3,3>
   */
  template<>
  class GridFactory< ALUSimplexGrid< 3, 3 > >
    : public ALU3dGridFactory< ALUSimplexGrid >
  {
    typedef GridFactory< ALUSimplexGrid< 3, 3 > > ThisType;
    typedef ALU3dGridFactory< ALUSimplexGrid > BaseType;

  public:
    typedef ALUSimplexGrid< 3, 3 > GridType;

    typedef MPIHelper :: MPICommunicator MPICommunicatorType;

  public:
    /** \brief Default constructor */
    explicit GridFactory ( const MPICommunicatorType &communicator
                             = MPIHelper :: getCommunicator() )
      : BaseType( communicator )
    {}

    /** \brief constructor taking filename */
    GridFactory ( const std::string filename,
                  const MPICommunicatorType &communicator
                    = MPIHelper :: getCommunicator() )
      : BaseType( filename, communicator )
    {}
  };



  /** \brief Specialization of the generic GridFactory for ALUCubeGrid<3,3>
   */
  template<>
  class GridFactory< ALUCubeGrid< 3, 3 > >
    : public ALU3dGridFactory< ALUCubeGrid >
  {
    typedef GridFactory< ALUCubeGrid< 3, 3 > > ThisType;
    typedef ALU3dGridFactory< ALUCubeGrid > BaseType;

  public:
    typedef ALUCubeGrid< 3, 3 > GridType;

    typedef MPIHelper :: MPICommunicator MPICommunicatorType;

  public:
    /** \brief Default constructor */
    explicit GridFactory ( const MPICommunicatorType &communicator
                             = MPIHelper :: getCommunicator() )
      : BaseType( communicator )
    {}

    /** \brief constructor taking filename */
    GridFactory ( const std::string filename,
                  const MPICommunicatorType &communicator
                    = MPIHelper :: getCommunicator() )
      : BaseType( filename, communicator )
    {}
  };

}

// This include is nasty, but we cannot incorporate 'alu3dgridfactory.cc' into
// the lib before HAVE_MPI behaves predictable
#include "alu3dgridfactory.cc"
#endif
