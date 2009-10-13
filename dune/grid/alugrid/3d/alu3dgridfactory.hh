// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU3DGRID_FACTORY_HH
#define DUNE_ALU3DGRID_FACTORY_HH

#ifdef ENABLE_ALUGRID

#include <dune/common/array.hh>
#include <dune/common/mpihelper.hh>

#include <dune/grid/common/genericreferenceelements.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/boundaryprojection.hh>

#include <dune/grid/alugrid/3d/alugrid.hh>

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
    typedef ALUGrid< 3, 3 > Grid;

    typedef MPIHelper::MPICommunicator MPICommunicatorType;

    //! \brief type of boundary projection class
    typedef DuneBoundaryProjection< 3 >  DuneBoundaryProjectionType;

  private:
    typedef typename Grid::ctype ctype;

    static const ALU3dGridElementType elementType = Grid::elementType;

    dune_static_assert( (elementType == tetra || elementType == hexa),
                        "ALU3dGridFactory supports only grids containing "
                        "tetrahedrons or hexahedrons exclusively." );

    static const unsigned int dimension = Grid::dimension;
    static const unsigned int dimensionworld = Grid::dimensionworld;

    static const unsigned int numCorners = EntityCount< elementType >::numVertices;
    static const unsigned int numFaces = EntityCount< elementType >::numFaces;
    static const unsigned int numFaceCorners = EntityCount< elementType >::numVerticesPerFace;

    typedef ElementTopologyMapping< elementType > ElementTopologyMappingType;
    typedef FaceTopologyMapping< elementType > FaceTopologyMappingType;

    typedef FieldVector< ctype, dimensionworld > VertexType;
    typedef std::vector< unsigned int > ElementType;
    typedef array< unsigned int, numFaceCorners > FaceType;

  private:
    struct FaceLess;

    typedef std::vector< VertexType > VertexVector;
    typedef std::vector< ElementType > ElementVector;
    typedef std::vector< std::pair< FaceType, int > > BoundaryIdVector;

    const std::string filename_;
    bool removeGeneratedFile_;
    MPICommunicatorType communicator_;
#if ALU3DGRID_PARALLEL
    int rank_;
#endif
    VertexVector vertices_;
    ElementVector elements_;
    BoundaryIdVector boundaryIds_;
    const DuneBoundaryProjectionType* bndPrjct_ ;

  public:
    /** \brief default constructor */
    explicit ALU3dGridFactory ( const MPICommunicatorType &communicator
                                  = MPIHelper::getCommunicator(),
                                bool removeGeneratedFile = true );

    /** \brief constructor taking filename for temporary outfile */
    explicit ALU3dGridFactory ( const std::string &filename,
                                const MPICommunicatorType &communicator
                                  = MPIHelper::getCommunicator() );

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

    /** \brief mark a face as boundary (and assign a boundary id)
     *
     *  \param[in]  element  index of the element, the face belongs to
     *  \param[in]  face     local number of the face within the element
     *  \param[in]  id       boundary id to assign to the face
     */
    virtual void insertBoundary ( const int element, const int face, const int id );

    /** \brief insert a boundary projection object, (a copy is made)
     *
     *  \param[in]  bndProjection instance of an ALUGridBoundaryProjection projecting vertices to a curved
     */
    virtual void insertBoundaryProjection ( const DuneBoundaryProjectionType& bndProjection );

    /** \brief finalize the grid creation and hand over the grid
     *
     *  The called takes responsibility for deleing the grid.
     */
    Grid *createGrid ();

    Grid *createGrid ( const bool addMissingBoundaries );

  private:
    template< class T >
    static void exchange ( T &x, T &y );

    void assertGeometryType( const GeometryType &geometry );
    static std::string temporaryFileName ();
    static void generateFace ( const ElementType &element, const int f, FaceType &face );
    void correctElementOrientation ();
    void recreateBoundaryIds ( const int defaultId = 1 );
  };


  template< template< int, int > class ALUGrid >
  struct ALU3dGridFactory< ALUGrid >::FaceLess
    : public std::binary_function< FaceType, FaceType, bool >
  {
    bool operator() ( const FaceType &a, const FaceType &b ) const
    {
      for( unsigned int i = 0; i < numFaceCorners; ++i )
      {
        if( a[ i ] != b[ i ] )
          return (a[ i ] < b[ i ]);
      }
      return false;
    }
  };


  template< template< int, int > class ALUGrid >
  template< class T >
  inline void ALU3dGridFactory< ALUGrid >::exchange ( T &x, T &y )
  {
    T dummy = x; x = y; y = dummy;
  }


  template< template< int, int > class ALUGrid >
  inline void ALU3dGridFactory< ALUGrid >
  ::assertGeometryType( const GeometryType &geometry )
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
   *  \ingroup GridFactory
   */
  template<>
  class GridFactory< ALUSimplexGrid< 3, 3 > >
    : public ALU3dGridFactory< ALUSimplexGrid >
  {
    typedef GridFactory< ALUSimplexGrid< 3, 3 > > ThisType;
    typedef ALU3dGridFactory< ALUSimplexGrid > BaseType;

  public:
    typedef ALUSimplexGrid< 3, 3 > Grid;

    typedef MPIHelper::MPICommunicator MPICommunicatorType;

  public:
    /** \brief Default constructor */
    explicit GridFactory ( const MPICommunicatorType &communicator
                             = MPIHelper::getCommunicator() )
      : BaseType( communicator )
    {}

    /** \brief constructor taking filename */
    GridFactory ( const std::string &filename,
                  const MPICommunicatorType &communicator
                    = MPIHelper::getCommunicator() )
      : BaseType( filename, communicator )
    {}
  };



  /** \brief Specialization of the generic GridFactory for ALUCubeGrid<3,3>
   *  \ingroup GridFactory
   */
  template<>
  class GridFactory< ALUCubeGrid< 3, 3 > >
    : public ALU3dGridFactory< ALUCubeGrid >
  {
    typedef GridFactory< ALUCubeGrid< 3, 3 > > ThisType;
    typedef ALU3dGridFactory< ALUCubeGrid > BaseType;

  public:
    typedef ALUCubeGrid< 3, 3 > Grid;

    typedef MPIHelper::MPICommunicator MPICommunicatorType;

  public:
    /** \brief Default constructor */
    explicit GridFactory ( const MPICommunicatorType &communicator
                             = MPIHelper::getCommunicator() )
      : BaseType( communicator )
    {}

    /** \brief constructor taking filename */
    GridFactory ( const std::string &filename,
                  const MPICommunicatorType &communicator
                    = MPIHelper::getCommunicator() )
      : BaseType( filename, communicator )
    {}
  };

}

// This include is nasty, but we cannot incorporate 'alu3dgridfactory.cc' into
// the lib before HAVE_MPI behaves predictable
#include "alu3dgridfactory.cc"

#endif // #ifdef ENABLE_ALUGRID

#endif
