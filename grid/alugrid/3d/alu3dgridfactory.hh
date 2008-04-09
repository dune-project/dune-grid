// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#if defined ENABLE_ALUGRID && !defined DUNE_ALUGRID_FACTORY_HH
#define DUNE_ALUGRID_FACTORY_HH

#include <dune/common/fixedarray.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/alugrid.hh>

namespace Dune
{

  /** \brief Specialization of the generic GridFactory for ALUSimplexGrid<3,3>
   */
  template<>
  class GridFactory< ALUSimplexGrid< 3, 3 > >
    : public GridFactoryInterface< ALUSimplexGrid< 3, 3 > >
  {
    typedef GridFactory< ALUSimplexGrid< 3, 3 > > ThisType;
    typedef GridFactoryInterface< ALUSimplexGrid< 3, 3 > > BaseType;

  public:
    typedef ALUSimplexGrid< 3, 3 > GridType;

    typedef MPIHelper :: MPICommunicator MPICommunicatorType;

  private:
    typedef GridType :: ctype ctype;

    enum { dimension = GridType :: dimension };
    enum { dimensionworld = GridType :: dimensionworld };

    enum { numCorners = EntityCount< tetra > :: numVertices };
    enum { numFaceCorners = EntityCount< tetra > :: numVerticesPerFace };

    typedef ElementTopologyMapping< tetra > ElementTopologyMappingType;
    typedef FaceTopologyMapping< tetra > FaceTopologyMappingType;

    typedef FieldVector< ctype, dimensionworld > VertexType;
    typedef array< unsigned int, numCorners > ElementType;
    typedef array< unsigned int, numFaceCorners > FaceType;

  private:
    MPICommunicatorType communicator_;
#if ALU3DGRID_PARALLEL
    int rank_;
#endif
    std :: vector< VertexType > vertices_;
    std :: vector< ElementType > elements_;
    std :: vector< std :: pair< FaceType, int > > boundaryIds_;

  public:
    /** \brief Default constructor */
    explicit GridFactory ( const MPICommunicatorType &communicator
                             = MPIHelper :: getCommunicator() );

    /** \brief Destructor */
    virtual ~GridFactory ();

    /** \brief Insert a vertex into the coarse grid
     *
     *  \param[in]  pos  position of the vertex
     */
    virtual void insertVertex ( const VertexType &pos );

    /** \brief Insert an element into the coarse grid
     *
     *  \param[in]  type      GeometryType of the new element
     *  \param[in]  vertices  vertices of the new element, using DUNE numbering
     */
    virtual void
    insertElement ( const GeometryType &type,
                    const std :: vector< unsigned int > &vertices );

    virtual void
    insertBoundaryId ( const GeometryType &faceType,
                       const std :: vector< unsigned int > &faceVertices,
                       int id );

    /** \brief Finalize the grid creation and hand over the grid
     *
     *  The called takes responsibility for deleing the grid.
     */
    virtual GridType *createGrid ();
  };



  /** \brief Specialization of the generic GridFactory for ALUCubeGrid<3,3>
   */
  template<>
  class GridFactory< ALUCubeGrid< 3, 3 > >
    : public GridFactoryInterface< ALUCubeGrid< 3, 3 > >
  {
    typedef GridFactory< ALUCubeGrid< 3, 3 > > ThisType;
    typedef GridFactoryInterface< ALUCubeGrid< 3, 3 > > BaseType;

  public:
    typedef ALUCubeGrid< 3, 3 > GridType;

    typedef MPIHelper :: MPICommunicator MPICommunicatorType;

  private:
    typedef GridType :: ctype ctype;

    enum { dimension = GridType :: dimension };
    enum { dimensionworld = GridType :: dimensionworld };

    enum { numCorners = EntityCount< hexa > :: numVertices };
    enum { numFaceCorners = EntityCount< hexa > :: numVerticesPerFace };

    typedef ElementTopologyMapping< hexa > ElementTopologyMappingType;
    typedef FaceTopologyMapping< hexa > FaceTopologyMappingType;

    typedef FieldVector< ctype, dimensionworld > VertexType;
    typedef array< unsigned int, numCorners > ElementType;
    typedef array< unsigned int, numFaceCorners > FaceType;

  private:
    MPICommunicatorType communicator_;
#if ALU3DGRID_PARALLEL
    int rank_;
#endif
    std :: vector< VertexType > vertices_;
    std :: vector< ElementType > elements_;
    std :: vector< std :: pair< FaceType, int > > boundaryIds_;

  public:
    /** \brief Default constructor */
    explicit GridFactory ( const MPICommunicatorType &communicator
                             = MPIHelper :: getCommunicator() );

    /** \brief Destructor */
    virtual ~GridFactory ();

    /** \brief Insert a vertex into the coarse grid
     *
     *  \param[in]  pos  position of the vertex
     */
    virtual void insertVertex ( const VertexType &pos );

    /** \brief Insert an element into the coarse grid
     *
     *  \param[in]  type      GeometryType of the new element
     *  \param[in]  vertices  vertices of the new element, using DUNE numbering
     */
    virtual void
    insertElement ( const GeometryType &type,
                    const std :: vector< unsigned int > &vertices );

    virtual void
    insertBoundaryId ( const GeometryType &faceType,
                       const std :: vector< unsigned int > &faceVertices,
                       int id );

    /** \brief Finalize the grid creation and hand over the grid
     *
     *  The called takes responsibility for deleing the grid.
     */
    virtual GridType *createGrid ();
  };

}

// This include is nasty, but we cannot incorporate 'alu3dgridfactory.cc' into
// the lib before HAVE_MPI behaves predictable
#include "alu3dgridfactory.cc"
#endif
