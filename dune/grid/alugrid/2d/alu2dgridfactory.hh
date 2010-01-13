// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRID_FACTORY_HH
#define DUNE_ALU2DGRID_FACTORY_HH

#ifdef ENABLE_ALUGRID

#include <dune/common/array.hh>
#include <dune/common/mpihelper.hh>

#include <dune/grid/common/genericreferenceelements.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/alugrid/2d/grid.hh>

namespace Dune
{

  /** \brief Factory class for 2d ALUGrids */
  template< template< int, int > class ALUGrid >
  class ALU2dGridFactory
    : public GridFactoryInterface< ALUGrid< 2, 2 > >
  {
  public:
    typedef ALUGrid< 2, 2 > Grid;

    //! \brief type of boundary projection class
    typedef DuneBoundaryProjection< 2 >  DuneBoundaryProjectionType;

    //! \brief type of boundary segment
    typedef Dune::BoundarySegment< 2, 2> BoundarySegmentType;

    template< int codim >
    struct Codim
    {
      typedef typename Grid::template Codim< codim >::Entity Entity;
    };

  private:
    typedef Dune::BoundarySegmentWrapper<2, 2> BoundarySegmentWrapperType;


    typedef ALU2dGridFactory< ALUGrid > ThisType;
    typedef GridFactoryInterface< Grid > BaseType;

    typedef typename Grid::ctype ctype;

    static const unsigned int dimension = Grid::dimension;
    static const unsigned int dimensionworld = Grid::dimensionworld;

    static const unsigned int numCorners = 3;
    static const unsigned int numFaces = 3;
    static const unsigned int numFaceCorners = 2;

    typedef FieldVector< ctype, dimensionworld > VertexType;
    typedef std::vector< unsigned int > ElementType;
    typedef array< unsigned int, numFaceCorners > FaceType;

  private:
    struct FaceLess;

    typedef std::vector< VertexType > VertexVector;
    typedef std::vector< ElementType > ElementVector;
    typedef std::vector< std::pair< FaceType, int > > BoundaryIdVector;

    typedef std::map< FaceType, const DuneBoundaryProjectionType* > BoundaryProjectionMap;
    typedef std::vector< const DuneBoundaryProjectionType* > BoundaryProjectionVector;

    const std::string filename_;
    bool removeGeneratedFile_;
    VertexVector vertices_;
    ElementVector elements_;
    BoundaryIdVector boundaryIds_;
    const DuneBoundaryProjectionType* globalProjection_ ;
    BoundaryProjectionMap boundaryProjections_;
    unsigned int numFacesInserted_;

    // copy vertex numbers and store smalled #dimension ones
    void copyAndSort(const std::vector<unsigned int>& vertices, FaceType& faceId) const
    {
      std::vector<unsigned int> tmp( vertices );
      std::sort( tmp.begin(), tmp.end() );

      // copy only the first dimension vertices (enough for key)
      for( size_t i = 0; i < faceId.size(); ++i ) faceId[ i ] = tmp[ i ];
    }

  public:
    /** \brief default constructor */
    explicit ALU2dGridFactory ( bool removeGeneratedFile = true );

    /** \brief constructor taking filename for temporary outfile */
    explicit ALU2dGridFactory ( const std::string &filename );

    /** \brief Destructor */
    virtual ~ALU2dGridFactory ();

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
     *  \param[in]  faceVertices vertices of the boundary element
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

    /** \brief insert a boundary projection into the macro grid
     *
     *  \param[in]  type        geometry type of boundary face
     *  \param[in]  vertices    vertices of the boundary face
     *  \param[in]  projection  boundary projection
     *
     *  \note The grid takes control of the projection object.
     */
    virtual void
    insertBoundaryProjection ( const GeometryType &type,
                               const std::vector< unsigned int > &vertices,
                               const DuneBoundaryProjectionType *projection );

    /** \brief insert a shaped boundary segment into the macro grid
     *
     *  \param[in]  vertices         vertex indices of boundary face
     *  \param[in]  boundarySegment  geometric realization of shaped boundary
     *
     *  \note The grid takes control of the boundary segment.
     */
    virtual void
    insertBoundarySegment ( const std::vector< unsigned int > vertices,
                            const BoundarySegmentType *boundarySegment ) ;

    /** \brief insert a boundary projection object, (a copy is made)
     *
     *  \param[in]  bndProjection instance of an ALUGridBoundaryProjection projecting vertices to a
     */
    virtual void insertBoundaryProjection ( const DuneBoundaryProjectionType& bndProjection );

    virtual unsigned int
    insertionIndex ( const typename Codim< 0 >::Entity &entity ) const
    {
      return Grid::getRealImplementation( entity ).getIndex();
    }
    virtual unsigned int
    insertionIndex ( const typename Codim< dimension >::Entity &entity ) const
    {
      return Grid::getRealImplementation( entity ).getIndex();
    }
    virtual unsigned int
    insertionIndex ( const typename Grid::LeafIntersection &intersection ) const
    {
      return intersection.boundarySegmentIndex();
    }
    virtual bool
    wasInserted ( const typename Grid::LeafIntersection &intersection ) const
    {
      return ( insertionIndex(intersection) < numFacesInserted_ );
    }

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
  struct ALU2dGridFactory< ALUGrid >::FaceLess
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
  inline void ALU2dGridFactory< ALUGrid >::exchange ( T &x, T &y )
  {
    T dummy = x; x = y; y = dummy;
  }


  template< template< int, int > class ALUGrid >
  inline void ALU2dGridFactory< ALUGrid >
  ::assertGeometryType( const GeometryType &geometry )
  {
    if( !geometry.isSimplex() )
      DUNE_THROW( GridError, "Only simplex geometries can be inserted into "
                  "ALUGrid< 2, 2 >." );
  }


  /** \brief Specialization of the generic GridFactory for ALUConformGrid<2,2>
   *  \ingroup GridFactory
   */
  template<>
  class GridFactory< ALUConformGrid<2,2> >
    : public ALU2dGridFactory<ALUConformGrid>
  {
    typedef GridFactory ThisType;
    typedef ALU2dGridFactory<ALUConformGrid> BaseType;

  public:
    typedef ALUConformGrid< 2, 2 > Grid;

  public:
    /** \brief Default constructor */
    explicit GridFactory ( )
      : BaseType( )
    {}

    /** \brief constructor taking filename */
    GridFactory ( const std::string &filename )
      : BaseType( filename )
    {}
  };
  /** \brief Specialization of the generic GridFactory for ALUSimplexGrid<2,2>
   *  \ingroup GridFactory
   */
  template<>
  class GridFactory< ALUSimplexGrid<2,2> >
    : public ALU2dGridFactory<ALUSimplexGrid>
  {
    typedef GridFactory ThisType;
    typedef ALU2dGridFactory<ALUSimplexGrid> BaseType;

  public:
    typedef ALUSimplexGrid< 2, 2 > Grid;

  public:
    /** \brief Default constructor */
    explicit GridFactory ( )
      : BaseType( )
    {}

    /** \brief constructor taking filename */
    GridFactory ( const std::string &filename )
      : BaseType( filename )
    {}
  };

}

// This include is nasty, but we cannot incorporate 'alu2dgridfactory.cc' into
// the lib before HAVE_MPI behaves predictable
#include "alu2dgridfactory.cc"

#endif // #ifdef ENABLE_ALUGRID

#endif
