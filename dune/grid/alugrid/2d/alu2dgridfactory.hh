// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_ALU2DGRID_FACTORY_HH
#define DUNE_ALU2DGRID_FACTORY_HH

#ifdef ENABLE_ALUGRID

#include <dune/common/container/array.hh>
#include <dune/common/mpihelper.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/gridfactory.hh>

#include <dune/grid/alugrid/common/persistentcontainer.hh>
#include <dune/grid/alugrid/common/transformation.hh>
#include <dune/grid/alugrid/2d/grid.hh>

namespace Dune
{

  /** \brief Factory class for 2d ALUGrids */
  template< class GridImp >
  class ALU2dGridFactory
    : public GridFactoryInterface< GridImp >
  {
  public:
    typedef GridImp Grid;

    static const int dimension = Grid::dimension;
    static const int dimensionworld = Grid::dimensionworld;

    //! \brief type of boundary projection class
    typedef DuneBoundaryProjection< dimensionworld >  DuneBoundaryProjectionType;

    template< int codim >
    struct Codim
    {
      typedef typename Grid::template Codim< codim >::Entity Entity;
    };

  private:
    typedef Dune::BoundarySegmentWrapper<2, dimensionworld> BoundarySegmentWrapperType;

    typedef ALU2dGridFactory< Grid > ThisType;
    typedef GridFactoryInterface< Grid > BaseType;

    typedef typename Grid::ctype ctype;

    static const ALU2DSPACE ElementType elementType = Grid::elementType;
    static const unsigned int numFaceCorners = 2;

    typedef FieldVector< ctype, dimensionworld > VertexType;
    typedef std::vector< unsigned int > ElementType;
    typedef array< unsigned int, numFaceCorners > FaceType;

    static const int periodicBndId = ALU2dImplTraits< dimensionworld, elementType >::HBndElType::general_periodic;

  public:
    typedef ALUGridTransformation< ctype, dimensionworld > Transformation;

    //! type of vector for world coordinates
    typedef typename Transformation::WorldVector WorldVector;
    //! type of matrix from world coordinates to world coordinates
    typedef typename Transformation::WorldMatrix WorldMatrix;

  private:
    struct FaceLess;

    typedef std::vector< VertexType > VertexVector;
    typedef std::vector< ElementType > ElementVector;
    typedef std::vector< std::pair< FaceType, int > > BoundaryIdVector;

    typedef std::map< FaceType, const DuneBoundaryProjectionType* > BoundaryProjectionMap;
    typedef std::vector< const DuneBoundaryProjectionType* > BoundaryProjectionVector;

    typedef std::pair< unsigned int, int > SubEntity;
    typedef std::map< FaceType, SubEntity, FaceLess > FaceMap;
    typedef std::vector< Transformation > FaceTransformationVector;
    typedef std::map< FaceType, unsigned int, FaceLess > PeriodicNeighborMap;

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

    /** \brief insert a boundary segment into the macro grid
     *
     *  \param[in]  vertices         vertex indices of boundary face
     */
    virtual void
    insertBoundarySegment ( const std::vector< unsigned int >& vertices ) ;

    /** \brief insert a shaped boundary segment into the macro grid
     *
     *  \param[in]  vertices         vertex indices of boundary face
     *  \param[in]  boundarySegment  geometric realization of shaped boundary
     */
    virtual void
    insertBoundarySegment ( const std::vector< unsigned int >& vertices,
                            const shared_ptr<BoundarySegment<2,dimensionworld> >& boundarySegment ) ;

    /** \brief insert a boundary projection object, (a copy is made)
     *
     *  \param[in]  bndProjection instance of an ALUGridBoundaryProjection projecting vertices to a
     */
    virtual void insertBoundaryProjection ( const DuneBoundaryProjectionType& bndProjection );

    /** \brief add a face transformation (for periodic identification)
     *
     *  A face transformation is an affine mapping T from world coordinates
     *  to world coordinates. The grid factory then glues two faces f and g
     *  if T( f ) = g or T( g ) = f.
     *
     *  \param[in]  matrix  matrix describing the linear part of T
     *  \param[in]  shift   vector describing T( 0 )
     */
    void insertFaceTransformation ( const WorldMatrix &matrix, const WorldVector &shift );

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
      return intersection.boundary() &&
             ( insertionIndex(intersection) < numFacesInserted_ );
    }

    /** \brief finalize the grid creation and hand over the grid
     *
     *  The called takes responsibility for deleing the grid.
     */
    Grid *createGrid ();

    Grid *createGrid ( const bool addMissingBoundaries, const std::string dgfName = "" );

    Grid *createGrid ( const bool addMissingBoundaries, bool temporary, const std::string dgfName = "" );

    void setTolerance ( const ctype &epsilon ) { epsilon_ = epsilon; }

  protected:
    //! virtual method for creating grid object (virtual to avoid compiling method into library)
    virtual Grid* createGridObj( const bool temporary,
                                 const std::string& filename,
                                 std::istream& inFile,
                                 BoundaryProjectionVector* bndProjections )
    {
      return ( temporary ) ?
             new Grid( filename, inFile, globalProjection_, bndProjections, grdVerbose_ ) :
             new Grid( filename, globalProjection_ , bndProjections, grdVerbose_ );
    }

    /** \brief set factory's verbosity
     *
     *  \param[in]  verbose  verbosity (true/flase)
     */
    void setVerbosity( const bool verbose ) { grdVerbose_ = verbose; }

  private:
    static void generateFace ( const ElementType &element, const int f, FaceType &face );
    void correctElementOrientation ();
    typename FaceMap::const_iterator findPeriodicNeighbor( const FaceMap &faceMap, const FaceType &key ) const;
    void reinsertBoundary ( const FaceMap &faceMap, const typename FaceMap::const_iterator &pos, const int id );
    void recreateBoundaryIds ( const int defaultId = 1 );

    VertexVector vertices_;
    ElementVector elements_;
    BoundaryIdVector boundaryIds_;
    const DuneBoundaryProjectionType* globalProjection_ ;
    BoundaryProjectionMap boundaryProjections_;
    unsigned int numFacesInserted_;
    bool grdVerbose_;
    FaceTransformationVector faceTransformations_;
    PeriodicNeighborMap periodicNeighborMap_;
    ctype epsilon_;
  };


  template< class GridImp >
  struct ALU2dGridFactory< GridImp >::FaceLess
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


  /** \brief Specialization of the generic GridFactory for ALUConformGrid<2,dimw>
   *  \ingroup GridFactory
   */
  template<int dimw>
  class GridFactory< ALUConformGrid<2,dimw> >
    : public ALU2dGridFactory<ALUConformGrid<2, dimw> >
  {
  public:
    typedef ALUConformGrid< 2, dimw > Grid;

  protected:
    typedef GridFactory ThisType;
    typedef ALU2dGridFactory< Grid > BaseType;

  public:
    /** \brief Default constructor */
    explicit GridFactory ( )
      : BaseType( )
    {}

    /** \brief constructor taking filename */
    GridFactory ( const std::string &filename )
      : BaseType( filename )
    {}

    /** \brief constructor setting verbosity level */
    GridFactory ( const bool verbose )
      : BaseType( )
    {
      this->setVerbosity( verbose );
    }
  };
  /** \brief Specialization of the generic GridFactory for ALUSimplexGrid<2,dimw>
   *  \ingroup GridFactory
   */
  template<int dimw>
  class GridFactory< ALUSimplexGrid<2,dimw> >
    : public ALU2dGridFactory<ALUSimplexGrid<2, dimw> >
  {
  public:
    typedef ALUSimplexGrid< 2, dimw > Grid;

  protected:
    typedef GridFactory ThisType;
    typedef ALU2dGridFactory< Grid > BaseType;

  public:
    /** \brief Default constructor */
    explicit GridFactory ( )
      : BaseType( )
    {}

    /** \brief constructor taking filename */
    GridFactory ( const std::string &filename )
      : BaseType( filename )
    {}

    /** \brief constructor setting verbosity flag
     *
     *  \param[in] verbose verbosity (true/false)
     */
    GridFactory ( const bool verbose )
      : BaseType( )
    {
      this->setVerbosity( verbose );
    }
  };

  /** \brief Specialization of the generic GridFactory for ALUCubeGrid<2,dimw>
   *  \ingroup GridFactory
   */
  template<int dimw>
  class GridFactory< ALUCubeGrid<2,dimw> >
    : public ALU2dGridFactory<ALUCubeGrid<2,dimw> >
  {
  public:
    typedef ALUCubeGrid< 2, dimw > Grid;

  protected:
    typedef GridFactory ThisType;
    typedef ALU2dGridFactory< Grid > BaseType;

  public:
    /** \brief Default constructor */
    explicit GridFactory ( )
      : BaseType( )
    {}

    /** \brief constructor taking filename */
    GridFactory ( const std::string &filename )
      : BaseType( filename )
    {}

    /** \brief constructor setting verbosity flag
     *
     *  \param[in] verbose verbosity (true/false)
     */
    GridFactory ( const bool verbose )
      : BaseType( )
    {
      this->setVerbosity( verbose );
    }
  };

  /** \brief Specialization of the generic GridFactory for
   *        ALUGrid<2,dimw, eltype, refinementtype, Comm >
   *  \ingroup GridFactory
   */
  template<int dimw, ALUGridElementType eltype, ALUGridRefinementType refinementtype, class Comm>
  class GridFactory< ALUGrid<2, dimw, eltype, refinementtype, Comm > >
    : public ALU2dGridFactory< ALUGrid<2, dimw, eltype, refinementtype, Comm > >
  {
  public:
    typedef ALUGrid<2, dimw, eltype, refinementtype, Comm > Grid;

  protected:
    typedef GridFactory ThisType;
    typedef ALU2dGridFactory< Grid > BaseType;

  public:
    /** \brief Default constructor */
    explicit GridFactory ( )
      : BaseType( )
    {}

    /** \brief constructor taking filename */
    GridFactory ( const std::string &filename )
      : BaseType( filename )
    {}

    /** \brief constructor setting verbose flag
     *
     *  \param[in] verbose verbosity (true/false)
     */
    GridFactory ( const bool verbose )
      : BaseType( )
    {
      this->setVerbosity( verbose );
    }
  };



  // Inline Implementations
  // ----------------------

  template< class GridImp >
  inline ALU2dGridFactory< GridImp >::ALU2dGridFactory ( bool removeGeneratedFile )
    : globalProjection_ ( 0 ),
      numFacesInserted_ ( 0 ),
      grdVerbose_( true ),
      epsilon_( 1e-8 )
  {}


  template< class GridImp >
  inline ALU2dGridFactory< GridImp >::ALU2dGridFactory ( const std::string &filename )
    : globalProjection_ ( 0 ),
      numFacesInserted_ ( 0 ),
      grdVerbose_( true ),
      epsilon_( 1e-8 )
  {}


  template< class GridImp >
  inline ALU2dGridFactory< GridImp >::~ALU2dGridFactory ()
  {}


  template< class GridImp >
  inline GridImp *ALU2dGridFactory< GridImp >::createGrid ()
  {
    return createGrid( true, true, "" );
  }


  template< class GridImp >
  inline GridImp *ALU2dGridFactory< GridImp >
  ::createGrid ( const bool addMissingBoundaries, const std::string dgfName )
  {
    return createGrid( addMissingBoundaries, true, dgfName );
  }

}

#endif // #ifdef ENABLE_ALUGRID

#endif // #ifndef DUNE_ALU2DGRID_FACTORY_HH
