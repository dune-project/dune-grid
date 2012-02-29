// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALU2DGRIDGRID_HH
#define DUNE_ALU2DGRIDGRID_HH

//- System includes
#include "alu2dinclude.hh"
#include <iostream>
#include <vector>

//- Dune includes
#include <dune/grid/utility/grapedataioformattypes.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/alugrid/common/interfaces.hh>
#include <dune/common/deprecated.hh>
#include <dune/common/static_assert.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/alugrid/common/declaration.hh>
#include <dune/grid/alugrid/common/defaultindexsets.hh>
#include <dune/grid/common/sizecache.hh>
#include <dune/grid/common/defaultgridview.hh>
#include <dune/common/mpihelper.hh>

#include <dune/grid/alugrid/common/intersectioniteratorwrapper.hh>

// bnd projection stuff
#include <dune/grid/common/boundaryprojection.hh>
#include <dune/grid/alugrid/2d/bndprojection.hh>

//- Local includes
#include "indexsets.hh"
#include <dune/grid/alugrid/common/objectfactory.hh>
#include "datahandle.hh"

namespace Dune {

  // Forward declarations
  template<int cd, int dim, class GridImp>
  class ALU2dGridEntity;
  template<int cd, PartitionIteratorType pitype, class GridImp >
  class ALU2dGridLevelIterator;
  template<int cd, class GridImp >
  class ALU2dGridEntityPointer;
  template<int cd, class GridImp >
  class ALU2dGridEntitySeed;
  template<int mydim, int coorddim, class GridImp>
  class ALU2dGridMakeableGeometry;
  template<int mydim, int cdim, class GridImp>
  class ALU2dGridGeometry;
  template<class GridImp>
  class ALU2dGridHierarchicIterator;
  template<class GridImp>
  class ALU2dGridIntersectionBase;
  template<class GridImp>
  class ALU2dGridLevelIntersectionIterator;
  template<class GridImp>
  class ALU2dGridLeafIntersectionIterator;
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class ALU2dGridLeafIterator;
  template <int mydim, int coorddim, class GridImp>
  class ALU2dGridMakeableEntity;
  template <class GridImp>
  class ALU2dGridFaceGeometryInfo;
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  class ALU2dGridLocalIdSet;
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  class ALU2dGridHierarchicIndexSet;
  template <class EntityImp>
  class ALUMemoryProvider;
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  class ALU2dGrid;
  template <class GridImp, class GeometryImp, int nChild>
  class ALULocalGeometryStorage;

  class ALU2dObjectStream;

  // Internal Forward Declarations
  // -----------------------------

  template < int dimw, class Comm >
  struct ALUGridBaseGrid< 2, dimw, cube, Comm >
  {
    typedef ALU2dGrid< 2, dimw, ALU2DSPACE quadrilateral >  BaseGrid ;
  };

  template < int dimw, class Comm >
  struct ALUGridBaseGrid< 2, dimw, simplex, Comm >
  {
    typedef ALU2dGrid< 2, dimw, ALU2DSPACE triangle >  BaseGrid ;
  };


  //
  // --ALU2dGrid
  // --Grid
  //
  //**********************************************************************
  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  struct ALU2dGridFamily
  {
    typedef ALU2dGrid< dim, dimworld, eltype > GridImp;

    //! Type of the global id set
    typedef ALU2dGridLocalIdSet<dim,dimworld,eltype> GlobalIdSetImp;

    //! Type of the local id set
    typedef ALU2dGridLocalIdSet<dim,dimworld,eltype> LocalIdSetImp;

    typedef int GlobalIdType;
    typedef int LocalIdType;

    struct Traits
    {
      typedef GridImp Grid;

      typedef Dune :: Intersection< const GridImp, LeafIntersectionWrapper > LeafIntersection;
      typedef Dune :: Intersection< const GridImp, LevelIntersectionWrapper > LevelIntersection;

      typedef Dune::IntersectionIterator<const GridImp, LeafIntersectionIteratorWrapper, LeafIntersectionWrapper > IntersectionIterator;
      typedef Dune::IntersectionIterator<const GridImp, LeafIntersectionIteratorWrapper, LeafIntersectionWrapper > LeafIntersectionIterator;
      typedef Dune::IntersectionIterator<const GridImp, LevelIntersectionIteratorWrapper, LevelIntersectionWrapper > LevelIntersectionIterator;

      typedef Dune::EntityIterator< 0, const GridImp, ALU2dGridHierarchicIterator< const GridImp > > HierarchicIterator;

      typedef DuneBoundaryProjection< dimworld > DuneBoundaryProjectionType;
      typedef std::vector< const DuneBoundaryProjectionType *> DuneBoundaryProjectionVector;

      template <int cd>
      struct Codim
      {
        // IMPORTANT: Codim<codim>::Geometry == Geometry<dim-codim,dimw>
        typedef ALU2dGridGeometry< dim-cd, dimworld, const GridImp > GeometryImpl;
        typedef ALU2dGridGeometry< dim-cd, dim, const GridImp > LocalGeometryImpl;
        typedef Dune::Geometry< dim-cd, dimworld, const GridImp, ALU2dGridGeometry > Geometry;
        typedef Dune::Geometry< dim-cd, dim, const GridImp, ALU2dGridGeometry > LocalGeometry;

        // we could - if needed - introduce an other struct for dimglobal of Geometry
        typedef Dune::Entity<cd, dim, const GridImp, ALU2dGridEntity> Entity;

        typedef ALU2dGridEntityPointer< cd, const GridImp > EntityPointerImpl;
        typedef Dune::EntityPointer< const GridImp, EntityPointerImpl > EntityPointer;

        // minimal information to generate entities
        typedef ALU2dGridEntitySeed< cd , const GridImp > EntitySeed ;

        template <PartitionIteratorType pitype>
        struct Partition
        {
          typedef Dune::EntityIterator< cd, const GridImp, ALU2dGridLevelIterator< cd, pitype, const GridImp > > LevelIterator;
          typedef Dune::EntityIterator< cd, const GridImp, ALU2dGridLeafIterator< cd, pitype, const GridImp > > LeafIterator;
        };

        typedef typename Partition< All_Partition >::LevelIterator LevelIterator;
        typedef typename Partition< All_Partition >::LeafIterator LeafIterator;

      };

      template <PartitionIteratorType pitype>
      struct Partition
      {
        typedef Dune::GridView<DefaultLevelGridViewTraits<const GridImp,pitype> >
        LevelGridView;
        typedef Dune::GridView<DefaultLeafGridViewTraits<const GridImp,pitype> >
        LeafGridView;
      };

      //! Type of the level index set
      typedef DefaultIndexSet< GridImp, typename Codim<0>::LevelIterator > LevelIndexSetImp;
      //! Type of the leaf index set
      typedef DefaultIndexSet< GridImp, typename Codim<0>::LeafIterator > LeafIndexSetImp;

      typedef IndexSet<GridImp,LevelIndexSetImp> LevelIndexSet;
      typedef LeafIndexSetImp LeafIndexSet;
      typedef IdSet<GridImp,GlobalIdSetImp,GlobalIdType> GlobalIdSet;
      typedef IdSet<GridImp,LocalIdSetImp,LocalIdType> LocalIdSet;

#if ALU2DGRID_PARALLEL
      typedef Dune :: CollectiveCommunication< MPI_Comm >
      CollectiveCommunication;
#else
      typedef Dune :: CollectiveCommunication< GridImp >
      CollectiveCommunication;
#endif
    };

    //! Type of the level index set implementation
    typedef typename Traits :: LevelIndexSetImp LevelIndexSetImp;

    //! Type of the leaf index set implementation
    typedef typename Traits :: LeafIndexSetImp LeafIndexSetImp;
  }; // end of ALU2dGridFamily


  /**
     \brief [<em> provides \ref Dune::Grid </em>]
     \brief 2D grid, will provide non conform grids
     The ALU2dGrid implements the Dune GridInterface for 2d meshes.
     This grid can be locally adapted and will provide non conform grids.

     @note
     Adaptive grid, written mainly by Bernard Schupp.
     This grid supports non conform grids.

     (see ALUGrid homepage: http://www.mathematik.uni-freiburg.de/IAM/Research/alugrid/)

   */

  template< int dim, int dimworld, ALU2DSPACE ElementType eltype >
  class ALU2dGrid
    : public GridDefaultImplementation< dim, dimworld, alu2d_ctype, ALU2dGridFamily< dim, dimworld, eltype > >,
      public HasObjectStream,
      public HasHierarchicIndexSet
  {
    typedef ALU2dGrid< dim, dimworld, eltype > ThisType;
    typedef GridDefaultImplementation< dim, dimworld, alu2d_ctype, ALU2dGridFamily< dim, dimworld, eltype > > BaseType;

    dune_static_assert( dim == 2, "ALU2dGrid only implemented for grid dim 2." );
    dune_static_assert( dimworld == 2 || dimworld == 3, "ALU2dGrid only implemented for world dim 2 or 3." );

  public:
    static const ALU2DSPACE ElementType elementType = eltype;

    typedef typename ALU2dGridFamily< dim, dimworld, elementType >::Traits Traits;

    // new intersection iterator is a wrapper which get itersectioniteratoimp as pointers
    typedef ALU2dGridLeafIntersectionIterator <const ThisType>  LeafIntersectionIteratorImp;
    typedef ALU2dGridLevelIntersectionIterator<const ThisType>  LevelIntersectionIteratorImp;

    typedef ALUGridObjectFactory< ThisType >  GridObjectFactoryType;

  private:

    typedef typename ALU2dImplTraits<dimworld, elementType >::HmeshType HmeshType ;
    typedef typename ALU2dImplTraits<dimworld, elementType >::HElementType HElementType ;
    typedef typename ALU2dImplTraits<dimworld, elementType >::ElementType ElementType ;

    template< class > friend class DGFBaseFactory;

    template< int, int, class > friend class ALU2dGridEntity;

    friend class ALU2dGridGeometry<0,dimworld,const ThisType>;
    friend class ALU2dGridGeometry<1,dimworld,const ThisType>;
    friend class ALU2dGridGeometry<dim,dimworld,const ThisType>;
    template< class, class, int > friend class ALULocalGeometryStorage;

    friend class ALU2dGridEntityPointer<0,const ThisType>;
    friend class ALU2dGridEntityPointer<1,const ThisType>;
    friend class ALU2dGridEntityPointer<dim,const ThisType>;

    friend class ALU2dGridHierarchicIndexSet<dim,dimworld,elementType>;

    friend class ALU2dGridIntersectionBase < const ThisType > ;
    friend class ALU2dGridLevelIntersectionIterator< const ThisType > ;
    friend class ALU2dGridLeafIntersectionIterator< const ThisType > ;

    //**********************************************************
    // The Interface Methods
    //**********************************************************
  protected:
    typedef MakeableInterfaceObject<typename Traits::template
        Codim<0>::Geometry> GeometryObject;
    friend class ALULocalGeometryStorage<const ThisType, GeometryObject, 4 >;
    friend class ALULocalGeometryStorage<const ThisType, GeometryObject, 2 >;

  public:

    //! dummy object stream
    typedef ALU2dGridObjectStream ObjectStreamType;
    typedef ObjectStreamType InStreamType ;
    typedef ObjectStreamType OutStreamType ;

    //! my Traits class
    typedef ALU2dGridFamily < dim, dimworld, eltype > GridFamily;

    //! Type of the hierarchic index set
    typedef ALU2dGridHierarchicIndexSet<dim,dimworld,elementType> HierarchicIndexSet;

    //! Type of the local id set
    typedef ALU2dGridLocalIdSet<dim,dimworld,elementType> LocalIdSetImp;
    typedef LocalIdSetImp GlobalIdSetImp;

    //! Type of the global id set
    typedef typename Traits :: GlobalIdSet GlobalIdSet;

    //! Type of the local id set
    typedef typename Traits :: LocalIdSet LocalIdSet;


    //! Type of the level index set
    typedef typename GridFamily :: LevelIndexSetImp LevelIndexSetImp;
    //! Type of the leaf index set
    typedef typename GridFamily :: LeafIndexSetImp LeafIndexSetImp;

    //! a standard leaf iterator
    typedef ALU2dGridLeafIterator<0, All_Partition, const ThisType> LeafIteratorImp;
    typedef typename Traits::template Codim<0>::LeafIterator LeafIteratorType;
    typedef typename Traits::template Codim<0>::LeafIterator LeafIterator;

    //! a standard leaf iterator
    typedef ALU2dGridLevelIterator<0, All_Partition, const ThisType> LevelIteratorImp;
    typedef typename Traits::template Codim<0>::LevelIterator LevelIteratorType;
    typedef typename Traits::template Codim<0>::LevelIterator LevelIterator;

    typedef ALU2dGridHierarchicIterator<ThisType> HierarchicIteratorImp;

    typedef typename Traits::CollectiveCommunication CollectiveCommunicationType;


    //! maximal number of levels
    enum {
      //! maximal number of levels is 64
      MAXL = 64
    };

    //! element chunk for refinement
    enum {
      //! \brief normal default number of new elements for new adapt method
      newElementsChunk_ = 100
    };

    //! upper estimate on number of elements that could be created when a new element is created
    enum {
      /** \brief if one element is refined then it
          causes apporximately not more than
          this number of new elements  */
      refineEstimate_ = 40
    };

    //! \brief boundary projection type
    typedef typename Traits :: DuneBoundaryProjectionType DuneBoundaryProjectionType;
    //! \brief boundary projection type
    typedef typename Traits :: DuneBoundaryProjectionVector DuneBoundaryProjectionVector;

#ifdef ALUGRID_VERTEX_PROJECTION
    //! type of ALUGrid Vertex Projection Interface
    typedef ALUGridSpace :: VertexProjection< dimworld > ALUGridVertexProjectionType;
#endif


  protected:

    friend class ALUGridBoundaryProjection< ThisType >;
    friend class ALU2dGridBoundaryProjection< ThisType >;
    // type of ALUGrid boundary projection wrapper
    typedef ALU2dGridBoundaryProjection< ThisType > ALUGridBoundaryProjectionType;

    //! Constructor which reads an ALU2dGrid Macro Triang file
    //! or given GridFile
    //- --constructor
    ALU2dGrid(const std::string macroTriangFilename,
              const int nrOfHangingNodes,
              const DuneBoundaryProjectionType*,
              const DuneBoundaryProjectionVector*,
              std::istream* macroFile = 0);

    // method creating mesh object
    HmeshType* createGrid(const std::string&,
                          const int,
                          std::istream* );

    //! Constructor which constructs an empty ALU2dGrid
    explicit ALU2dGrid( int );

  public:
    //! Desctructor
    ~ALU2dGrid();

    //! Return maximum level defined in this grid. Levels are numbered
    //! 0 ... maxLevel with 0 the coarsest level.
    int maxLevel() const;

    //! --Leveliterator
    //! Iterator to first entity of given codim on level
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LevelIterator
    lbegin (int level) const;

    //! one past the end on this level
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LevelIterator
    lend (int level) const;

    //! Iterator to first entity of given codim on level
    template<int cd>
    typename Traits::template Codim<cd>::
    template Partition<All_Partition>::LevelIterator
    lbegin (int level) const;

    //! one past the end on this level
    template<int cd>
    typename Traits::template Codim<cd>::
    template Partition<All_Partition>::LevelIterator
    lend (int level) const;

    //! Iterator to first entity of codim 0 on level
    LevelIteratorType lbegin (int level) const;

    //! last entity of codim 0 on level
    LevelIteratorType lend (int level) const;

    //! General definiton for a leaf iterator
    template <int codim, PartitionIteratorType pitype>
    typename Traits::template Codim<codim>::template Partition<pitype>::LeafIterator
    leafbegin() const;

    //! General definition for an end iterator on leaf level
    template <int codim, PartitionIteratorType pitype>
    typename Traits::template Codim<codim>::template Partition<pitype>::LeafIterator
    leafend() const;

    //! General definiton for a leaf iterator
    template <int codim>
    typename Traits::template Codim<codim>::LeafIterator
    leafbegin() const;

    //! General definition for an end iterator on leaf level
    template <int codim>
    typename Traits::template Codim<codim>::LeafIterator
    leafend() const;

  private:
    //! Iterator to first entity of codim 0 on leaf level (All_Partition)
    LeafIteratorType leafbegin () const;

    //! one past the end on this leaf level (codim 0 and All_Partition)
    LeafIteratorType leafend () const;

  public:
    //! number of grid entities per level and codim
    int size (int level, int cd) const;

    //! number of leaf entities per codim in this process
    int size (int codim) const;

    //! number of entities per level, codim and geometry type in this process
    int size (int level, GeometryType type) const;

    //! number of leaf entities per codim and geometry type in this process
    int size (GeometryType type) const;

    //! deliver all geometry types used in this grid
    const std::vector<GeometryType>& geomTypes (int codim) const { return geomTypes_[codim]; }

    //****************************************************************
    // index and id sets
    //****************************************************************

    //! get global id set of grid
    const GlobalIdSet & globalIdSet () const;

    //! get global id set of grid
    const LocalIdSet & localIdSet () const;

    //! number of grid entities in the entire grid for given codim
    int hierSetSize (int cd) const;

    //! get hierarchic index set of the grid
    const HierarchicIndexSet & hierarchicIndexSet () const ;

    //! get leaf index set of the grid
    const typename Traits :: LeafIndexSet & leafIndexSet () const;

    //! get level index set of the grid
    const typename Traits :: LevelIndexSet & levelIndexSet (int level) const;


    //**********************************************************
    // End of Interface Methods
    //**********************************************************

    // return reference to org ALU2dGrid
    // private method, but otherwise we have to friend class all possible
    // types of LevelIterator ==> later
    HmeshType & myGrid();
    HmeshType & myGrid() const;

    //! refine grid refCount times
    void globalRefine ( int refCount );

    template< class GridImp, class DataHandle >
    void globalRefine ( int refCount, AdaptDataHandleInterface< GridImp, DataHandle > &hamdle );

    //! returns if a least one entity was marked for coarsening
    bool preAdapt ( );

    //! clear all entity new markers
    void postAdapt ( );

    /**! refine all positive marked leaf entities,
       return true if a least one entity was refined
     */
    bool adapt ();

    template< class GridImp, class DataHandle >
    bool adapt ( AdaptDataHandleInterface< GridImp, DataHandle > &handle );

    // refine grid
    bool refineGrid();

    //! @copydoc Dune::Grid::getMark
    int getMark(const typename Traits::template Codim<0>::Entity & e) const;

    //! @copydoc Dune::Grid::mark
    bool mark( int refCount , const typename Traits::template Codim<0>::Entity & e);

    //! return dummy communication
    const CollectiveCommunicationType & comm() const;
  private:
    CollectiveCommunicationType comm_;

    void updateStatus();

    void calcMaxlevel();

    void calcExtras();

    // clear refinement marker of element and all children
    void hierarchicClear( HElementType *el );
  public:
    typedef MakeableInterfaceObject<typename Traits::template Codim<0>::Entity> EntityObject;
    typedef MakeableInterfaceObject<typename Traits::template Codim<1>::Entity> FaceObject;
    typedef MakeableInterfaceObject<typename Traits::template Codim<2>::Entity> VertexObject;

  protected:
    // create GeomTypes
    void makeGeomTypes ();

    friend class Conversion<ALU2dGrid<dim, dimworld,eltype>, HasObjectStream>;
    friend class Conversion<const ALU2dGrid<dim, dimworld,eltype>, HasObjectStream>;

    friend class Conversion<ALU2dGrid<dim, dimworld,eltype>, HasHierarchicIndexSet>;
    friend class Conversion<const ALU2dGrid<dim, dimworld,eltype>, HasHierarchicIndexSet>;

  private:
    //! Copy constructor should not be used
    ALU2dGrid( const ThisType & g );

    //! assignment operator should not be used
    ThisType & operator = (const ThisType & g);

  protected:
    // check macro file and return const char * to filename
    const char * checkMacroGridFile(const std::string & filename);

    //! the real grid
    mutable HmeshType* mygrid_;

    // return reference to grid
    HmeshType& mesh() const {
      assert(mygrid_);
      return *mygrid_;
    }

#ifdef USE_SMP_PARALLEL
    std::vector< GridObjectFactoryType > factoryVec_;
#else
    GridObjectFactoryType factory_;
#endif

    //! the hierarchic index set
    HierarchicIndexSet hIndexSet_;

    //! out global id set
    LocalIdSetImp localIdSet_;

    //! the level index set ( default type )
    mutable std::vector < LevelIndexSetImp * > levelIndexVec_;

    // at the moment the number of different geom types is 1
    enum { numberOfGeomTypes = 1 };
    std::vector< std::vector<GeometryType> > geomTypes_;

    //! the leaf index set
    mutable LeafIndexSetImp * leafIndexSet_;

    int maxLevel_;
    int refineMarked_ , coarsenMarked_;
    const int nrOfHangingNodes_;

    //! the type of our size cache
    typedef SizeCache<ThisType> SizeCacheType;
    SizeCacheType * sizeCache_;

    // flag to make sure postAdapt is called after adapt
    bool lockPostAdapt_;

    // pointer to Dune boundary projection
    const DuneBoundaryProjectionType* bndPrj_;

    // pointer to Dune boundary projection
    const DuneBoundaryProjectionVector* bndVec_;

    // boundary projection for vertices
    ALUGridBoundaryProjectionType* vertexProjection_ ;

    const DuneBoundaryProjectionType *globalProjection () const
    {
      return bndPrj_;
    }

    //! return boudanry projection for given segment Id
    const DuneBoundaryProjectionType* boundaryProjection(const int segmentIndex) const
    {
      if( bndPrj_ )
      {
        return bndPrj_;
      }
      else
      {
        // note pointer can be zero (identity mapping)
        assert( bndVec_ );
        assert( segmentIndex < (int) bndVec_->size() );
        return (*bndVec_)[ segmentIndex ];
      }
    }

  public:
    //! return number of macro boundary segments
    size_t numBoundarySegments () const
    {
#ifdef ALUGRID_VERTEX_PROJECTION
      return myGrid().numMacroBndSegments();
#else
      derr << "Method available in any version of ALUGrid > 1.14 \n";
      return 0;
#endif
    }

    //! return true if boudanry projection is set
    bool hasBoundaryProjection() const
    {
      return (vertexProjection_ != 0);
    }

    using BaseType :: getRealImplementation ;

  public:
    template< class IntersectionType >
    const typename BaseType
    :: template ReturnImplementationType< IntersectionType>
    :: ImplementationType &
    getRealIntersection ( const IntersectionType &intersection ) const
    {
      return this->getRealImplementation( intersection );
    }

    const GridObjectFactoryType& factory() const {
#ifdef USE_SMP_PARALLEL
      assert( (int) factoryVec_.size() > GridObjectFactoryType :: threadNumber() );
      return factoryVec_[ GridObjectFactoryType :: threadNumber() ];
#else
      return factory_;
#endif
    }

  protected:
    // max level of grid
    int maxlevel_;
    friend class IntersectionIteratorWrapper < const ThisType, LeafIntersectionIteratorImp > ;
    friend class IntersectionIteratorWrapper < const ThisType, LevelIntersectionIteratorImp > ;
    friend class LeafIntersectionIteratorWrapper < const ThisType > ;
    friend class LevelIntersectionIteratorWrapper< const ThisType > ;

    mutable ALU2dGridMarkerVector marker_[MAXL];
  public:
    typedef ALU2dGridLeafMarkerVector ALU2dGridLeafMarkerVectorType;
  private:
    // always update this marker!!!
    mutable ALU2dGridLeafMarkerVectorType leafMarker_;

  public:
    template < class EntitySeed >
    typename Traits :: template Codim< EntitySeed :: codimension > :: EntityPointer
    entityPointer( const EntitySeed& seed ) const
    {
      enum { codim = EntitySeed :: codimension };
      typedef typename Traits :: template Codim< codim > :: EntityPointer EntityPointer;
      typedef ALU2dGridEntityPointer < codim, const ThisType > ALUPointer ;
      return ALUPointer( factory(), seed ) ;
    }

    //! return reference to vector telling on which element a face is
    //! visted for this level
    ALU2dGridMarkerVector & getMarkerVector(int level) const
    {
      assert( level >= 0);
      assert( level <= MAXL);
      return marker_[level];
    }

    //! return reference to vector determing on which element a face is
    //! visited
    ALU2dGridLeafMarkerVectorType & getLeafMarker() const
    {
      return leafMarker_;
    }

    /** \brief write Grid to file in specified FileFormatType
     */
    template <GrapeIOFileFormatType ftype>
    bool writeGrid( const std::string filename, alu2d_ctype time ) const ;

    bool writeGrid_Xdr( const std::string filename, alu2d_ctype time ) const ;
    bool writeGrid_Ascii( const std::string filename, alu2d_ctype time ) const ;

    /** \brief read Grid from file filename and store time of mesh in time
     */
    template <GrapeIOFileFormatType ftype>
    bool readGrid( const std::string filename, alu2d_ctype & time );

  protected:
    //! return true if grid allows hanging nodes on leaf level
    //! i.e. returns true for ALUSimplexGrid and returns false for ALUConformGrid
    bool nonConform () const
    {
      return (nrOfHangingNodes_ > 0);
    }

#if ALU2DGRID_PARALLEL
    typedef RankManager<ThisType> RankManagerType;
    RankManagerType rankManager_;
  public:
    const RankManagerType& rankManager() const
    {
      return rankManager_;
    }
#endif

  public:
    /** \brief @copydoc Dune::Grid::communicate */
    template<class DataHandleImp,class DataTypeImp>
    void communicate (CommDataHandleIF<DataHandleImp,DataTypeImp> & data,
                      InterfaceType iftype, CommunicationDirection dir, int level) const;

    /** \brief Communicate information on distributed entities on the leaf grid.
       Template parameter is a model of Dune::CommDataHandleIF.
     */
    template<class DataHandleImp,class DataTypeImp>
    void communicate (CommDataHandleIF<DataHandleImp,DataTypeImp> & data,
                      InterfaceType iftype, CommunicationDirection dir) const;

    int ghostSize ( int level, int codim ) const
    {
      return ghostSize( codim );
    }

    int ghostSize ( int codim ) const
    {
#if ALU2DGRID_PARALLEL
      return 1;
#else
      return 0;
#endif
    }

    /** \brief @copydoc Dune::Grid::loadBalance */
    bool loadBalance() ;

    /** \brief @copydoc Dune::Grid::loadBalance */
    template<class DataHandle>
    bool loadBalance(DataHandle& data) ;

    void checkManager() {
#if ALU2DGRID_PARALLEL
      rankManager_.notifyMarking () ;
#endif
    }

  }; // end class ALU2dGrid

  namespace Capabilities
  {
    template<int dim, int dimw, ALU2DSPACE ElementType eltype, int cdim>
    struct hasEntity<ALU2dGrid<dim,dimw,eltype>, cdim >
    {
      static const bool v = true;
    };

    template<int dim, int dimw, ALU2DSPACE ElementType eltype>
    struct isLevelwiseConforming< ALU2dGrid<dim,dimw,eltype> >
    {
      static const bool v = false;
    };

  } // end namespace Capabilities

} // end namespace Dune

#include "entity.hh"
#include "geometry.hh"
#include <dune/grid/alugrid/2d/intersection.hh>
#include <dune/grid/alugrid/2d/iterator.hh>

#include "grid_imp.cc"

#endif
