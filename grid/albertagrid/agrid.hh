// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTAGRID_IMP_HH
#define DUNE_ALBERTAGRID_IMP_HH

#if HAVE_ALBERTA

#include <iostream>
#include <fstream>
#include <dune/common/deprecated.hh>

#include <vector>
#include <assert.h>
#include <algorithm>

/** @file
   @author Robert Kloefkorn
   @brief Provides base classes for AlbertaGrid
 **/

// Dune includes
#include <dune/common/misc.hh>
#include <dune/common/interfaces.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/common/stdstreams.hh>

#if HAVE_MPI
#include <dune/common/mpicollectivecommunication.hh>
#else
#include <dune/common/collectivecommunication.hh>
#endif

#include <dune/common/exceptions.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/defaultindexsets.hh>
#include <dune/grid/common/sizecache.hh>
#include <dune/grid/common/intersectioniteratorwrapper.hh>
#include <dune/grid/common/defaultgridview.hh>

// stack for index management
#include <dune/grid/common/indexstack.hh>

//- Local includes
// some cpp defines and include of alberta.h
#include "albertaheader.hh"

// grape data io
#include <dune/grid/utility/grapedataioformattypes.hh>

// IndexManager defined in indexstack.hh
// 10000 is the size of the finite stack used by IndexStack
typedef Dune::IndexStack<int,100000> IndexManagerType;

//#define CALC_COORD
// some extra functions for handling the Albert Mesh
#include "albertaextra.hh"

#include <dune/grid/albertagrid/exceptions.hh>

// contains a simple memory management for some componds of this grid
#include "agmemory.hh"

namespace Dune
{
  // i.e. double or float
  typedef ALBERTA REAL albertCtype;
}

#include "referencetopo.hh"
#include "indexsets.hh"
#include "geometry.hh"
#include "entity.hh"
#include "entitypointer.hh"
#include "hierarchiciterator.hh"
#include "treeiterator.hh"

namespace Dune
{
  template<int codim, int dim, class GridImp> class AlbertaGridEntity;
  template<int codim, PartitionIteratorType pitype, class GridImp> class AlbertaGridTreeIterator;
  template<int codim, PartitionIteratorType pitype, class GridImp> class AlbertaGridLeafIterator;
  template<int cd, class GridImp> class AlbertaGridEntityPointer;

  template<class GridImp>         class AlbertaGridHierarchicIterator;
  template<class GridImp>         class AlbertaGridIntersectionIterator;
  template<int dim, int dimworld> class AlbertaGrid;
  template<int dim, int dimworld> class AlbertaGridHierarchicIndexSet;

  template <int codim, int dim, class GridImp>
  struct SelectEntityImp
  {
    typedef AlbertaGridEntity<codim,dim,GridImp> EntityImp;
    typedef Dune::Entity<codim, dim, const GridImp, AlbertaGridEntity> Entity;
    typedef MakeableInterfaceObject<Entity> EntityObject;
  };

  //**********************************************************************
  //
  // --AlbertaGridIntersectionIterator
  // --IntersectionIterator
  /*!
     Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
     a neighbor is an entity of codimension 0 which has a common entity of codimension 1
     These neighbors are accessed via a IntersectionIterator. This allows the implementation of
     non-matching meshes. The number of neigbors may be different from the number of faces
     of an element!
   */
  template<class GridImp>
  class AlbertaGridIntersectionIterator
  {
    enum { dim      = GridImp::dimension };
    enum { dimworld = GridImp::dimensionworld };

    friend class AlbertaGridEntity<0,dim,GridImp>;
    typedef AlbertaGridIntersectionIterator<GridImp> ThisType;

  public:
    typedef Dune::Intersection< GridImp, Dune::AlbertaGridIntersectionIterator >
    Intersection;
    typedef ThisType ImplementationType;

    typedef AGMemoryProvider< ThisType > StorageType;
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;

    typedef typename SelectEntityImp<0,dim,GridImp>::EntityImp EntityImp;

    //typedef AlbertaGridMakeableGeometry<dim-1,dimworld,GridImp> LocalGeometryImp;
    typedef AlbertaGridGeometry<dim-1,dimworld,GridImp> LocalGeometryImp;
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;

    //! know your own dimension
    enum { dimension=dim };
    //! know your own dimension of world
    enum { dimensionworld=dimworld };
    //! define type used for coordinates in grid module
    typedef typename GridImp::ctype ctype;

    const Intersection &dereference () const
    {
      return reinterpret_cast< const Intersection & >( *this );
    }

    //! equality
    bool equals (const AlbertaGridIntersectionIterator<GridImp> & i) const;

    //! increment
    void increment();

    //! equality
    bool operator==(const AlbertaGridIntersectionIterator<GridImp>& i) const;

    //! access neighbor
    EntityPointer outside() const;

    //! access element where IntersectionIterator started
    EntityPointer inside() const;

    //! The default Constructor
    AlbertaGridIntersectionIterator(const GridImp & grid,
                                    int level);

    //! The Constructor
    AlbertaGridIntersectionIterator(const GridImp & grid,
                                    int level,
                                    ALBERTA EL_INFO *elInfo,
                                    bool leafIt );
    //! The copy constructor
    AlbertaGridIntersectionIterator(const AlbertaGridIntersectionIterator<GridImp> & org);

    //! assignment operator, implemented because default does not the right thing
    void assign (const AlbertaGridIntersectionIterator<GridImp> & org);

    //! The Destructor
    //~AlbertaGridIntersectionIterator();

    //! return true if intersection is with boundary.
    bool boundary () const;

    //! return true if across the edge an neighbor on this level exists
    bool neighbor () const;

    //! return information about the Boundary
    int boundaryId () const;

    //! return true if intersection is conform.
    bool conforming () const;

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    const LocalGeometry& intersectionSelfLocal () const;
    /*! intersection of codimension 1 of this neighbor with element where iteration started.
       Here returned element is in LOCAL coordinates of neighbor
     */
    const LocalGeometry& intersectionNeighborLocal () const;
    /*! intersection of codimension 1 of this neighbor with element where iteration started.
       Here returned element is in GLOBAL coordinates of the element where iteration started.
     */
    const Geometry& intersectionGlobal () const;

    //! local number of codim 1 entity in self where intersection is contained in
    int numberInSelf () const;
    //! local number of codim 1 entity in neighbor where intersection is contained in
    int numberInNeighbor () const;

    //! twist of the face seen from the inner element
    int twistInSelf() const;

    //! twist of the face seen from the outer element
    int twistInNeighbor() const;

    //! return unit outer normal, this should be dependent on local
    //! coordinates for higher order boundary
    typedef FieldVector<albertCtype, GridImp::dimensionworld> NormalVecType;
    typedef FieldVector<albertCtype, GridImp::dimension-1> LocalCoordType;

    const NormalVecType & unitOuterNormal (const LocalCoordType & local) const;

    //! return outer normal, this should be dependent on local
    //! coordinates for higher order boundary
    const NormalVecType & outerNormal (const LocalCoordType & local) const;

    //! return outer normal, this should be dependent on local
    //! coordinates for higher order boundary
    const NormalVecType & integrationOuterNormal (const LocalCoordType & local) const;

    //! return level of inside entity
    int level () const;

    //**********************************************************
    //  private methods
    //**********************************************************

    // reset IntersectionIterator
    template <class EntityType>
    void first(const EntityType & en, int level );

    // calls EntityPointer done and sets done_ to true
    void done ();

  private:
    // returns true if actual neighbor has same level
    bool neighborHasSameLevel () const;

    //! make Iterator set to begin of actual entitys intersection Iterator
    void makeBegin (const GridImp & grid,
                    int level,
                    ALBERTA EL_INFO * elInfo ) const;

    //! set Iterator to end of actual entitys intersection Iterator
    void makeEnd (const GridImp & grid,int level ) const;

    // put objects on stack
    void freeObjects () const;

    //! setup the virtual neighbor
    void setupVirtEn () const;

    //! calculate normal to current face
    void calcOuterNormal () const;

    // return whether the iterator was called from a LeafIterator entity or
    // LevelIterator entity
    bool leafIt () const { return leafIt_; }
    ////////////////////////////////////////////////
    // private member variables
    ////////////////////////////////////////////////

    //! know the grid were im coming from
    const GridImp& grid_;

    //! the actual level
    mutable int level_;

    //! count on which neighbor we are lookin' at
    mutable int neighborCount_;

    //! implement with virtual element
    //! Most of the information can be generated from the ALBERTA EL_INFO
    //! therefore this element is only created on demand.
    mutable bool builtNeigh_;

    bool leafIt_;

    //! pointer to the EL_INFO struct storing the real element information
    mutable ALBERTA EL_INFO * elInfo_;

    typedef MakeableInterfaceObject<LocalGeometry> LocalGeometryObject;

    // the objects holding the real implementations
    mutable LocalGeometryObject fakeNeighObj_;
    mutable LocalGeometryObject fakeSelfObj_;
    mutable LocalGeometryObject neighGlobObj_;

    //! pointer to element holding the intersectionNeighbourLocal information.
    //! This element is created on demand.
    mutable LocalGeometryImp & fakeNeigh_;
    //! pointer to element holding the intersectionSelfLocal information.
    //! This element is created on demand.
    mutable LocalGeometryImp & fakeSelf_;
    //! pointer to element holding the neighbor_global and neighbor_local
    //! information. This element is created on demand.
    mutable LocalGeometryImp & neighGlob_;

    //! EL_INFO th store the information of the neighbor if needed
    mutable ALBERTA EL_INFO neighElInfo_;

    mutable NormalVecType outNormal_;
    mutable NormalVecType unitNormal_;

    // tmp memory for normal calculation
    mutable FieldVector<albertCtype, dimworld> tmpU_;
    mutable FieldVector<albertCtype, dimworld> tmpV_;

    // twist seen from the neighbor
    mutable int twist_;

    //! is true when iterator finished
    bool done_;
  };


  //**********************************************************************
  //
  // --AlbertaGridTreeIterator
  // --LevelIterator
  // --TreeIterator
  //


  //! --LevelIterator
  //! the same as TreeIterator
  template<int cd, PartitionIteratorType pitype, class GridImp>
  class AlbertaGridLevelIterator
    : public AlbertaGridTreeIterator<cd,pitype,GridImp>
  {
  public:
    typedef typename GridImp::template Codim<cd>::Entity Entity;

    //! Constructor making end iterator
    AlbertaGridLevelIterator(const GridImp & grid, int level, int proc) :
      AlbertaGridTreeIterator<cd,pitype,GridImp> (grid,level,proc)
    {}

    //! Constructor making begin iterator
    AlbertaGridLevelIterator(const GridImp & grid,
                             const AlbertaMarkerVector * vec, int level, int proc) :
      AlbertaGridTreeIterator<cd,pitype,GridImp> (grid,vec,level,proc)
    {}

    //! increment the iterator
    void increment ()
    {
      AlbertaGridTreeIterator<cd,pitype,GridImp>::increment();
    }
  };

  //**********************************************************************
  //
  //  AlbertaGridLeafIterator
  //  --LeafIterator
  //
  //**********************************************************************
  //! LeafIterator which is just a hull for the LevelIterator
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class AlbertaGridLeafIterator
    : public AlbertaGridTreeIterator<codim, pitype, GridImp>
  {
  public:
    typedef typename GridImp::template Codim<codim>::Entity Entity;

    //! Constructor making end iterator
    AlbertaGridLeafIterator(const GridImp & grid, int level, int proc) :
      AlbertaGridTreeIterator<codim,pitype,GridImp> (grid,level,proc,true)
    {}

    //! Constructor making begin iterator
    AlbertaGridLeafIterator(const GridImp & grid,
                            const AlbertaMarkerVector * vec, int level, int proc) :
      AlbertaGridTreeIterator<codim, pitype, GridImp> (grid,vec,level,proc,true)
    {}

    //! increment the iterator
    void increment ()
    {
      AlbertaGridTreeIterator<codim, pitype, GridImp>::increment();
    }
  };

  //**********************************************************************
  //
  // --AlbertaGrid
  // --Grid
  //
  //**********************************************************************

  template <int dim, int dimworld>
  struct AlbertaGridFamily
  {
    typedef AlbertaGrid<dim,dimworld> GridImp;

    typedef DefaultLevelIndexSet< AlbertaGrid<dim,dimworld> > LevelIndexSetImp;
    typedef DefaultLeafIndexSet< AlbertaGrid<dim,dimworld> > LeafIndexSetImp;

    typedef AlbertaGridIdSet<dim,dimworld> IdSetImp;
    typedef int IdType;

    struct Traits
    {
      typedef GridImp Grid;

      typedef Dune :: Intersection< const GridImp, LeafIntersectionIteratorWrapper > LeafIntersection;
      typedef Dune :: Intersection< const GridImp, LeafIntersectionIteratorWrapper > LevelIntersection;
      typedef Dune::IntersectionIterator<const GridImp, LeafIntersectionIteratorWrapper, LeafIntersectionIteratorWrapper > LeafIntersectionIterator;
      typedef Dune::IntersectionIterator<const GridImp, LeafIntersectionIteratorWrapper, LeafIntersectionIteratorWrapper > LevelIntersectionIterator;

      typedef Dune::HierarchicIterator<const GridImp, AlbertaGridHierarchicIterator> HierarchicIterator;

      typedef IdType GlobalIdType;
      typedef IdType LocalIdType;

      template <int cd>
      struct Codim
      {
        // IMPORTANT: Codim<codim>::Geometry == Geometry<dim-codim,dimw>
        typedef Dune::Geometry<dim-cd, dimworld, const GridImp, AlbertaGridGeometry> Geometry;
        typedef Dune::Geometry<dim-cd, dim, const GridImp, AlbertaGridGeometry> LocalGeometry;
        // we could - if needed - introduce an other struct for dimglobal of Geometry

        typedef typename SelectEntityImp<cd,dim,GridImp>::Entity Entity;

        typedef Dune::LevelIterator<cd,All_Partition,const GridImp,AlbertaGridLevelIterator> LevelIterator;

        typedef Dune::LeafIterator<cd,All_Partition,const GridImp,AlbertaGridLeafIterator> LeafIterator;

        typedef AlbertaGridEntityPointer< cd, const GridImp > EntityPointerImpl;
        typedef Dune::EntityPointer< const GridImp, EntityPointerImpl > EntityPointer;

        template <PartitionIteratorType pitype>
        struct Partition
        {
          typedef Dune::LevelIterator<cd,pitype,const GridImp,AlbertaGridLevelIterator> LevelIterator;
          typedef Dune::LeafIterator<cd,pitype,const GridImp,AlbertaGridLeafIterator> LeafIterator;
        };

      };

      template <PartitionIteratorType pitype>
      struct Partition
      {
        typedef Dune::GridView<DefaultLevelGridViewTraits<const GridImp,pitype> >
        LevelGridView;
        typedef Dune::GridView<DefaultLeafGridViewTraits<const GridImp,pitype> >
        LeafGridView;
      };

      typedef IndexSet<GridImp,LevelIndexSetImp,DefaultLevelIteratorTypes<GridImp> > LevelIndexSet;
      typedef IndexSet<GridImp,LeafIndexSetImp,DefaultLeafIteratorTypes<GridImp> > LeafIndexSet;
      typedef IdSet<GridImp,IdSetImp,IdType> GlobalIdSet;
      typedef IdSet<GridImp,IdSetImp,IdType> LocalIdSet;

      //#if HAVE_MPI
      // use collective communciation with MPI
      //      typedef CollectiveCommunication<MPI_Comm> CollectiveCommunication;
      //#else
      // use dummy collective communication
      typedef Dune :: CollectiveCommunication< GridImp >
      CollectiveCommunication;
      //#endif
    };
  };

  /**
     \class AlbertaGrid

     \brief [<em> provides \ref Dune::Grid </em>]
     \brief Provides the simplicial meshes of the finite element tool box
     \brief ALBERTA (http://www.alberta-fem.de/)
     \brief written by Kunibert Siebert and Alfred Schmidt.
     \ingroup GridImplementations
     \ingroup AlbertaGrid

     This is one implementation of the grid interface using the
     the finite elemente tool box ALBERTA ( ALBERTA was written
     by Alfred Schmidt and Kunibert G. Siebert, see http://www.alberta-fem.de/).
     ALBERTA provides simplex meshes in 1d, 2d, and 3d space dimensions.
     Also the ALBERTA meshes can be dynamically adapted using
     bisection algorithms.

     The actual version of ALBERTA 1.2 can be
     downloaded at http://www.alberta-fem.de/.
     After installing the lib to a path of your choise (the PATH_TO_ALBERTA)
     you can use %Dune with this package and you have only to deliver
     the path with the --with-alberta option to %Dune.

     For using %Dune with ALBERTA you must tell %Dune where to find ALBERTA,
     which dimension to use and which dimension your world should have, i.e:

     <tt> ./autogen.sh [OPTIONS]
     --with-alberta=PATH_TO_ALBERTA and
     --with-alberta-dim=DIM --with-alberta-world-dim=DIMWORLD
     </tt>

     Now you must use the AlbertaGrid with DIM and DIMWORLD, otherwise
     unpredictable results may occur. The alberta-dim value is defaulted to 2
     if --with-alberta-dim is not provided and --with-alberta-world-dim is
     defaulted to alberta-dim if not provided. If the --with-grid-dim (see
     DGF Parser's gridtype.hh) is provided the alberta-dim,
     if not provided, is defaulted to the value of grid-dim.


     \e Note: Although ALBERTA supports different combination of DIM <= DIMWORLD,
        so far only the combinations \c DIM=DIMWORLD=2 and \c DIM=DIMWORLD=3
        are supported.

     For installation instructions see http://www.dune-project.org/doc/contrib-software.html#alberta .
   */
  template <int dim, int dimworld>
  class AlbertaGrid :
    public GridDefaultImplementation <dim,dimworld,albertCtype, AlbertaGridFamily<dim,dimworld> >,
    public HasObjectStream ,
    public HasHierarchicIndexSet
  {
    friend class AlbertaGridEntity <0,dim,const AlbertaGrid<dim,dimworld> >;
    friend class AlbertaGridEntity <1,dim,const AlbertaGrid<dim,dimworld> >;
    friend class AlbertaGridEntity <2,dim,const AlbertaGrid<dim,dimworld> >;
    friend class AlbertaGridEntity <dim,dim,const AlbertaGrid<dim,dimworld> >;

    friend class AlbertaGridEntityPointer <0,const AlbertaGrid<dim,dimworld> >;
    friend class AlbertaGridEntityPointer <1,const AlbertaGrid<dim,dimworld> >;
    friend class AlbertaGridEntityPointer <2,const AlbertaGrid<dim,dimworld> >;
    friend class AlbertaGridEntityPointer <3,const AlbertaGrid<dim,dimworld> >;

    // friends because of fillElInfo
    friend class AlbertaGridTreeIterator<0,All_Partition,AlbertaGrid<dim,dimworld> >;
    friend class AlbertaGridTreeIterator<1,All_Partition,AlbertaGrid<dim,dimworld> >;
    friend class AlbertaGridTreeIterator<2,All_Partition,AlbertaGrid<dim,dimworld> >;
    friend class AlbertaGridTreeIterator<3,All_Partition,AlbertaGrid<dim,dimworld> >;

    friend class AlbertaGridHierarchicIterator<AlbertaGrid<dim,dimworld> >;

    friend class AlbertaGridIntersectionIterator<AlbertaGrid<dim,dimworld> >;
    friend class AlbertaGridIntersectionIterator<const AlbertaGrid<dim,dimworld> >;

    //! AlbertaGrid is only implemented for 2 and 3 dimension
    //! for 1d use SGrid or SimpleGrid
    //CompileTimeChecker<dimworld != 1>   Do_not_use_AlbertaGrid_for_1d_Grids;

    typedef AlbertaGrid<dim,dimworld> MyType;

    friend class AlbertaMarkerVector;
    friend class AlbertaGridHierarchicIndexSet<dim,dimworld>;

    // minimum number of elements assumed to be created during adaption
    enum { defaultElementChunk_ = 100 };

    //**********************************************************
    // The Interface Methods
    //**********************************************************
  public:
    //! the grid family of AlbertaGrid
    typedef AlbertaGridFamily<dim,dimworld> GridFamily;
    typedef GridDefaultImplementation <dim,dimworld,albertCtype,
        AlbertaGridFamily<dim,dimworld> > BaseType;

    //! the Traits
    typedef typename AlbertaGridFamily<dim,dimworld> :: Traits Traits;

    //! type of hierarchic index set
    typedef AlbertaGridHierarchicIndexSet<dim,dimworld> HierarchicIndexSet;

  private:
    //! type of communication class
    typedef typename Traits:: CollectiveCommunication CollectiveCommunicationType;

    //! type of LeafIterator
    typedef typename Traits::template Codim<0>::LeafIterator LeafIterator;

    //! impl types of iterators
    typedef typename GridFamily:: LevelIndexSetImp LevelIndexSetImp;
    typedef typename GridFamily:: LeafIndexSetImp LeafIndexSetImp;

    //! type of leaf index set
    typedef typename Traits :: LeafIndexSet LeafIndexSet;

    //! id set impl
    typedef AlbertaGridIdSet<dim,dimworld> IdSetImp;
    typedef typename Traits :: GlobalIdSet GlobalIdSet;
    typedef typename Traits :: LocalIdSet LocalIdSet;

  public:
    //! type of leaf data
    typedef typename ALBERTA AlbertHelp::AlbertLeafData<dimworld,dim+1> LeafDataType;

    class ObjectStream
    {
    public:
      class EOFException {} ;
      template <class T>
      void readObject (T &) {}
      void readObject (int) {}
      void readObject (double) {}
      template <class T>
      void writeObject (T &) {}
      void writeObject (int) {}
      void writeObject (double) {}

      template <class T>
      void read (T &) const {}
      template <class T>
      void write (const T &) {}
    };

    typedef ObjectStream ObjectStreamType;

    enum {
      //! \brief we always have dim+1 codimensions
      numCodim = dim+1
    };

    enum {
      //! \brief max number of allowed levels is 64
      MAXL = 64
    };

    /*
       levInd = true means that a consecutive level index is generated
       if levInd == true the the element number of first macro element is
       set to 1 so hasLevelIndex_ can be identified we grid is read from
       file */
    /** \brief Constructor which reads an ALBERTA macro triangulation file.*/
    AlbertaGrid(const std::string macroTriangFilename);

    /* (for internal use only)
       Constructor which reads an ALBERTA macro triangulation file
       or given GridFile , proc is the number of domain ,
       levInd = true means that a consecutive level index is generated
       if levInd == true the the element number of first macro element is
       set to 1 so hasLevelIndex_ can be identified we grid is read from
       file */
    //AlbertaGrid(AlbertaGrid<dim,dimworld> & oldGrid, int proc);

    //! \brief empty Constructor
    AlbertaGrid();

    //! \brief Desctructor
    ~AlbertaGrid();

    //! Return maximum level defined in this grid. Levels are numbered
    //! 0 ... maxLevel with 0 the coarsest level.
    int maxLevel() const;

    //! Iterator to first entity of given codim on level
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LevelIterator
    lbegin (int level) const;

    //! one past the end on this level
    template<int cd, PartitionIteratorType pitype>
    typename Traits::template Codim<cd>::template Partition<pitype>::LevelIterator
    lend (int level) const;

    //! Iterator to first entity of given codim on level
    template<int cd>  typename Traits::template Codim<cd>::
    template Partition<All_Partition>::LevelIterator
    lbegin (int level) const;

    //! one past the end on this level
    template<int cd>  typename Traits::template Codim<cd>::
    template Partition<All_Partition>::LevelIterator
    lend (int level) const;

    //! return LeafIterator which points to first leaf entity
    template <int codim, PartitionIteratorType pitype>
    typename Traits::template Codim<codim>::template Partition<pitype>::LeafIterator
    leafbegin () const;

    //! return LeafIterator which points to first leaf entity
    template <int codim>
    typename Traits::template Codim<codim>::LeafIterator
    leafbegin () const;

    //! return LeafIterator which points behind last leaf entity
    template <int codim, PartitionIteratorType pitype>
    typename Traits::template Codim<codim>::template Partition<pitype>::LeafIterator
    leafend   () const;

    //! return LeafIterator which points behind last leaf entity
    template <int codim>
    typename Traits::template Codim<codim>::LeafIterator
    leafend   () const;

  private:
    //! return LeafIterator which points to first leaf entity
    template <int codim, PartitionIteratorType pitype>
    typename Traits::template Codim<codim>::template Partition<pitype>::LeafIterator
    leafbegin ( int maxlevel, int proc = -1 ) const;

    //! return LeafIterator which points to first leaf entity
    template <int codim>
    typename Traits::template Codim<codim>::LeafIterator
    leafbegin ( int maxlevel, int proc = -1 ) const;

    //! return LeafIterator which points behind last leaf entity
    template <int codim, PartitionIteratorType pitype>
    typename Traits::template Codim<codim>::template Partition<pitype>::LeafIterator
    leafend   ( int maxlevel, int proc = -1 ) const;

    //! return LeafIterator which points behind last leaf entity
    template <int codim>
    typename Traits::template Codim<codim>::LeafIterator
    leafend   ( int maxlevel, int proc = -1 ) const;

    //! return LeafIterator which points to first leaf entity
    LeafIterator leafbegin ( int maxlevel, int proc = -1 ) const;

    //! return LeafIterator which points behind last leaf entity
    LeafIterator leafend   ( int maxlevel, int proc = -1 ) const;

    //! return LeafIterator which points to first leaf entity
    LeafIterator leafbegin () const;

    //! return LeafIterator which points behind last leaf entity
    LeafIterator leafend   () const;

  public:
    /** \brief Number of grid entities per level and codim
     * because lbegin and lend are none const, and we need this methods
     * counting the entities on each level, you know.
     */
    int size (int level, int codim) const;

    //! number of entities per level and geometry type in this process
    int size (int level, GeometryType type) const;

    //! number of leaf entities per codim in this process
    int size (int codim) const;

    //! number of leaf entities per geometry type in this process
    int size (GeometryType type) const;

  public:
    //***************************************************************
    //  Interface for Adaptation
    //***************************************************************
    //! @copydoc Dune::Grid::mark
    bool mark( int refCount , const typename Traits::template Codim<0>::EntityPointer & en ) const DUNE_DEPRECATED;

    //! @copydoc Dune::Grid::getMark
    int getMark( const typename Traits::template Codim<0>::EntityPointer & ) const DUNE_DEPRECATED;

    //! @copydoc Dune::Grid::getMark
    int getMark( const typename Traits::template Codim<0>::Entity & ) const;

    //! @copydoc Dune::Grid::mark
    bool mark( int refCount , const typename Traits::template Codim<0>::Entity & en ) const;

  public:
    //! uses the interface, mark on entity and refineLocal
    bool globalRefine(int refCount);

    /*! \brief refine all positive marked leaf entities,
    *  coarsen all negative marked entities if possible,
    *  return true if a least one element was refined */
    bool adapt ( );

    //! adapt method with DofManager
    template <class DofManagerType, class RestrictProlongOperatorType>
    bool adapt (DofManagerType &, RestrictProlongOperatorType &, bool verbose=false );

    //! returns true, if a least one element is marked for coarsening
    bool preAdapt ();

    //! clean up some markers
    bool postAdapt();

    /** \brief return reference to collective communication, if MPI found
     * this is specialisation for MPI */
    const CollectiveCommunicationType & comm () const
    {
      return comm_;
    }

    /** \brief return name of the grid */
    std::string name () const { return "AlbertaGrid"; };

    //**********************************************************
    // End of Interface Methods
    //**********************************************************
    /** \brief write Grid to file in specified GrapeIOFileFormatType */
    template <GrapeIOFileFormatType ftype>
    bool writeGrid( const std::basic_string<char> filename, albertCtype time ) const;

    /** \brief read Grid from file filename and store time of mesh in time */
    template <GrapeIOFileFormatType ftype>
    bool readGrid( const std::basic_string<char> filename, albertCtype & time );

    /* returns size of mesh include all levels
       max Index of grid entities with given codim
       for outside the min index is 0, the shift has to done inside
       the grid which is of minor cost
     */
    int global_size (int codim) const;

    // return number of my processor
    int myRank () const { return myRank_; };

    //! transform grid N = scalar * x + trans
    void setNewCoords(const FieldVector<albertCtype, dimworld> & trans, const albertCtype scalar);

    // return hierarchic index set
    const HierarchicIndexSet & hierarchicIndexSet () const { return hIndexSet_; }

    //! return level index set for given level
    const typename Traits :: LevelIndexSet & levelIndexSet (int level) const;

    //! return leaf index set
    const typename Traits :: LeafIndexSet & leafIndexSet () const;

    //! return global IdSet
    const GlobalIdSet & globalIdSet () const { return globalIdSet_; }

    //! return local IdSet
    const LocalIdSet & localIdSet () const { return globalIdSet_; }

    // access to mesh pointer, needed by some methods
    ALBERTA MESH* getMesh () const { return mesh_; };

    // return real entity implementation
    template <int cd>
    AlbertaGridEntity<cd,dim,const AlbertaGrid<dim,dimworld> >&
    getRealEntity(typename Traits::template Codim<cd>::Entity& entity)
    {
      return this->getRealImplementation(entity);
    }

  private:
    //! return real entity implementation
    template <int cd>
    const AlbertaGridEntity<cd,dim,const AlbertaGrid<dim,dimworld> >&
    getRealEntity(const typename Traits::template Codim<cd>::Entity& entity) const
    {
      return this->getRealImplementation(entity);
    }

  public:
    //! returns geometry type vector for codimension
    const std::vector < GeometryType > & geomTypes (int codim) const
    {
      assert( codim >= 0 );
      assert( codim < dim+1 );
      return geomTypes_[codim];
    }

  private:
    friend class Conversion<AlbertaGrid<dim, dimworld>, HasObjectStream>;
    friend class Conversion<const AlbertaGrid<dim, dimworld>, HasObjectStream>;

    friend class Conversion<AlbertaGrid<dim, dimworld>, HasHierarchicIndexSet >;
    friend class Conversion<const AlbertaGrid<dim, dimworld>, HasHierarchicIndexSet>;

    // do not use copy constructor
    AlbertaGrid(const MyType& other);
    // do not use assigment
    MyType& operator=(const MyType& other);

  private:
    typedef std::vector<int> ArrayType;

    ArrayType ghostFlag_; // store ghost information

    // initialize of some members
    void initGrid(int proc);

    // make the calculation of indexOnLevel and so on.
    // extra method because of Reihenfolge
    void calcExtras();

    // write ALBERTA mesh file
    bool writeGridXdr  ( const std::basic_string<char> filename, albertCtype time ) const;

    //! reads ALBERTA mesh file
    bool readGridXdr   ( const std::basic_string<char> filename, albertCtype & time );

    //! reads ALBERTA macro file
    bool readGridAscii ( const std::basic_string<char> filename, albertCtype & time );

    // delete mesh and all vectors
    void removeMesh();

    // pointer to an Albert Mesh, which contains the data
    ALBERTA MESH *mesh_;

    // object of collective communication
    CollectiveCommunicationType comm_;

    // number of maxlevel of the mesh
    int maxlevel_;

    // true if grid was refined or coarsend
    bool wasChanged_;

    // help vector for setNewCoords
    mutable ArrayType macroVertices_;

  public:
    // this method is new fill_elinfo from ALBERTA but here the neighbor
    // relations are calced diffrent, on ervery level there are neighbor
    // realtions ( in ALBERTA only on leaf level ), so we needed a new
    // fill_elinfo.
    void fillElInfo(int ichild, int actLevel ,const ALBERTA EL_INFO *elinfo_old,
                    ALBERTA EL_INFO *elinfo, bool hierachical, bool leaf=false ) const;

    // calc the neigh[0]
    void firstNeigh(const int ichild,const ALBERTA EL_INFO *elinfo_old,
                    ALBERTA EL_INFO *elinfo, const bool leafLevel) const;

    // calc the neigh[1]
    void secondNeigh(const int ichild, const ALBERTA EL_INFO *elinfo_old,
                     ALBERTA EL_INFO *elinfo, const bool leafLevel) const;

    // calc the neigh[2]
    void thirdNeigh(const int ichild, const ALBERTA EL_INFO *elinfo_old,
                    ALBERTA EL_INFO *elinfo, const bool leafLevel) const;

  private:
    // needed for VertexIterator, mark on which element a vertex is treated
    mutable AlbertaMarkerVector vertexMarkerLeaf_;

    // needed for VertexIterator, mark on which element a vertex is treated
    mutable AlbertaMarkerVector vertexMarkerLevel_[MAXL];

    //***********************************************************************
    //  MemoryManagement for Entitys and Geometrys
    //**********************************************************************
    typedef typename SelectEntityImp<0,dim,const MyType>::EntityObject EntityObject;

  public:
    typedef AGMemoryProvider< EntityObject > EntityProvider;

    typedef AlbertaGridIntersectionIterator< const MyType > IntersectionIteratorImp;
    typedef IntersectionIteratorImp LeafIntersectionIteratorImp;
    typedef AGMemoryProvider< LeafIntersectionIteratorImp > LeafIntersectionIteratorProviderType;
    friend class LeafIntersectionIteratorWrapper< const MyType > ;

    typedef LeafIntersectionIteratorWrapper<const MyType >
    AlbertaGridIntersectionIteratorType;

    LeafIntersectionIteratorProviderType & leafIntersetionIteratorProvider() const { return leafInterItProvider_; }

  private:
    mutable EntityProvider entityProvider_;
    mutable LeafIntersectionIteratorProviderType leafInterItProvider_;

  public:
    template< class IntersectionInterfaceType >
    const typename BaseType
    :: template ReturnImplementationType< IntersectionInterfaceType >
    :: ImplementationType & DUNE_DEPRECATED
    getRealIntersectionIterator ( const IntersectionInterfaceType &iterator ) const
    {
      return this->getRealImplementation( iterator );
    }

    template< class IntersectionType >
    const typename BaseType
    :: template ReturnImplementationType< IntersectionType >
    :: ImplementationType &
    getRealIntersection ( const IntersectionType &intersection ) const
    {
      return this->getRealImplementation( intersection );
    }

    // (for internal use only) return obj pointer to EntityImp
    template <int codim>
    typename SelectEntityImp<codim,dim,const MyType>::EntityObject *
    getNewEntity (int level, bool leafIt ) const;

    // (for internal use only) free obj pointer of EntityImp
    template <int codim>
    void freeEntity (typename SelectEntityImp<codim,dim,const MyType>::EntityObject * en) const;

  private:
    //*********************************************************************
    // organisation of the global index
    //*********************************************************************
    // provides the indices for the elements
    IndexManagerType indexStack_[AlbertHelp::numOfElNumVec];

    // the DOF_INT_VECs we need
    // * change to mutable here
    mutable ALBERTA AlbertHelp::DOFVEC_STACK dofvecs_;

    const ALBERTA DOF_ADMIN * elAdmin_;

    const ALBERTA REAL_D * coordsVec_;

    // pointer to vec of elNumbers_
    const int * elNewVec_;

    // for access in the elNewVec and ownerVec
    const int nv_;
    const int dof_;

  public:
    // make some shortcuts
    void arrangeDofVec();

    // return true if el is new
    bool checkElNew ( const ALBERTA EL * el ) const;

    // read global element number from elNumbers_
    const ALBERTA REAL_D & getCoord ( const ALBERTA EL_INFO * elInfo, int vx ) const
    {
      assert( vx>= 0);
      assert( vx < dim+1 );
#ifndef CALC_COORD
      assert(coordsVec_);
      return coordsVec_[elInfo->el->dof[vx][0]];
#else
      return elInfo->coord[vx];
#endif
    }

    // read level from elNewCehck vector
    int getLevelOfElement ( const ALBERTA EL * el ) const;

    // read global element number from elNumbers_
    int getElementNumber ( const ALBERTA EL * el ) const;

    // read global element number from elNumbers_
    int getEdgeNumber ( const ALBERTA EL * el, int edge ) const;

    // read global element number from elNumbers_
    int getFaceNumber ( const ALBERTA EL * el, int face ) const;

    // read global element number from elNumbers_
    int getVertexNumber ( const ALBERTA EL * el, int vx ) const;

  private:
    // rank of my thread, i.e. number of my processor
    const int myRank_;

    // the hierarchical numbering of AlbertaGrid, unique per codim and processor
    AlbertaGridHierarchicIndexSet<dim,dimworld> hIndexSet_;

    // the id set of this grid
    IdSetImp globalIdSet_;

    // the level index set, is generated from the HierarchicIndexSet
    // is generated, when accessed
    mutable std::vector < LevelIndexSetImp * > levelIndexVec_;

    // the leaf index set, is generated from the HierarchicIndexSet
    // is generated, when accessed
    mutable LeafIndexSetImp* leafIndexSet_;

    //! stores geometry types of this grid
    std::vector < std::vector< GeometryType > > geomTypes_;

    // creates geomType_ vector
    void makeGeomTypes();

    // stack for storing BOUNDARY objects created during mesh creation
    std::stack < BOUNDARY * > bndStack_;

    typedef SingleTypeSizeCache<MyType> SizeCacheType;
    SizeCacheType * sizeCache_;

    // count how much elements where marked
    mutable int coarsenMarked_;
    mutable int refineMarked_;

    mutable bool lockPostAdapt_;
  }; // end class AlbertaGrid


  namespace Capabilities
  {

    /** \struct isParallel
       \ingroup AlbertaGrid
     */

    /** \struct IsUnstructured
       \ingroup AlbertaGrid
     */

    /** \brief AlbertaGrid has entities for all codimension
       \ingroup AlbertaGrid
     */
    template<int dim, int dimw, int cdim>
    struct hasEntity<AlbertaGrid<dim,dimw>, cdim >
    {
      static const bool v = true;
    };

    /** \brief AlbertaGrid is not levelwise conforming (since it uses bisection)
       \ingroup AlbertaGrid
     */
    template<int dim, int dimw>
    struct isLevelwiseConforming< AlbertaGrid<dim,dimw> >
    {
      static const bool v = false;
    };

    /** \brief AlbertaGrid is leafwise conforming
       \ingroup AlbertaGrid
     */
    template<int dim, int dimw>
    struct isLeafwiseConforming< AlbertaGrid<dim,dimw> >
    {
      static const bool v = true;
    };

    /** \brief AlbertaGrid does not support hanging nodes
       \ingroup AlbertaGrid
     */
    template<int dim, int dimw>
    struct hasHangingNodes< AlbertaGrid<dim,dimw> >
    {
      static const bool v = false;
    };

    /** \brief AlbertaGrid has backup and restore facilities
       \ingroup AlbertaGrid
     */
    template<int dim, int dimw>
    struct hasBackupRestoreFacilities< AlbertaGrid<dim,dimw> >
    {
      static const bool v = true;
    };

  } // end namespace Capabilities

} // namespace Dune

#include "agmemory.hh"
#include "albertagrid.cc"

// undef all dangerous defines
#undef DIM
#undef DIM_OF_WORLD
#undef CALC_COORD
#include "alberta_undefs.hh"

#endif // HAVE_ALBERTA

#endif
