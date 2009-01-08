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

#include <dune/grid/albertagrid/misc.hh>
#include <dune/grid/albertagrid/capabilities.hh>

// contains a simple memory management for some componds of this grid
#include "agmemory.hh"

#include "elementinfo.hh"

namespace Dune
{
  // i.e. double or float
  typedef Alberta::Real albertCtype;
}

#include "referencetopo.hh"
#include "indexsets.hh"
#include "geometry.hh"
#include "entity.hh"
#include "entitypointer.hh"
#include "hierarchiciterator.hh"
#include "treeiterator.hh"
#include "leveliterator.hh"
#include "leafiterator.hh"
#include "intersection.hh"

namespace Dune
{
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

    typedef AlbertaGridIdSet< dim, dimworld > IdSetImp;
    typedef unsigned int IdType;

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

        typedef typename SelectEntityImp<cd,dim,GridImp>::Entity Entity;

        typedef AlbertaGridEntityPointer< cd, const GridImp > EntityPointerImpl;
        typedef Dune::EntityPointer< const GridImp, EntityPointerImpl > EntityPointer;

        template <PartitionIteratorType pitype>
        struct Partition
        {
          typedef Dune::LevelIterator<cd,pitype,const GridImp,AlbertaGridLevelIterator> LevelIterator;
          typedef Dune::LeafIterator<cd,pitype,const GridImp,AlbertaGridLeafIterator> LeafIterator;
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

      typedef IndexSet<GridImp,LevelIndexSetImp,DefaultLevelIteratorTypes<GridImp> > LevelIndexSet;
      typedef IndexSet<GridImp,LeafIndexSetImp,DefaultLeafIteratorTypes<GridImp> > LeafIndexSet;
      typedef IdSet<GridImp,IdSetImp,IdType> GlobalIdSet;
      typedef IdSet<GridImp,IdSetImp,IdType> LocalIdSet;

      typedef Dune::CollectiveCommunication< int > CollectiveCommunication;
    };
  };

  /** \class AlbertaGrid
   *  \brief [<em> provides \ref Dune::Grid </em>]
   *  \brief simplicial grid imlementation from the ALBERTA finite element
   *         toolbox
   *  \ingroup GridImplementations
   *  \ingroup AlbertaGrid
   *
   *  AlbertaGrid provides access to the grid from the ALBERTA finite element
   *  toolbox through the %Dune interface.
   *
   *  ALBERTA is a finite element toolbox written by Alfred Schmidt and
   *  Kunibert G. Siebert (see http://www.alberta-fem.de). It contains a
   *  simplicial mesh in 1, 2 and 3 space dimensions that can be dynamically
   *  adapted by a bisection algorithm.
   *
   *  Supported ALBERTA versions include 1.2 and 2.0. Both versions can be
   *  downloaded from the ALBERTA website (www.alberta-fem.de). After
   *  installing ALBERTA, just configure DUNE with the --with-alberta option
   *  and provide the path to ALBERTA. You also have to specify which
   *  dimensions of grid and world to use. For example, your %Dune configure
   *  options could contain the following settings
   *  \code
   *  --with-alberta=ALBERTAPATH
   *  --with-alberta-dim=DIMGRID
   *  --with-alberta-world-dim=DIMWORLD
   *  \endcode
   *  The default values are <tt>DIMGRID</tt>=2 and
   *  <tt>DIMWORLD</tt>=<tt>DIMGRID</tt>.
   *  If the <tt>--with-grid-dim</tt> (see DGF Parser's gridtype.hh) is
   *  provided, <tt>DIMGRID</tt> will default to this value.
   *  You can then use <tt>AlbertaGrid< DIMGRID, DIMWORLD ></tt>.
   *  Using other template parameters might result in unpredictable behavior.
   *
   *  Further installation instructions can be found here:
   *  http://www.dune-project.org/doc/contrib-software.html#alberta
   *
   *  \note Although ALBERTA supports different combinations of
   *        <tt>DIMGRID</tt><=<tt>DIMWORLD</tt>, so far only the
   *        case <tt>DIMGRID</tt>=<tt>DIMWORLD</tt> is supported.
   */
  template< int dim, int dimworld = DIM_OF_WORLD >
  class AlbertaGrid
    : public GridDefaultImplementation
      < dim, dimworld, Alberta::Real, AlbertaGridFamily< dim, dimworld > >,
      public HasObjectStream,
      public HasHierarchicIndexSet
  {
    typedef AlbertaGrid< dim, dimworld > This;
    typedef GridDefaultImplementation
    < dim, dimworld, Alberta::Real, AlbertaGridFamily< dim, dimworld > >
    Base;

    // make Conversion a friend
    template< class, class > friend class Conversion;

    friend class AlbertaGridEntity <0,dim,const AlbertaGrid<dim,dimworld> >;
    friend class AlbertaGridEntity <1,dim,const AlbertaGrid<dim,dimworld> >;
    friend class AlbertaGridEntity <2,dim,const AlbertaGrid<dim,dimworld> >;
    friend class AlbertaGridEntity <dim,dim,const AlbertaGrid<dim,dimworld> >;

    friend class AlbertaGridEntityPointer <0,const AlbertaGrid<dim,dimworld> >;
    friend class AlbertaGridEntityPointer <1,const AlbertaGrid<dim,dimworld> >;
    friend class AlbertaGridEntityPointer <2,const AlbertaGrid<dim,dimworld> >;
    friend class AlbertaGridEntityPointer <3,const AlbertaGrid<dim,dimworld> >;

    friend class AlbertaGridHierarchicIterator<AlbertaGrid<dim,dimworld> >;

    friend class AlbertaGridIntersectionIterator<AlbertaGrid<dim,dimworld> >;
    friend class AlbertaGridIntersectionIterator<const AlbertaGrid<dim,dimworld> >;

    //! AlbertaGrid is only implemented for 2 and 3 dimension
    //! for 1d use SGrid or SimpleGrid
    //CompileTimeChecker<dimworld != 1>   Do_not_use_AlbertaGrid_for_1d_Grids;

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

    // the Traits
    typedef typename AlbertaGridFamily< dim, dimworld >::Traits Traits;

    typedef Alberta::Real ctype;

    static const int dimension = dim;
    static const int dimensionworld = dimworld;

    //! type of hierarchic index set
    typedef AlbertaGridHierarchicIndexSet<dim,dimworld> HierarchicIndexSet;

    //! type of collective communication
    typedef typename Traits::CollectiveCommunication CollectiveCommunication;

  private:
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

  private:
    typedef Alberta::MeshPointer< dimension > MeshPointer;

    // forbid copying and assignment
    AlbertaGrid ( const This & );
    This &operator= ( const This & );

  public:
    /*
       levInd = true means that a consecutive level index is generated
       if levInd == true the the element number of first macro element is
       set to 1 so hasLevelIndex_ can be identified we grid is read from
       file */
    /** \brief Constructor which reads an ALBERTA macro triangulation file.*/
    AlbertaGrid ( const std::string &macroGridFileName,
                  const std::string &gridName = "AlbertaGrid" );

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
    template< int codim >
    typename Traits::template Codim< codim >::LevelIterator
    lbegin ( int level ) const;

    //! one past the end on this level
    template< int codim >
    typename Traits::template Codim< codim >::LevelIterator
    lend ( int level ) const;

    //! return LeafIterator which points to first leaf entity
    template< int codim, PartitionIteratorType pitype >
    typename Traits
    ::template Codim< codim >::template Partition< pitype >::LeafIterator
    leafbegin () const;

    //! return LeafIterator which points behind last leaf entity
    template< int codim, PartitionIteratorType pitype >
    typename Traits
    ::template Codim< codim >::template Partition< pitype >::LeafIterator
    leafend () const;

    //! return LeafIterator which points to first leaf entity
    template< int codim >
    typename Traits::template Codim< codim >::LeafIterator
    leafbegin () const;

    //! return LeafIterator which points behind last leaf entity
    template< int codim >
    typename Traits::template Codim< codim >::LeafIterator
    leafend () const;

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
    using Base::getMark;
    using Base::mark;

    /** \copydoc Dune::Grid::getMark(const typename Codim<0>::Entity &e) const */
    int getMark ( const typename Traits::template Codim< 0 >::Entity &e ) const;

    /** \copydoc Dune::Grid::mark(int refCount,const typename Codim<0>::Entity &e) */
    bool mark ( int refCount, const typename Traits::template Codim< 0 >::Entity &e );

    //! uses the interface, mark on entity and refineLocal
    bool globalRefine ( int refCount );

    /** \copydoc Dune::Grid::adapt() */
    bool adapt ();

    //! adapt method with DofManager
    template <class DofManagerType, class RestrictProlongOperatorType>
    bool adapt (DofManagerType &, RestrictProlongOperatorType &, bool verbose=false );

    //! returns true, if a least one element is marked for coarsening
    bool preAdapt ();

    //! clean up some markers
    bool postAdapt();

    /** \brief return reference to collective communication, if MPI found
     * this is specialisation for MPI */
    const CollectiveCommunication &comm () const
    {
      return comm_;
    }

    static std::string typeName ()
    {
      std::ostringstream s;
      s << "AlbertaGrid< " << dim << ", " << dimworld << " >";
      return s.str();
    }

    /** \brief return name of the grid */
    std::string name () const
    {
      return mesh_.name();
    };

    //**********************************************************
    // End of Interface Methods
    //**********************************************************
    /** \brief write Grid to file in specified GrapeIOFileFormatType */
    template< GrapeIOFileFormatType ftype >
    bool writeGrid( const std::string &filename, ctype time ) const;

    /** \brief read Grid from file filename and store time of mesh in time */
    template< GrapeIOFileFormatType ftype >
    bool readGrid( const std::string &filename, ctype &time );

    //! transform grid N = scalar * x + trans
    void setNewCoords(const FieldVector<albertCtype, dimworld> & trans, const albertCtype scalar);

    // return hierarchic index set
    const HierarchicIndexSet & hierarchicIndexSet () const { return hIndexSet_; }

    //! return level index set for given level
    const typename Traits :: LevelIndexSet & levelIndexSet (int level) const;

    //! return leaf index set
    const typename Traits :: LeafIndexSet & leafIndexSet () const;

    //! return global IdSet
    const GlobalIdSet &globalIdSet () const
    {
      return idSet_;
    }

    //! return local IdSet
    const LocalIdSet &localIdSet () const
    {
      return idSet_;
    }

    // access to mesh pointer, needed by some methods
    ALBERTA MESH* getMesh () const
    {
      return mesh_;
    };

    const MeshPointer &meshPointer () const
    {
      return mesh_;
    }

  private:
    using Base::getRealImplementation;

    typedef std::vector<int> ArrayType;

    // initialize of some members
    void initGrid ();

    // make the calculation of indexOnLevel and so on.
    // extra method because of Reihenfolge
    void calcExtras();

    // write ALBERTA mesh file
    bool writeGridXdr ( const std::string &filename, ctype time ) const;

    //! reads ALBERTA mesh file
    bool readGridXdr ( const std::string &filename, ctype &time );

#if 0
    //! reads ALBERTA macro file
    bool readGridAscii ( const std::basic_string<char> filename, albertCtype & time );
#endif

    // delete mesh and all vectors
    void removeMesh();

    // pointer to an Albert Mesh, which contains the data
    MeshPointer mesh_;

    // collective communication
    CollectiveCommunication comm_;

    // number of maxlevel of the mesh
    int maxlevel_;

    // true if grid was refined or coarsend
    bool wasChanged_;

    // help vector for setNewCoords
    mutable ArrayType macroVertices_;

    // needed for VertexIterator, mark on which element a vertex is treated
    mutable AlbertaMarkerVector vertexMarkerLeaf_;

    // needed for VertexIterator, mark on which element a vertex is treated
    mutable AlbertaMarkerVector vertexMarkerLevel_[MAXL];

    //***********************************************************************
    //  MemoryManagement for Entitys and Geometrys
    //**********************************************************************
    typedef typename SelectEntityImp< 0, dim, const This >::EntityObject
    EntityObject;

  public:
    typedef AGMemoryProvider< EntityObject > EntityProvider;

    typedef AlbertaGridIntersectionIterator< const This > IntersectionIteratorImp;
    typedef IntersectionIteratorImp LeafIntersectionIteratorImp;
    typedef AGMemoryProvider< LeafIntersectionIteratorImp > LeafIntersectionIteratorProviderType;
    friend class LeafIntersectionIteratorWrapper< const This >;

    typedef LeafIntersectionIteratorWrapper< const This >
    AlbertaGridIntersectionIteratorType;

    LeafIntersectionIteratorProviderType & leafIntersetionIteratorProvider() const { return leafInterItProvider_; }

  private:
    mutable EntityProvider entityProvider_;
    mutable LeafIntersectionIteratorProviderType leafInterItProvider_;

  public:
    template< class IntersectionInterfaceType >
    const typename Base
    :: template ReturnImplementationType< IntersectionInterfaceType >
    :: ImplementationType & DUNE_DEPRECATED
    getRealIntersectionIterator ( const IntersectionInterfaceType &iterator ) const
    {
      return this->getRealImplementation( iterator );
    }

    template< class IntersectionType >
    const typename Base
    :: template ReturnImplementationType< IntersectionType >
    :: ImplementationType &
    getRealIntersection ( const IntersectionType &intersection ) const
    {
      return this->getRealImplementation( intersection );
    }

    // (for internal use only) return obj pointer to EntityImp
    template< int codim >
    typename SelectEntityImp< codim, dim, const This >::EntityObject *
    getNewEntity () const;

    // (for internal use only) free obj pointer of EntityImp
    template <int codim>
    void freeEntity ( typename SelectEntityImp< codim, dim, const This >::EntityObject *en ) const;

  public:
    // make some shortcuts
    void arrangeDofVec();

    // return true if el is new
    bool checkElNew ( const Alberta::Element *element ) const;

    // read global element number from elNumbers_
    const Alberta::GlobalVector &
    getCoord ( const Alberta::ElementInfo< dimension > &elementInfo,
               int vertex ) const;

    // read level from elNewCehck vector
    int getLevelOfElement ( const Alberta::Element *element ) const;

    // read global element number from elNumbers_
    int getElementNumber ( const ALBERTA EL * el ) const;

    // read global element number from elNumbers_
    int getEdgeNumber ( const ALBERTA EL * el, int edge ) const;

    // read global element number from elNumbers_
    int getFaceNumber ( const ALBERTA EL * el, int face ) const;

    // read global element number from elNumbers_
    int getVertexNumber ( const ALBERTA EL * el, int vx ) const;

  private:
    Alberta::HierarchyDofNumbering< dimension > dofNumbering_;

    Alberta::DofVectorPointer< int > elNewCheck_;
    struct ElNewCheckInterpolation;
    class SetLocalElementLevel;

#ifndef CALC_COORD
    typedef Alberta::DofVectorPointer< Alberta::GlobalVector > CoordVectorPointer;
    CoordVectorPointer coords_;
    class SetLocalCoords;
#endif

    // for access in the elNewVec and ownerVec
    const int nv_;
    const int dof_;

    // hierarchical numbering of AlbertaGrid, unique per codim
    AlbertaGridHierarchicIndexSet<dim,dimworld> hIndexSet_;

    // the id set of this grid
    IdSetImp idSet_;

    // the level index set, is generated from the HierarchicIndexSet
    // is generated, when accessed
    mutable std::vector< LevelIndexSetImp * > levelIndexVec_;

    // the leaf index set, is generated from the HierarchicIndexSet
    // is generated, when accessed
    mutable LeafIndexSetImp* leafIndexSet_;

    typedef SingleTypeSizeCache< This > SizeCacheType;
    SizeCacheType * sizeCache_;

    // count how much elements where marked
    int coarsenMarked_;
    int refineMarked_;

    mutable bool lockPostAdapt_;
  }; // end class AlbertaGrid

} // namespace Dune

#include "agmemory.hh"
#include "albertagrid.cc"

// undef all dangerous defines
#undef DIM
#undef DIM_OF_WORLD
#undef CALC_COORD

#ifdef _ABS_NOT_DEFINED_
#undef ABS
#endif

#ifdef _MIN_NOT_DEFINED_
#undef MIN
#endif

#ifdef _MAX_NOT_DEFINED_
#undef MAX
#endif

#if DUNE_ALBERTA_VERSION >= 0x201
#ifdef obstack_chunk_alloc
#undef obstack_chunk_alloc
#endif
#ifdef obstack_chunk_free
#undef obstack_chunk_free
#endif
#include <dune/grid/albertagrid/undefine-2.1.hh>
#elif DUNE_ALBERTA_VERSION == 0x200
#include <dune/grid/albertagrid/undefine-2.0.hh>
#else
#include <dune/grid/albertagrid/undefine-1.2.hh>
#endif

#define _ALBERTA_H_

#endif // HAVE_ALBERTA

#endif
