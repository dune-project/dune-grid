// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTAGRID_IMP_HH
#define DUNE_ALBERTAGRID_IMP_HH

/** \file
 *  \author Robert Kloefkorn and Martin Nolte
 *  \brief  provides the AlbertaGrid class
 */

#if HAVE_ALBERTA

#include <iostream>
#include <fstream>
#include <dune/common/deprecated.hh>

#include <vector>
#include <assert.h>
#include <algorithm>

// Dune includes
#include <dune/common/misc.hh>
#include <dune/common/interfaces.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/stdstreams.hh>
#include <dune/common/collectivecommunication.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/adaptcallback.hh>
#include <dune/grid/common/sizecache.hh>
#include <dune/grid/common/defaultgridview.hh>

//- Local includes
// some cpp defines and include of alberta.h
#include "albertaheader.hh"

// grape data io
#include <dune/grid/utility/grapedataioformattypes.hh>

// calculate coordinates on the fly (forbidden for ALBERTA 2.0 and earlier)?
#define CALC_COORD 0

#include <dune/grid/albertagrid/misc.hh>
#include <dune/grid/albertagrid/capabilities.hh>

// contains a simple memory management for some componds of this grid
#include "agmemory.hh"

#include <dune/grid/albertagrid/coordcache.hh>
#include <dune/grid/albertagrid/level.hh>
#include <dune/grid/albertagrid/intersection.hh>
#include <dune/grid/albertagrid/intersectioniterator.hh>

#include "indexsets.hh"
#include "geometry.hh"
#include "entity.hh"
#include "entitypointer.hh"
#include "hierarchiciterator.hh"
#include "treeiterator.hh"
#include "leveliterator.hh"
#include "leafiterator.hh"

namespace Dune
{

  // InternalForward Declarations
  // ----------------------------

  template< int dim, int dimworld >
  class AlbertaGrid;



  // AlbertaGridFamily
  // -----------------

  template <int dim, int dimworld>
  struct AlbertaGridFamily
  {
    typedef AlbertaGrid<dim,dimworld> GridImp;

    typedef AlbertaGridIndexSet< dim, dimworld > LevelIndexSetImp;
    typedef AlbertaGridIndexSet< dim, dimworld > LeafIndexSetImp;
    //typedef DefaultLevelIndexSet< AlbertaGrid<dim,dimworld> > LevelIndexSetImp;
    //typedef DefaultLeafIndexSet< AlbertaGrid<dim,dimworld> > LeafIndexSetImp;

    typedef AlbertaGridIdSet< dim, dimworld > IdSetImp;
    typedef unsigned int IdType;

    struct Traits
    {
      typedef GridImp Grid;

      typedef Dune::Intersection< const GridImp, AlbertaGridLeafIntersection > LeafIntersection;
      typedef Dune::Intersection< const GridImp, AlbertaGridLeafIntersection > LevelIntersection;
      typedef Dune::IntersectionIterator
      < const GridImp, AlbertaGridLeafIntersectionIterator, AlbertaGridLeafIntersection >
      LeafIntersectionIterator;
      typedef Dune::IntersectionIterator
      < const GridImp, AlbertaGridLeafIntersectionIterator, AlbertaGridLeafIntersection >
      LevelIntersectionIterator;

      typedef Dune::HierarchicIterator<const GridImp, AlbertaGridHierarchicIterator> HierarchicIterator;

      typedef IdType GlobalIdType;
      typedef IdType LocalIdType;

      template <int cd>
      struct Codim
      {
        // IMPORTANT: Codim<codim>::Geometry == Geometry<dim-codim,dimw>
        typedef Dune::Geometry<dim-cd, dimworld, const GridImp, AlbertaGridGeometry> Geometry;
        typedef Dune::Geometry<dim-cd, dim, const GridImp, AlbertaGridGeometry> LocalGeometry;

        typedef Dune::Entity< cd, dim, const GridImp, AlbertaGridEntity > Entity;

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

      typedef IndexSet< GridImp, LevelIndexSetImp > LevelIndexSet;
      typedef IndexSet< GridImp, LeafIndexSetImp > LeafIndexSet;
      typedef AlbertaGridHierarchicIndexSet< dim, dimworld > HierarchicIndexSet;
      typedef IdSet<GridImp,IdSetImp,IdType> GlobalIdSet;
      typedef IdSet<GridImp,IdSetImp,IdType> LocalIdSet;

      typedef Dune::CollectiveCommunication< AlbertaGrid<dim,dimworld> > CollectiveCommunication;
    };
  };



  // AlbertaGrid
  // -----------

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
   *  \endcode
   *  The default value for <tt>DIMGRID</tt> is obtained from
   *  <tt>--with-world-dim</tt>.  If <tt>--with-world-dim</tt> was not given,
   *  the default is 2.
   *
   *  Further installation instructions can be found here:
   *  http://www.dune-project.org/external_libraries/install_alberta.html
   *
   *  \note Although ALBERTA supports different combinations of
   *        <tt>DIMGRID</tt><=<tt>DIMWORLD</tt>, so far only the
   *        case <tt>DIMGRID</tt>=<tt>DIMWORLD</tt> is supported.
   *
   *  If you use automake and want to compile a program maude, the following
   *  <tt>Makefile.am</tt> snippet might help:
   *  \code
   *  bin_PROGRAMS = maude
   *
   *  maude_SOURCES = maude.cc
   *  maude_CPPFLAGS = $(AM_CPPFLAGS) $(ALBERTA_CPPFLAGS)
   *  maude_LDFLAGS = $(AM_LDFLAGS) $(ALBERTA_LDFLAGS) $(DUNE_LDFLAGS)
   *  maude_LDADD = $(ALBERTA_LIBS) $(DUNE_LIBS)
   *  \endcode
   *  This will compile and link your program with the default ALBERTA
   *  dimension set with the <tt>--with-alberta-dim</tt> parameter to
   *  <tt>configure</tt>.  If you want to use a non-default dimension (or
   *  different dimension in different programs, you can use the following
   *  snippet:
   *  \code
   *  bin_PROGRAMS = maude2d maude3d
   *
   *  maude2d_SOURCES = maude.cc
   *  maude2d_CPPFLAGS = $(AM_CPPFLAGS) \
   *    $(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=2 -DENABLE_ALBERTA
   *  maude2d_LDFLAGS = $(AM_LDFLAGS) $(ALBERTA_LDFLAGS) $(DUNE_LDFLAGS)
   *  maude2d_LDADD = -lalberta_2d $(ALBERTA_BASE_LIBS) $(DUNE_LIBS)
   *
   *  maude3d_SOURCES = maude.cc
   *  maude3d_CPPFLAGS = $(AM_CPPFLAGS) \
   *    $(ALBERTA_INCLUDE_CPPFLAGS) -DALBERTA_DIM=3 -DENABLE_ALBERTA
   *  maude3d_LDFLAGS = $(AM_LDFLAGS) $(ALBERTA_LDFLAGS) $(DUNE_LDFLAGS)
   *  maude3d_LDADD = -lalberta_3d $(ALBERTA_BASE_LIBS) $(DUNE_LIBS)
   *  \endcode
   *  It is not possible to use alberta grids with different world dimensions
   *  in the same binary however.
   *
   *  In both cases you have in your program the preprocessor defines
   *  <tt>HAVE_ALBERTA</tt> which tells you whether alberta was found by
   *  configure and <tt>ALBERTA_DIM</tt> which tells you the dimension of
   *  alberta <em>for this program</em>.
   *
   *  For further details look into the <tt>alberta.m4</tt> autoconf snippet.
   */
  template< int dim, int dimworld = Alberta::dimWorld >
  class AlbertaGrid
    : public GridDefaultImplementation
      < dim, dimworld, Alberta::Real, AlbertaGridFamily< dim, dimworld > >,
      public HasHierarchicIndexSet
  {
    typedef AlbertaGrid< dim, dimworld > This;
    typedef GridDefaultImplementation
    < dim, dimworld, Alberta::Real, AlbertaGridFamily< dim, dimworld > >
    Base;

    // make Conversion a friend
    template< class, class > friend class Conversion;

    template< int, int, class > friend class AlbertaGridEntity;

    friend class AlbertaGridHierarchicIterator< This >;
    friend class AlbertaGridLeafIntersection< const This >;

    friend class AlbertaMarkerVector< dim, dimworld >;
    friend class AlbertaGridHierarchicIndexSet< dim, dimworld >;

  public:
    typedef Alberta::Real ctype;

    static const int dimension = dim;
    static const int dimensionworld = dimworld;

    //! the grid family of AlbertaGrid
    typedef AlbertaGridFamily< dim, dimworld > GridFamily;

    // the Traits
    typedef typename AlbertaGridFamily< dim, dimworld >::Traits Traits;

    //! type of leaf index set
    typedef typename Traits::LeafIndexSet LeafIndexSet;
    //! type of level index sets
    typedef typename Traits::LevelIndexSet LevelIndexSet;

    //! type of hierarchic index set
    typedef typename Traits::HierarchicIndexSet HierarchicIndexSet;

    //! type of global id set
    typedef typename Traits::GlobalIdSet GlobalIdSet;
    //! type of local id set
    typedef typename Traits::LocalIdSet LocalIdSet;

    //! type of collective communication
    typedef typename Traits::CollectiveCommunication CollectiveCommunication;

  private:
    //! type of LeafIterator
    typedef typename Traits::template Codim<0>::LeafIterator LeafIterator;

    //! id set impl
    typedef AlbertaGridIdSet<dim,dimworld> IdSetImp;

    struct AdaptationState;

    template< class DataHandler >
    struct AdaptationCallback;

    // max number of allowed levels is 64
    static const int MAXL = 64;

    typedef Alberta::ElementInfo< dimension > ElementInfo;
    typedef Alberta::MeshPointer< dimension > MeshPointer;
    typedef Alberta::HierarchyDofNumbering< dimension > DofNumbering;
    typedef AlbertaGridLevelProvider< dimension > LevelProvider;

    // forbid copying and assignment
    AlbertaGrid ( const This & );
    This &operator= ( const This & );

  public:
    /** \brief create an empty grid */
    AlbertaGrid ();

    /** \brief create a grid from an ALBERTA macro data structure
     *
     *  \param[in]  macroData  macro data to create grid from
     *  \param[in]  gridName   name of the grid (defaults to "AlbertaGrid")
     */
    AlbertaGrid ( const Alberta::MacroData< dimension > &macroData,
                  const std::string &gridName = "AlbertaGrid",
                  const DuneBoundaryProjection< dimensionworld > *projection = 0 );

    /** \brief create a grid from an ALBERTA macro grid file
     *
     *  \param[in]  macroGridFileName  name of the macro grid file
     *  \param[in]  gridName           name of the grid (defaults to "AlbertaGrid")
     */
    AlbertaGrid ( const std::string &macroGridFileName,
                  const std::string &gridName = "AlbertaGrid" );

    /** \brief desctructor */
    ~AlbertaGrid ();

    //! Return maximum level defined in this grid. Levels are numbered
    //! 0 ... maxLevel with 0 the coarsest level.
    int maxLevel () const;

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
    void globalRefine ( int refCount );

    template< class DataHandle >
    void globalRefine ( int refCount, AdaptDataHandleInterface< This, DataHandle > &handle );

    /** \copydoc Dune::Grid::adapt() */
    bool adapt ();

    //! callback adapt method with AdaptDataHandleInterface
    template< class DataHandle >
    bool adapt ( AdaptDataHandleInterface< This, DataHandle > &handle );

    //! returns true, if a least one element is marked for coarsening
    bool preAdapt ();

    //! clean up some markers
    void postAdapt();

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

    const DofNumbering &dofNumbering () const
    {
      return dofNumbering_;
    }

    const LevelProvider &levelProvider () const
    {
      return levelProvider_;
    }

    int dune2alberta ( int codim, int i ) const
    {
      return numberingMap_.dune2alberta( codim, i );
    }

    int alberta2dune ( int codim, int i ) const
    {
      return numberingMap_.alberta2dune( codim, i );
    }

    int generic2alberta ( int codim, int i ) const
    {
      return genericNumberingMap_.dune2alberta( codim, i );
    }

    int alberta2generic ( int codim, int i ) const
    {
      return genericNumberingMap_.alberta2dune( codim, i );
    }

  private:
    using Base::getRealImplementation;

    typedef std::vector<int> ArrayType;

    void setup ();

    // make the calculation of indexOnLevel and so on.
    // extra method because of Reihenfolge
    void calcExtras();

    // write ALBERTA mesh file
    bool writeGridXdr ( const std::string &filename, ctype time ) const;

    //! reads ALBERTA mesh file
    bool readGridXdr ( const std::string &filename, ctype &time );

#if 0
    //! reads ALBERTA macro file
    bool readGridAscii ( const std::string &filename, ctype &time );
#endif

    // delete mesh and all vectors
    void removeMesh();

    // pointer to an Albert Mesh, which contains the data
    MeshPointer mesh_;

    // collective communication
    CollectiveCommunication comm_;

    // number of maxlevel of the mesh
    int maxlevel_;

    //***********************************************************************
    //  MemoryManagement for Entitys and Geometrys
    //**********************************************************************
    typedef MakeableInterfaceObject< typename Traits::template Codim< 0 >::Entity >
    EntityObject;

  public:
    typedef AGMemoryProvider< EntityObject > EntityProvider;

    friend class AlbertaGridLeafIntersectionIterator< const This >;

  private:
    mutable EntityProvider entityProvider_;

  public:
    template< int codim >
    static int
    getTwist ( const typename Traits::template Codim< codim >::Entity &entity )
    {
      return getRealImplementation( entity ).twist();
    }

    template< int codim >
    static int
    getTwist ( const typename Traits::template Codim< 0 >::Entity &entity, int subEntity )
    {
      return getRealImplementation( entity ).template twist< codim >( subEntity );
    }

    static int
    getTwistInInside ( const typename Traits::LeafIntersection &intersection )
    {
      return getRealImplementation( intersection ).twistInSelf();
    }

    static int
    getTwistInOutside ( const typename Traits::LeafIntersection &intersection )
    {
      return getRealImplementation( intersection ).twistInNeighbor();
    }

    const AlbertaGridLeafIntersection< const This > &
    getRealIntersection ( const typename Traits::LeafIntersection &intersection ) const
    {
      return getRealImplementation( intersection );
    }

    template< class Intersection >
    DUNE_DEPRECATED
    const typename Base
    ::template ReturnImplementationType< Intersection >::ImplementationType &
    getRealIntersection ( const Intersection &intersection ) const
    {
      return getRealImplementation( intersection );
    }

    // (for internal use only) return obj pointer to EntityImp
    template< int codim >
    MakeableInterfaceObject< typename Traits::template Codim< codim >::Entity > *
    getNewEntity () const;

    // (for internal use only) free obj pointer of EntityImp
    template <int codim>
    void freeEntity ( MakeableInterfaceObject< typename Traits::template Codim< codim >::Entity > *entity ) const;

  public:
    // read global element number from elNumbers_
    const Alberta::GlobalVector &
    getCoord ( const ElementInfo &elementInfo, int vertex ) const;

  private:
    // map between ALBERTA and DUNE numbering
    Alberta::NumberingMap< dimension > numberingMap_;
    Alberta::NumberingMap< dimension, Alberta::Generic2AlbertaNumbering > genericNumberingMap_;

    DofNumbering dofNumbering_;

    LevelProvider levelProvider_;

    // hierarchical numbering of AlbertaGrid, unique per codim
    HierarchicIndexSet hIndexSet_;

    // the id set of this grid
    IdSetImp idSet_;

    // the level index set, is generated from the HierarchicIndexSet
    // is generated, when accessed
    mutable std::vector< typename GridFamily::LevelIndexSetImp * > levelIndexVec_;

    // the leaf index set, is generated from the HierarchicIndexSet
    // is generated, when accessed
    mutable typename GridFamily::LeafIndexSetImp* leafIndexSet_;

    typedef SingleTypeSizeCache< This > SizeCacheType;
    SizeCacheType * sizeCache_;

    typedef AlbertaMarkerVector< dim, dimworld > MarkerVector;

    // needed for VertexIterator, mark on which element a vertex is treated
    mutable MarkerVector leafMarkerVector_;

    // needed for VertexIterator, mark on which element a vertex is treated
    mutable std::vector< MarkerVector > levelMarkerVector_;

#if !CALC_COORD
    Alberta::CoordCache< dimension > coordCache_;
#endif

    // current state of adaptation
    AdaptationState adaptationState_;
  };



  // AlbertaGrid::AdaptationState
  // ----------------------------

  template< int dim, int dimworld >
  struct AlbertaGrid< dim, dimworld >::AdaptationState
  {
    enum Phase { ComputationPhase, PreAdaptationPhase, PostAdaptationPhase };

  private:
    Phase phase_;
    int coarsenMarked_;
    int refineMarked_;

  public:
    AdaptationState ()
      : phase_( ComputationPhase ),
        coarsenMarked_( 0 ),
        refineMarked_( 0 )
    {}

    void mark ( int count )
    {
      if( count < 0 )
        ++coarsenMarked_;
      if( count > 0 )
        refineMarked_ += (2 << count);
    }

    void unmark ( int count )
    {
      if( count < 0 )
        --coarsenMarked_;
      if( count > 0 )
        refineMarked_ -= (2 << count);
    }

    bool coarsen () const
    {
      return (coarsenMarked_ > 0);
    }

    int refineMarked () const
    {
      return refineMarked_;
    }

    void preAdapt ()
    {
      if( phase_ != ComputationPhase )
        error( "preAdapt may only be called in computation phase." );
      phase_ = PreAdaptationPhase;
    }

    void adapt ()
    {
      if( phase_ != PreAdaptationPhase )
        error( "adapt may only be called in preadapdation phase." );
      phase_ = PostAdaptationPhase;
    }

    void postAdapt ()
    {
      if( phase_ != PostAdaptationPhase )
        error( "postAdapt may only be called in postadaptation phase." );
      phase_ = ComputationPhase;

      coarsenMarked_ = 0;
      refineMarked_ = 0;
    }

  private:
    void error ( const std::string &message )
    {
      DUNE_THROW( InvalidStateException, message );
    }
  };

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

#if DUNE_ALBERTA_VERSION >= 0x300
#ifdef obstack_chunk_alloc
#undef obstack_chunk_alloc
#endif
#ifdef obstack_chunk_free
#undef obstack_chunk_free
#endif
#include <dune/grid/albertagrid/undefine-3.0.hh>
#else
#include <dune/grid/albertagrid/undefine-2.0.hh>
#endif

#define _ALBERTA_H_

#endif // HAVE_ALBERTA

#endif
