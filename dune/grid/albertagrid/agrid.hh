// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTAGRID_IMP_HH
#define DUNE_ALBERTAGRID_IMP_HH

/** \file
 *  \author Robert Kloefkorn and Martin Nolte
 *  \brief  provides the AlbertaGrid class
 */

#if HAVE_ALBERTA || DOXYGEN

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>

// Dune includes
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/stdstreams.hh>
#include <dune/common/parallel/communication.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/adaptcallback.hh>
#include <dune/grid/common/sizecache.hh>

//- Local includes
// some cpp defines and include of alberta.h
#include "albertaheader.hh"

#include <dune/grid/albertagrid/misc.hh>
#include <dune/grid/albertagrid/capabilities.hh>
#include <dune/grid/albertagrid/backuprestore.hh>

#include <dune/grid/albertagrid/coordcache.hh>
#include <dune/grid/albertagrid/gridfamily.hh>
#include <dune/grid/albertagrid/level.hh>
#include <dune/grid/albertagrid/intersection.hh>
#include <dune/grid/albertagrid/intersectioniterator.hh>
#include <dune/grid/albertagrid/datahandle.hh>
#include <dune/grid/albertagrid/entityseed.hh>

#include "indexsets.hh"
#include "geometry.hh"
#include "entity.hh"
#include "hierarchiciterator.hh"
#include "treeiterator.hh"
#include "leveliterator.hh"
#include "leafiterator.hh"

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class Grid >
  struct DGFGridFactory;



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
   *  Supported ALBERTA versions include 3.0 or higher.
   *  It can be downloaded from the ALBERTA website
   *  (https://gitlab.mathematik.uni-stuttgart.de/ians-nmh/alberta/alberta3).
   *
   *  After installing ALBERTA, just configure DUNE with the cmake option
   *  <tt>-DAlberta_ROOT=[path/to/alberta]</tt> and provide the path to ALBERTA.
   *
   *  Each program linking to ALBERTA only supports a fixed dimension of world.
   *  This is obtained from the <tt>ALBERTA_DIM</tt> preprocessor variable. This
   *  variable and all other necessary link flags and libraries are set by the
   *  cmake function <tt>add_dune_alberta_flags(WORLDDIM [N] [targets]...)</tt>
   *  with N the world dimension.
   *
   *  CMake sets the config.h variables <tt>HAVE_ALBERTA</tt> that tells you whether
   *  ALBERTA was found, and <tt>ALBERTA_DIM</tt> compile flags which tells you the
   *  dimension of world <em>for this program</em>.
   */
  template< int dim, int dimworld = Alberta::dimWorld >
  class AlbertaGrid
    : public GridDefaultImplementation
      < dim, dimworld, Alberta::Real, AlbertaGridFamily< dim, dimworld > >
  {
    typedef AlbertaGrid< dim, dimworld > This;
    typedef GridDefaultImplementation
    < dim, dimworld, Alberta::Real, AlbertaGridFamily< dim, dimworld > >
    Base;

    template< int, int, class > friend class AlbertaGridEntity;
    template< class > friend class AlbertaLevelGridView;
    template< class > friend class AlbertaLeafGridView;
    template< int, class, bool > friend class AlbertaGridTreeIterator;
    template< class > friend class AlbertaGridHierarchicIterator;

    friend class GridFactory< This >;
    friend struct DGFGridFactory< This >;

    friend class AlbertaGridIntersectionBase< const This >;
    friend class AlbertaGridLeafIntersection< const This >;

    friend class AlbertaMarkerVector< dim, dimworld >;
#if (__GNUC__ < 4) && !(defined __ICC)
    // add additional friend decls for gcc 3.4
    friend struct AlbertaMarkerVector< dim, dimworld >::MarkSubEntities<true>;
    friend struct AlbertaMarkerVector< dim, dimworld >::MarkSubEntities<false>;
#endif
    friend class AlbertaGridIndexSet< dim, dimworld >;
    friend class AlbertaGridHierarchicIndexSet< dim, dimworld >;

    template< class, class >
    friend class Alberta::AdaptRestrictProlongHandler;

  public:
    //! the grid family of AlbertaGrid
    typedef AlbertaGridFamily< dim, dimworld > GridFamily;

    typedef typename GridFamily::ctype ctype;

    static const int dimension = GridFamily::dimension;
    static const int dimensionworld = GridFamily::dimensionworld;

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

    //! type of communication
    typedef typename Traits::Communication Communication;

    /**
     * \deprecated Use Communication instead! Will be removed after Dune 2.9.
     */
    [[deprecated("Use Communication instead!")]]
    typedef Communication CollectiveCommunication;

  private:
    //! type of LeafIterator
    typedef typename Traits::template Codim<0>::LeafIterator LeafIterator;

    //! id set impl
    typedef AlbertaGridIdSet<dim,dimworld> IdSetImp;

    //! AdaptationState
    struct AdaptationState
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

    template< class DataHandler >
    struct AdaptationCallback;

    // max number of allowed levels is 64
    static const int MAXL = 64;

    typedef Alberta::ElementInfo< dimension > ElementInfo;
    typedef Alberta::MeshPointer< dimension > MeshPointer;
    typedef Alberta::HierarchyDofNumbering< dimension > DofNumbering;
    typedef AlbertaGridLevelProvider< dimension > LevelProvider;

  public:
    AlbertaGrid ( const This & ) = delete;
    This &operator= ( const This & ) = delete;

    /** \brief create an empty grid */
    AlbertaGrid ();

    /** \brief create a grid from an ALBERTA macro data structure
     *
     *  \param[in]  macroData   macro data to create grid from
     *  \param[in]  projection  shared pointer to a global boundary projection (defaults to 0)
     */
    AlbertaGrid ( const Alberta::MacroData< dimension > &macroData,
                  const std::shared_ptr< DuneBoundaryProjection< dimensionworld > > &projection
                    = std::shared_ptr< DuneBoundaryProjection< dimensionworld > >() );

    template< class Proj, class Impl >
    AlbertaGrid ( const Alberta::MacroData< dimension > &macroData,
                  const Alberta::ProjectionFactoryInterface< Proj, Impl > &projectionFactory );

    /** \brief create a grid from an ALBERTA macro grid file
     *
     *  \param[in]  macroGridFileName  name of the macro grid file
     */
    AlbertaGrid ( const std::string &macroGridFileName );

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

    //! number of boundary segments within the macro grid
    std::size_t numBoundarySegments () const
    {
      return numBoundarySegments_;
    }

    //! View for a grid level for All_Partition
    typename Traits::LevelGridView levelGridView ( int level ) const
    {
      typedef typename Traits::LevelGridView View;
      typedef typename View::GridViewImp ViewImp;
      return View( ViewImp( *this, level ) );
    }

    //! View for the leaf grid for All_Partition
    typename Traits::LeafGridView leafGridView () const
    {
      typedef typename Traits::LeafGridView View;
      typedef typename View::GridViewImp ViewImp;
      return View( ViewImp( *this ) );
    }

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

    /** \brief return reference to communication, if MPI found
     * this is specialisation for MPI */
    const Communication &comm () const
    {
      return comm_;
    }

    static std::string typeName ()
    {
      std::ostringstream s;
      s << "AlbertaGrid< " << dim << ", " << dimworld << " >";
      return s.str();
    }

    /** \brief obtain Entity from EntitySeed. */
    template< class EntitySeed >
    typename Traits::template Codim< EntitySeed::codimension >::Entity
    entity ( const EntitySeed &seed ) const
    {
      typedef typename Traits::template Codim< EntitySeed::codimension >::EntityImpl EntityImpl;
      return EntityImpl( *this, seed.impl().elementInfo( meshPointer() ), seed.impl().subEntity() );
    }

    //**********************************************************
    // End of Interface Methods
    //**********************************************************
    /** \brief write Grid to file in Xdr */
    bool writeGrid( const std::string &filename, ctype time ) const;

    /** \brief read Grid from file filename and store time of mesh in time */
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
    typedef std::vector<int> ArrayType;

    void setup ();

    // make the calculation of indexOnLevel and so on.
    // extra method because of Reihenfolge
    void calcExtras();

  private:
    // delete mesh and all vectors
    void removeMesh();

    //***********************************************************************
    //  MemoryManagement for Entitys and Geometrys
    //**********************************************************************
    typedef MakeableInterfaceObject< typename Traits::template Codim< 0 >::Entity >
    EntityObject;

  public:
    friend class AlbertaGridLeafIntersectionIterator< const This >;

    template< int codim >
    static int
    getTwist ( const typename Traits::template Codim< codim >::Entity &entity )
    {
      return entity.impl().twist();
    }

    template< int codim >
    static int
    getTwist ( const typename Traits::template Codim< 0 >::Entity &entity, int subEntity )
    {
      return entity.impl().template twist< codim >( subEntity );
    }

    static int
    getTwistInInside ( const typename Traits::LeafIntersection &intersection )
    {
      return intersection.impl().twistInInside();
    }

    static int
    getTwistInOutside ( const typename Traits::LeafIntersection &intersection )
    {
      return intersection.impl().twistInOutside();
    }

  public:
    // read global element number from elNumbers_
    const Alberta::GlobalVector &
    getCoord ( const ElementInfo &elementInfo, int vertex ) const;

  private:
    // pointer to an Albert Mesh, which contains the data
    MeshPointer mesh_;

    // communication
    Communication comm_;

    // maximum level of the mesh
    int maxlevel_;

    // number of boundary segments within the macro grid
    size_t numBoundarySegments_;

    // map between ALBERTA and DUNE numbering
    Alberta::NumberingMap< dimension, Alberta::Dune2AlbertaNumbering > numberingMap_;
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

    SizeCache< This > sizeCache_;

    typedef AlbertaMarkerVector< dim, dimworld > MarkerVector;

    // needed for VertexIterator, mark on which element a vertex is treated
    mutable MarkerVector leafMarkerVector_;

    // needed for VertexIterator, mark on which element a vertex is treated
    mutable std::vector< MarkerVector > levelMarkerVector_;

#if DUNE_ALBERTA_CACHE_COORDINATES
    Alberta::CoordCache< dimension > coordCache_;
#endif

    // current state of adaptation
    AdaptationState adaptationState_;
  };

} // namespace Dune

#include "albertagrid.cc"

// undef all dangerous defines
#undef DIM
#undef DIM_OF_WORLD

#ifdef _ABS_NOT_DEFINED_
#undef ABS
#endif

#ifdef _MIN_NOT_DEFINED_
#undef MIN
#endif

#ifdef _MAX_NOT_DEFINED_
#undef MAX
#endif

#ifdef obstack_chunk_alloc
#undef obstack_chunk_alloc
#endif
#ifdef obstack_chunk_free
#undef obstack_chunk_free
#endif
#include <dune/grid/albertagrid/undefine-3.0.hh>

// We use MEM_ALLOC, so undefine it here.
#undef MEM_ALLOC

// We use MEM_REALLOC, so undefine it here.
#undef MEM_REALLOC

// We use MEM_CALLOC, so undefine it here.
#undef MEM_CALLOC

// We use MEM_FREE, so undefine it here.
#undef MEM_FREE

// Macro ERROR may be defined by alberta_util.h. If so, undefine it.
#ifdef ERROR
#undef ERROR
#endif // #ifdef ERROR

// Macro ERROR_EXIT may be defined by alberta_util.h. If so, undefine it.
#ifdef ERROR_EXIT
#undef ERROR_EXIT
#endif // #ifdef ERROR_EXIT

// Macro WARNING may be defined by alberta_util.h. If so, undefine it.
#ifdef WARNING
#undef WARNING
#endif // #ifdef WARNING

// Macro TEST may be defined by alberta_util.h. If so, undefine it.
#ifdef TEST
#undef TEST
#endif // #ifdef TEST

// Macro TEST_EXIT may be defined by alberta_util.h. If so, undefine it.
#ifdef TEST_EXIT
#undef TEST_EXIT
#endif // #ifdef TEST_EXIT

// Macro DEBUG_TEST may be defined by alberta_util.h. If so, undefine it.
#ifdef DEBUG_TEST
#undef DEBUG_TEST
#endif // #ifdef DEBUG_TEST

// Macro DEBUG_TEST_EXIT may be defined by alberta_util.h. If so, undefine it.
#ifdef DEBUG_TEST_EXIT
#undef DEBUG_TEST_EXIT
#endif // #ifdef DEBUG_TEST_EXIT

// Macro INFO may be defined by alberta_util.h. If so, undefine it.
#ifdef INFO
#undef INFO
#endif // #ifdef INFO

// Macro PRINT_INFO may be defined by alberta_util.h. If so, undefine it.
#ifdef PRINT_INFO
#undef PRINT_INFO
#endif // #ifdef PRINT_INFO

// Macro PRINT_INT_VEC may be defined by alberta_util.h. If so, undefine it.
#ifdef PRINT_INT_VEC
#undef PRINT_INT_VEC
#endif // #ifdef PRINT_INT_VEC

// Macro PRINT_REAL_VEC may be defined by alberta_util.h. If so, undefine it.
#ifdef PRINT_REAL_VEC
#undef PRINT_REAL_VEC
#endif // #ifdef PRINT_REAL_VEC

// Macro WAIT may be defined by alberta_util.h. If so, undefine it.
#ifdef WAIT
#undef WAIT
#endif // #ifdef WAIT

// Macro WAIT_REALLY may be defined by alberta_util.h. If so, undefine it.
#ifdef WAIT_REALLY
#undef WAIT_REALLY
#endif // #ifdef WAIT_REALLY

// Macro GET_WORKSPACE may be defined by alberta_util.h. If so, undefine it.
#ifdef GET_WORKSPACE
#undef GET_WORKSPACE
#endif // #ifdef GET_WORKSPACE

// Macro FREE_WORKSPACE may be defined by alberta_util.h. If so, undefine it.
#ifdef FREE_WORKSPACE
#undef FREE_WORKSPACE
#endif // #ifdef FREE_WORKSPACE

// Macro MAT_ALLOC may be defined by alberta_util.h. If so, undefine it.
#ifdef MAT_ALLOC
#undef MAT_ALLOC
#endif // #ifdef MAT_ALLOC

// Macro MAT_FREE may be defined by alberta_util.h. If so, undefine it.
#ifdef MAT_FREE
#undef MAT_FREE
#endif // #ifdef MAT_FREE

// Macro NAME may be defined by alberta_util.h. If so, undefine it.
#ifdef NAME
#undef NAME
#endif // #ifdef NAME

// Macro GET_STRUCT may be defined by alberta_util.h. If so, undefine it.
#ifdef GET_STRUCT
#undef GET_STRUCT
#endif // #ifdef GET_STRUCT

// Macro ADD_PARAMETER may be defined by alberta_util.h. If so, undefine it.
#ifdef ADD_PARAMETER
#undef ADD_PARAMETER
#endif // #ifdef ADD_PARAMETER

// Macro GET_PARAMETER may be defined by alberta_util.h. If so, undefine it.
#ifdef GET_PARAMETER
#undef GET_PARAMETER
#endif // #ifdef GET_PARAMETER

#define _ALBERTA_H_

#endif // HAVE_ALBERTA || DOXYGEN

#endif
