// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_GRID_HH
#define DUNE_GRID_COMMON_GRID_HH

/** \file
    \brief Different resources needed by all grid implementations
 */
// system includes
#include <iostream>
#include <string>
#include <vector>

// dune-common includes
#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/typeutilities.hh>

// dune-geometry includes
#include <dune/geometry/type.hh>

// local includes
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/exceptions.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridview.hh>
#include <dune/grid/common/defaultgridview.hh>
#include <dune/grid/common/entityseed.hh>

// include this file after all other, because other files might undef the
// macros that are defined in that file
#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune {

  /**
         @addtogroup Grid Grid

         The Dune Grid module defines a general interface to a parallel, in general
     nonconforming, locally refined and hierarchical
         finite element mesh.  The interface is independent of dimension and
         element type.

         @section Grid1 Terminology


         @subsection subs1 Entity

         An entity is a geometric object that is part of a grid. It is
     generalized polytope that has the same dimensionality as the grid
     or a lower dimension.


         @subsection subs20 Dimension

         A grid has a fixed dimension \f$d\f$ which is the number of coordinates
     required to specify any point in the grid. The dimension is a template parameter
     of a grid.


         @subsection subs21 Codimension of an entity

         Each entity has a codimension \f$c\f$ where \f$0 \leq c \leq d\f$ (the dimension of the grid).
         An entity with codimension \f$ c\f$ in a grid of dimension \f$ d\f$ is a \f$d-c\f$-dimensional
         object.


         @subsection subs5 Subentity

         Entities are hierarchically constructed in the sense that entities of
         codimension 0 are made up of entities of codimension 1 which are themselves
         made up of entities of codimension 2 etc. until entities of codimension \f$d-1\f$
         which consist of entities of codimension \f$ d\f$.


         @subsection subs3 Element

         An element is an entity of codimension 0.


         @subsection subs4 Vertex

         A vertex is an entity of codimension \f$ d\f$ (the same as the grid's dimension).


         @subsection subs22 World dimension

         Each grid has a world dimension \f$ w\f$ with \f$ w\geq d\f$. This is the number
     of coordinates of the positions of the grid's vertices.


         @subsection subs33 Hierarchical grid

         The %Dune grid interface describes not only a single grid but a sequence of
         grids with different resolution. This is achieved by beginning with an
         intentionally coarse grid, the so-called macro grid. Then each
     element may be individually subdivided to yield new (smaller) elements.
         This construction is recursive such that each macro element and
         all the elements that resulted from subdividing it form a tree structure.

         @subsection subs33333 Grid refinement

         The grid can only be modified in special phases, the so-called refinement phase.
         In between refinement phases the entities of the grid can not be modified in any way.
         During refinement currently only the hierachic subdivision can be modified.


         @subsection subs3333 Grid level

         All elements of the macro grid form level 0 of the grid structure. All
         elements that are obtained from an \f$ l\f$-fold subdivision of a macro
         element form level \f$ l\f$ of the grid structure.

         @subsection subs333 Leaf grid

         All elements of a grid that are not subdivided any further make up
         the leaf grid. The leaf grid is the mesh with the finest resolution.

         @subsection subs6 Assignable

         A type is said to be assignable if it has a (public) copy constructor and
         assignment operator. Note that this definition requires always both methods.


         @subsection subs7 Default-constructible

         A type is said to be default-constructible if it has a constructor without arguments.


         @subsection subs8 Copy-constructible from type X

         A type is said to be copy constructible from some other type X if it has
         a copy constructor that takes a reference to an object of type X.


         @subsection subs9 Equality-comparable

         A type is said to be equality-comparable if it has an operator==.


         @subsection subs10 LessThan-comparable

         A type is lessthan-comparable if it has an operator<.


         @subsection subs11 Dereferenceable

         A type is dereferenceable if it has an operator* that delivers
         a reference to a value type.

         @subsection subs11 Iterator

         An iterator is a type that can be dereferenced to yield an object of
         its value type, i.e. it behaves like a pointer, and that can be incremented to
         point to the next element in a linear sequence. In that respect it is comparable to
         ForwardIterator in the Standard Template Library.


         @subsection subs12 Mutable iterator

         An iterator is called mutable if the value it refers to can be changed, i.e. it is
         assignable.


         @subsection subs13 Immutable iterator

         An iterator is called immutable if the value referenced by the iterator can not
         be changed, i. e. the value is not assignable and only methods marked const on the value
         can be called.


         @subsection subs14 Model

         A type M is called a model of another type X if it implements all the methods
         of X with the intended semantics. Typically X is a type that describes an interface.


         @section Grid3 Types common to all grid implementations

         - Dune::ReferenceElement describes the topology and geometry of standard entities.
         Any given entity of the grid can be completely specified by a reference element
         and a map from this reference element to world coordinate space.

         - Dune::GeometryType defines names for the reference elements.

         - Dune::Communication defines an interface to global communication
         operations in a portable and transparent way. In particular also for sequential grids.



         @section Grid2 Types making up a grid implementation

     Each implementation of the Dune grid interface consist of a number of related types which
         together form a model of the grid interface. These types are the following:

         - %Grid which is a model of Dune::Grid where the template parameters are at least the
         dimension and the world dimension. It is a container of entities that allows to access
         these entities and that knows the number of entities.

         - %Entity which is a model of Dune::Entity. This class is parametrized by dimension and
         codimension. The entity encapsulates the topological part of an entity, i.e. its hierarchical
         construction from subentities and the relation to other entities. Entities cannot
         be created, copied or modified by the user. They can only be read-accessed through
         immutable iterators.

         - %Geometry which is a model of Dune::Geometry. This class encapsulates the geometric part
         of an entity by mapping local coordinates in a reference element to world coordinates.

         - %LevelIterator which is a model of Dune::LevelIterator is an immutable iterator
     that provides access to all entities of a given codimension and level of the
         grid.
         - %LeafIterator which is a model of Dune::LeafIterator is an immutable iterator
         that provides access to all entities of a given codimension of the leaf grid.

         - %HierarchicIterator which is a model of Dune::HierarchicIterator is an immutable
         iterator that provides access to all entities of codimension 0 that resulted from subdivision
         of a given entity of codimension 0.

         - %Intersection which is a model of Dune::Intersection
         provides access an intersection of codimension 1 of two entity of codimension 0
     or one entity and the boundary. In a conforming mesh this
         is a face of an element. For two entities with a common intersection
         the %Intersection also provides information about the geometric location
         of the intersection. Furthermore it also provides information about intersections
         of an entity with the internal or external boundaries.

         - %IntersectionIterator which is a model of Dune::IntersectionIterator
         provides access to all intersections of a given entity of codimension 0.

         - %LevelIndexSet and %LeafIndexSet which are both models of Dune::IndexSet are
         used to attach any kind of user-defined data to (subsets of) entities of the grid.
         This data is supposed to be stored in one-dimensional arrays for reasons
         of efficiency.

         - %LocalIdSet and %GlobalIdSet which are both models of Dune::IdSet are used to
         save user data during a grid refinement phase and during dynamic load balancing
         in the parallel case.


         @section Grid22 Overview of basic capabilities of the types

         <TABLE>
         <TR>
         <TD>Class</TD>
         <TD>Assignable</TD>
         <TD>DefaultConstructible</TD>
         <TD>EqualityComparable</TD>
         <TD>LessThanComparable</TD>
         </TR>
         <TR>
         <TD>Grid</TD>
         <TD>no</TD>
         <TD>no</TD>
         <TD>no</TD>
         <TD>no</TD>
         </TR>
         <TR>
         <TD>Entity</TD>
         <TD>no</TD>
         <TD>no</TD>
         <TD>no</TD>
         <TD>no</TD>
         </TR>
         <TR>
         <TD>GeometryType</TD>
         <TD>yes</TD>
         <TD>yes</TD>
         <TD>yes</TD>
         <TD>yes</TD>
         </TR>
         <TR>
         <TD>Geometry</TD>
         <TD>no</TD>
         <TD>no</TD>
         <TD>no</TD>
         <TD>no</TD>
         </TR>
         <TR>
         <TD>LevelIterator</TD>
         <TD>yes</TD>
         <TD>no</TD>
         <TD>yes</TD>
         <TD>no</TD>
         </TR>
         <TR>
         <TD>LeafIterator</TD>
         <TD>yes</TD>
         <TD>no</TD>
         <TD>yes</TD>
         <TD>no</TD>
         </TR>
         <TR>
         <TD>HierarchicIterator</TD>
         <TD>yes</TD>
         <TD>no</TD>
         <TD>yes</TD>
         <TD>no</TD>
         </TR>
         <TR>
         <TD>Intersection</TD>
         <TD>yes</TD>
         <TD>no</TD>
         <TD>yes</TD>
         <TD>no</TD>
         </TR>
         <TR>
         <TD>IntersectionIterator</TD>
         <TD>yes</TD>
         <TD>no</TD>
         <TD>yes</TD>
         <TD>no</TD>
         </TR>
         <TR>
         <TD>IndexSet</TD>
         <TD>no</TD>
         <TD>no</TD>
         <TD>no</TD>
         <TD>no</TD>
         </TR>
         <TR>
         <TD>IdSet</TD>
         <TD>no</TD>
         <TD>no</TD>
         <TD>no</TD>
         <TD>no</TD>
         </TR>
         </TABLE>
   */



  // Forward Declarations
  // --------------------

  template<int mydim, int cdim, class GridImp,template<int,int,class> class GeometryImp> class Geometry;
  template< int mydim, int cdim, class GridImp > class GlobalGeometryReference;
  template< int mydim, int cdim, class GridImp > class LocalGeometryReference;
  // dim is necessary because Entity will be specialized for codim==0 _and_ codim==dim
  // EntityImp gets GridImp as 3rd template parameter to distinguish between const and mutable grid
  template<int codim, int dim, class GridImp,template<int,int,class> class EntityImp> class Entity;
  template< int codim, class Grid, class IteratorImp > class EntityIterator;
  template<class GridImp, class EntitySeedImp> class EntitySeed;
  template< class GridImp, class IntersectionImp > class Intersection;
  template< class GridImp, class IntersectionIteratorImp, class IntersectionImp > class IntersectionIterator;
  template< class GridImp, class IndexSetImp, class IndexTypeImp = unsigned int, class TypesImp = std::vector< GeometryType > > class IndexSet;
  template<class GridImp, class IdSetImp, class IdTypeImp> class IdSet;


  //************************************************************************
  // G R I D
  //************************************************************************

  /**
     \brief Grid abstract base class
     @ingroup GIGrid

     This class is the base class for all grid implementations. Although
     no virtual functions are used we call it abstract since its
     methods do not contain an implementation but forward to the methods of
     the derived class via the Barton-Nackman trick.

     \tparam dim specifies the dimension of the grid.
     \tparam dimworld specifies the dimension of the surrounding space, this can be
       different from dim, if the grid is defined on a manifold .
     \tparam ct field type of the world vector space.
     \tparam GridFamily traits class providing all types
       associated with the grid implementation.

     \nosubgrouping
   */
  template< int dim, int dimworld, class ct, class GridFamily>
  class Grid {
    typedef typename GridFamily::Traits::Grid GridImp;
    typedef Grid<dim,dimworld,ct,GridFamily> ThisType;
  public:

    //===========================================================
    /** @name Exported constants
     */
    //@{
    //===========================================================

    //! \brief The dimension of the grid
    constexpr static int dimension = dim;

    //! \brief The dimension of the world the grid lives in.
    constexpr static int dimensionworld = dimworld;
    //@}

    //===========================================================
    /** @name Exported types
     */
    //@{
    //===========================================================

    /** \brief type of view for leaf grid */
    typedef typename GridFamily::Traits::LeafGridView LeafGridView;
    /** \brief type of view for level grid */
    typedef typename GridFamily::Traits::LevelGridView LevelGridView;


    /** \brief A Traits struct that collects all associated types of one implementation

       \tparam cd codimension. Note that not all types in this struct depend on this template parameter.
     */
    template <int cd>
    struct Codim
    {
      //! A type that is a model of Dune::Geometry<dim-cd,dimworld>.
      typedef typename GridFamily::Traits::template Codim<cd>::Geometry Geometry;

      //! A type that is a model of Dune::Geometry<dim-cd,dim>.
      typedef typename GridFamily::Traits::template Codim<cd>::LocalGeometry LocalGeometry;

      //! A type that is a model of a Dune::Entity<cd,dim,...>.
      typedef typename GridFamily::Traits::template Codim<cd>::Entity Entity;

      //! A type that is a model (not yet) of Dune::EntitySeed<cd,dim,...>.
      typedef typename GridFamily::Traits::template Codim<cd>::EntitySeed EntitySeed;

      //! A struct collecting all types depending on the partition iterator type.
      template <PartitionIteratorType pitype>
      struct Partition
      {
        /*! \brief A type that is a model of Dune::LevelIterator<cd,pitype,...>
              which is s type of iterator that may be used to examine, but not to modify, the
              entities of codimension cd with partition type
              pitype  on a certain level of the grid, i. e. the increment of
              the iterator adjusts it to the next entity on that level.
         */
        typedef typename GridFamily::Traits::template Codim<cd>::template Partition<pitype>::LevelIterator LevelIterator;
        /*! \brief A type that is a model of Dune::LeafIterator<cd,pitype,...>
              which is a type of iterator that may be used to examine, but not to modify, the
              entities of codimension cd with partition type
              pitype in the leaf grid, i. e. the increment of
              the iterator adjusts it to the next entity in the leaf grid.
         */
        typedef typename GridFamily::Traits::template Codim<cd>::template Partition<pitype>::LeafIterator LeafIterator;
      };

      /*! \brief A type that is a model of Dune::LevelIterator with partition type All_Partition
       */
      typedef typename GridFamily::Traits::template Codim<cd>::LevelIterator LevelIterator;

      /*! \brief A type that is a model of Dune::LeafIterator with partition type All_Partition
       */
      typedef typename GridFamily::Traits::template Codim<cd>::LeafIterator LeafIterator;
    };

    /*! \brief A type that is a model of Dune::Intersection, an
       intersections of two codimension 1 of two codimension 0 entities in the leaf view.
     */
    typedef typename GridFamily::Traits::LeafIntersection LeafIntersection;

    /*! \brief A type that is a model of Dune::Intersection, an
       intersections of two codimension 1 of two codimension 0 entities in a level view.
     */
    typedef typename GridFamily::Traits::LevelIntersection LevelIntersection;

    /*! \brief A type that is a model of Dune::IntersectionIterator
       which is an iterator that allows to examine, but not to modify, the
       intersections of codimension 1 of an leaf element (entity of codimension 0)
       with other leaf elements.
     */
    typedef typename GridFamily::Traits::LeafIntersectionIterator LeafIntersectionIterator;

    /*! \brief A type that is a model of Dune::IntersectionIterator
       which is an iterator that allows to examine, but not to modify, the
       intersections of codimension 1 of an element (entity of codimension 0)
       with other elements on the same level.
     */
    typedef typename GridFamily::Traits::LevelIntersectionIterator LevelIntersectionIterator;

    /*! \brief A type that is a model of Dune::HierarchicIterator
       A type of iterator that allows to examine, but not to modify, entities
       of codimension 0 that result from refinement of an entity of
       codimension 0.
     */
    typedef typename GridFamily::Traits::HierarchicIterator HierarchicIterator;

    /*!  \brief A type that is a model of Dune::IndexSet
       which provides a consecutive, but non persistent, numbering for
       entities on a grid level.
     */
    typedef typename GridFamily::Traits::LevelIndexSet LevelIndexSet;

    /*! \brief A type that is a model of Dune::IndexSet
       which provides a consecutive, but non persistent, numbering for
       entities in the leaf grid.
     */
    typedef typename GridFamily::Traits::LeafIndexSet LeafIndexSet;

    /*!  \brief A type that is a model of Dune::IdSet
       which provides a unique and persistent numbering for
       all entities in the grid. The numbering is unique over all processes
       over which the grid is partitioned. The numbering is not necessarily
       consecutive.
     */
    typedef typename GridFamily::Traits::GlobalIdSet GlobalIdSet;

    /*! \brief A type that is a model of Dune::IdSet
       which provides a unique and persistent numbering for
       all entities in the grid. The numbering is only unique in a single process
       and it is not necessarily consecutive.
     */
    typedef typename GridFamily::Traits::LocalIdSet LocalIdSet;

  protected:
    template <class T>
    using Communication_t = typename T::Communication;
    template <class T>
    using DeprecatedCollectiveCommunication_t = typename T::CollectiveCommunication;

  public:
    /*! \brief A type that is a model of Dune::Communication.
       It provides a portable way for communication on the set
       of processes used by the grid.
     */
    // if this line produces a warning then the Communication typedef is missing
    // in the GridFamily::Traits
    using Communication = detected_or_fallback_t< DeprecatedCollectiveCommunication_t,
                                                  Communication_t, typename GridFamily::Traits>;

    /** \deprecated Use Communication instead! Will be removed at some point in the future.
     */
    using CollectiveCommunication [[deprecated("CollectiveCommunication is deprecated, use Communication instead!")]] = Communication;

    //! Define type used for coordinates in grid module
    typedef ct ctype;
    //@}


    //===========================================================
    /** @name Size methods
     */
    //@{
    //===========================================================

    /*! \brief Return maximum level defined in this grid. Levels are numbered
       0 ... maxLevel with 0 the coarsest level.
     */
    int maxLevel() const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().maxLevel());
      return asImp().maxLevel();
    }

    //! Return number of grid entities of a given codim on a given level in this process
    int size (int level, int codim) const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().size(level,codim));
      return asImp().size(level,codim);
    }

    //! Return number of leaf entities of a given codim in this process
    int size (int codim) const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().size(codim));
      return asImp().size(codim);
    }

    //! Return number of entities per level and geometry type in this process
    int size (int level, GeometryType type) const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().size(level,type));
      return asImp().size(level,type);
    }

    //! Return number of leaf entities per geometry type in this process
    int size (GeometryType type) const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().size(type));
      return asImp().size(type);
    }
    //@}


    /** \brief returns the number of boundary segments within the macro grid
     *
     *  \returns number of boundary segments within the macro grid
     */
    size_t numBoundarySegments () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().numBoundarySegments());
      return asImp().numBoundarySegments();
    }

    //===========================================================
    /** @name Views
     */
    //@{
    //===========================================================

    //! View for a grid level for All_Partition
    LevelGridView levelGridView(int level) const {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().levelGridView(level)));
      return asImp().levelGridView(level);
    }

    //! View for the leaf grid for All_Partition
    LeafGridView leafGridView() const {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().leafGridView()));
      return asImp().leafGridView();
    }

    //@}


    //===========================================================
    /** @name Access to index and id sets
     */
    //@{
    //===========================================================

    //! return const reference to the grids global id set
    const GlobalIdSet &globalIdSet () const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().globalIdSet());
      return asImp().globalIdSet();
    }

    //! return const reference to the grids local id set
    const LocalIdSet &localIdSet () const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().localIdSet());
      return asImp().localIdSet();
    }

    //! return const reference to the grids level index set for level level
    const LevelIndexSet &levelIndexSet ( int level ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().levelIndexSet(level));
      return asImp().levelIndexSet(level);
    }

    //! return const reference to the grids leaf index set
    const LeafIndexSet &leafIndexSet () const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().leafIndexSet());
      return asImp().leafIndexSet();
    }
    //@}


    //===========================================================
    /** @name Adaptivity and grid refinement
     */
    //@{
    //===========================================================

    /** \brief Refine the grid refCount times using the default refinement rule.
     *
     * This behaves like marking all elements for refinement and then calling preAdapt, adapt and postAdapt.
     * The state after globalRefine is comparable to the state after postAdapt.
     */
    void globalRefine (int refCount)
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().globalRefine(refCount));
      return;
    }

    /** \brief Marks an entity to be refined/coarsened in a subsequent adapt.

       \param[in] refCount Number of subdivisions that should be applied. Negative value means coarsening.
       \param[in] e        Entity that should be marked

       \return true if Entity was marked, false otherwise.
     */
    bool mark( int refCount, const typename Codim<0>::Entity & e )
    {
      return asImp().mark(refCount,e);
    }

    /** \brief returns adaptation mark for given entity

       \param[in] e   Entity for which adaptation mark should be determined

       \return int adaptation mark currently set for given Entity e
     */
    int getMark(const typename Codim<0>::Entity & e) const
    {
      return asImp().getMark(e);
    }

    /** \brief To be called after entities have been marked and before adapt() is called.
     *
     * This sets the mightVanish flags of the elements for the next adapt call.
     *
     * \return true if an entity may be coarsened during a subsequent adapt(), false otherwise.
     */
    bool preAdapt ()
    {
      return asImp().preAdapt();
    }

    /** \brief Refine all positive marked leaf entities,
        coarsen all negative marked entities if possible

            \return true if a least one entity was refined

            The complete adaptation process works as follows:

            - mark entities with the mark() method
            - call preAdapt()
            - if preAdapt() returned true: possibly save current solution
            - call adapt()
            - if adapt() returned true: possibly interpolate the (saved) solution
            - call postAdapt()
     */
    bool adapt ()
    {
      return asImp().adapt();
    }

    /** \brief To be called after grid has been adapted and information left over by the adaptation has been processed.
     *
     * This removes the isNew flags of the elements from the last adapt call.
     */
    void postAdapt()
    {
      return asImp().postAdapt();
    }
    //@}


    //===========================================================
    /** @name Parallel data distribution and communication
     */
    //@{
    //===========================================================

    //! return const reference to a communication object. The return type is a model of Dune::Communication.
    const Communication &comm () const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().comm());
      return asImp().comm();
    }
    //@}

    /*! \brief Re-balances the load each process has to handle for a parallel grid,
     *  \return True if the grid has changed.
     */
    bool loadBalance()
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().loadBalance());
      return asImp().loadBalance();
    }

    /*! \brief Re-balances the load each process has to handle for a parallel grid and moves the data.
     * \param data A data handle telling the method which data should be communicated
     * and how. Has to adhere to the interface describe by CommDataHandleIf
     * just like the data handle for the communicate
     * methods.
     * \return True if the grid has changed.
     */
    template<class DataHandle>
    bool loadBalance (DataHandle& data)
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().loadBalance(data));
      return asImp().loadBalance(data);
    }

    /** \brief obtain Entity from EntitySeed. */
    template < class EntitySeed >
    typename Codim< EntitySeed :: codimension > :: Entity
    entity( const EntitySeed& seed ) const
    {
      //CHECK_INTERFACE_IMPLEMENTATION( asImp().entity( seed ) );
      return asImp().entity( seed );
    }
  protected:
    //!  Barton-Nackman trick
    GridImp& asImp () {return static_cast<GridImp &> (*this);}
    //!  Barton-Nackman trick
    const GridImp& asImp () const {return static_cast<const GridImp &>(*this);}
  };

#undef CHECK_INTERFACE_IMPLEMENTATION
#undef CHECK_AND_CALL_INTERFACE_IMPLEMENTATION

  /**
   * \addtogroup GIGridView
   *  @{
   */

  /**
   * \brief level grid view for the given grid and level.
   *
   * Identical to the method in the Grid interface, but provided
   * as a free function.
   *
   * \see Grid::levelGridView
   *
   * \param grid Grid to obtain the level grid view for
   * \param level level of the level grid view
   */
  template<int dim, int dimworld, class ct, class GridFamily>
  typename Grid<dim, dimworld, ct, GridFamily>::LevelGridView
  levelGridView(const Grid<dim, dimworld, ct, GridFamily>& grid, int level)
  {
    return grid.levelGridView(level);
  }

  /**
   * \brief leaf grid view for the given grid
   *
   * Identical to the method in the Grid interface, but provided
   * as a free function.
   *
   * \see Grid::leafGridView
   *
   * \param grid Grid to obtain the leaf grid view for
   */
  template<int dim, int dimworld, class ct, class GridFamily>
  typename Grid<dim, dimworld, ct, GridFamily>::LeafGridView
  leafGridView(const Grid<dim, dimworld, ct, GridFamily>& grid)
  {
    return grid.leafGridView();
  }

  /** @} */

  /**
     \ingroup GridDevel
     @{

     A Grid is a container of grid entities. Given a dimension dim
     these entities have a codimension codim with 0 <= codim <= dim.

     The Grid is assumed to be hierachically refined and nested. It
     enables iteration over entities of a given level and codimension.

     The grid can be non-matching.

     All information is provided to allocate degrees of freedom in
     appropriate vector data structures (which are not part of this
     module).

     Template class Grid defines a "base class" for all grids.

     \par Classes implementing the Grid Interface
     \li Dune::AlbertaGrid <br>
         <i> Provides the simplicial meshes of the finite element tool box
             ALBERTA (http://www.alberta-fem.de/)
             written by Kunibert Siebert and Alfred Schmidt.</i>
     \li Dune::OneDGrid <br>
         <i> Onedimensional adaptive grid</i>
     \li Dune::UGGrid <br>
         <i> Provides the meshes of the finite element toolbox UG3.
             (described in https://doi.org/10.1007/s007910050003).</i>
     \li Dune::YaspGrid (Yet Another Structured Parallel Grid) <br>
         <i> Provides a distributed structured cube mesh.</i>
   */
  template<int dim,
      int dimworld,
      class ct,
      class GridFamily>
  class GridDefaultImplementation : public Grid <dim,dimworld,ct,GridFamily>
  {
    typedef typename GridFamily::Traits::Grid GridImp;

  public:
    /**
     * \brief The traits of this class.
     *
     * Presents the typedefs as described in GridTraits.
     */
    typedef typename GridFamily::Traits Traits;

    //! View for a grid level for All_Partition
    typename Traits::LevelGridView levelGridView(int level) const
    {
      typedef typename Traits::LevelGridView View;
      typedef typename View::GridViewImp ViewImp;
      return View(ViewImp(asImp(),level));
    }

    //! View for the leaf grid for All_Partition
    typename Traits::LeafGridView leafGridView() const
    {
      typedef typename Traits::LeafGridView View;
      typedef typename View::GridViewImp ViewImp;
      return View(ViewImp(asImp()));
    }

    //***************************************************************
    //  Interface for Adaptation
    //***************************************************************

    /** \brief Marks an entity to be refined/coarsened in a subsequent adapt.

       \param[in] refCount Number of subdivisions that should be applied. Negative value means coarsening.
       \param[in] e        Entity to Entity that should be refined

       \return true if Entity was marked, false otherwise.

       \note
          -  \b default \b implementation is: return false; for grids with no
             adaptation.
          -  for the grid programmer:
             this method is implemented as a template method, because the
             Entity type is not defined when the class is instantiated
             You won't need this trick in the implementation.
             In your implementation you should use it as
             \code
             bool mark( int refCount,
                        typename Traits::template Codim<0>::Entity & e ).
             \endcode
             This template method will vanish due to the inheritance
             rules.
     */
    bool mark( int refCount, const typename Traits :: template Codim<0>::Entity & e )
    {
      return false;
    }

    /** \brief returns adaptation mark for given entity, i.e. here the
     * default implementation returns 0.

       \param[in] e   Entity for which adaptation mark should be determined

       \return int adaptation mark, here the default value 0 is returned
     */
    int getMark ( const typename Traits::template Codim< 0 >::Entity &e ) const
    {
      return 0;
    }

    /** \brief Refine all positive marked leaf entities
        coarsen all negative marked entities if possible
        \return true if a least one entity was refined

        \note this default implementation always returns false
          so grid with no adaptation doesn't need to implement these methods
     */
    bool adapt ()    { return false; }

    //! returns true, if at least one entity is marked for adaption
    bool preAdapt () { return false; }

    //! clean up some markers
    void postAdapt() {}

    /*! \brief default implementation of load balance does nothing and returns false */
    bool loadBalance()
    {
      return false;
    }

    /*! \brief default implementation of load balance does nothing and returns false */
    template<class DataHandle>
    bool loadBalance ([[maybe_unused]] DataHandle& data)
    {
      return false;
    }

  protected:
    using Grid< dim, dimworld, ct, GridFamily >::asImp;
  };

  /** @} */


  /**
     \brief A traits struct that collects all associated types of one grid model
     @ingroup GIMiscellaneous

     \tparam dim Grid dimension
     \tparam dimw Dimension of the world that the grid is embedded in
     \tparam GIDType Type used for global ids
     \tparam LIDType Type used for local ids
     \tparam CCType Communication implementation class
   */
  template <int dim, int dimw, class GridImp,
      template<int,int,class> class GeometryImp,
      template<int,int,class> class EntityImp,
      template<int,PartitionIteratorType,class> class LevelIteratorImp,
      template<class> class LeafIntersectionImp,
      template<class> class LevelIntersectionImp,
      template<class> class LeafIntersectionIteratorImp,
      template<class> class LevelIntersectionIteratorImp,
      template<class> class HierarchicIteratorImp,
      template<int,PartitionIteratorType,class> class LeafIteratorImp,
      class LevelIndexSetImp, class LeafIndexSetImp,
      class GlobalIdSetImp, class GIDType, class LocalIdSetImp, class LIDType, class CCType,
      template<class> class LevelGridViewTraits,
      template<class> class LeafGridViewTraits,
      template<int,class> class EntitySeedImp,
      template<int,int,class> class LocalGeometryImp = GeometryImp
      >
  struct GridTraits
  {
    /** \brief The type that implements the grid. */
    typedef GridImp Grid;

    /** \brief The type of the intersection at the leafs of the grid. */
    typedef Dune::Intersection< const GridImp, LeafIntersectionImp< const GridImp > >  LeafIntersection;
    /** \brief The type of the intersection at the levels of the grid. */
    typedef Dune::Intersection< const GridImp, LevelIntersectionImp< const GridImp > > LevelIntersection;
    /** \brief The type of the intersection iterator at the leafs of the grid. */
    typedef Dune::IntersectionIterator< const GridImp, LeafIntersectionIteratorImp< const GridImp >, LeafIntersectionImp< const GridImp > > LeafIntersectionIterator;
    /** \brief The type of the intersection iterator at the levels of the grid. */
    typedef Dune::IntersectionIterator< const GridImp, LevelIntersectionIteratorImp< const GridImp >, LevelIntersectionImp< const GridImp > > LevelIntersectionIterator;

    /** \brief The type of the  hierarchic iterator. */
    typedef Dune::EntityIterator< 0, const GridImp, HierarchicIteratorImp< const GridImp > > HierarchicIterator;

    /**
     * \brief Traits associated with a specific codim.
     * \tparam cd The codimension.
     */
    template <int cd>
    struct Codim
    {
    public:
      typedef GeometryImp<dim-cd, dimw, const GridImp> GeometryImpl;
      typedef LocalGeometryImp<dim-cd, dim, const GridImp> LocalGeometryImpl;
      //! IMPORTANT: Codim<codim>::Geometry == Geometry<dim-codim,dimw>
      /** \brief The type of the geometry associated with the entity.*/
      typedef Dune::Geometry<dim-cd, dimw, const GridImp, GeometryImp> Geometry;
      /** \brief The type of the local geometry associated with the entity.*/
      typedef Dune::Geometry<dim-cd, dim, const GridImp, LocalGeometryImp> LocalGeometry;
      /** \brief The type of the entity. */
      // we could - if needed - introduce another struct for dimglobal of Geometry
      typedef Dune::Entity<cd, dim, const GridImp, EntityImp> Entity;

      /** \brief The type of the entity seed of this codim.*/
      typedef Dune::EntitySeed<const GridImp, EntitySeedImp<cd, const GridImp> > EntitySeed;

      /**
       * \brief Traits associated with a specific grid partition type.
       * \tparam pitype The type of the grid partition.
       */
      template <PartitionIteratorType pitype>
      struct Partition
      {
        /** \brief The type of the iterator over the level entities of this codim on this partition. */
        typedef Dune::EntityIterator< cd, const GridImp, LevelIteratorImp< cd, pitype, const GridImp > > LevelIterator;
        /** \brief The type of the iterator over the leaf entities of this codim on this partition. */
        typedef Dune::EntityIterator< cd, const GridImp, LeafIteratorImp< cd, pitype, const GridImp > > LeafIterator;
      };

      /** \brief The type of the iterator over all leaf entities of this codim. */
      typedef typename Partition< All_Partition >::LeafIterator LeafIterator;

      /** \brief The type of the entity pointer for entities of this codim.*/
      typedef typename Partition< All_Partition >::LevelIterator LevelIterator;

    private:
      friend class Dune::Entity<cd, dim, const GridImp, EntityImp>;
    };

    /** \brief type of view for leaf grid */
    typedef Dune::GridView< LeafGridViewTraits< const GridImp > > LeafGridView;
    /** \brief type of view for level grid */
    typedef Dune::GridView< LevelGridViewTraits< const GridImp > > LevelGridView;

    /** \brief The type of the level index set. */
    typedef IndexSet<const GridImp,LevelIndexSetImp> LevelIndexSet;
    /** \brief The type of the leaf index set. */
    typedef IndexSet<const GridImp,LeafIndexSetImp> LeafIndexSet;
    /** \brief The type of the global id set. */
    typedef IdSet<const GridImp,GlobalIdSetImp,GIDType> GlobalIdSet;
    /** \brief The type of the local id set. */
    typedef IdSet<const GridImp,LocalIdSetImp,LIDType> LocalIdSet;

    /** \brief The type of the communication. */
    typedef CCType Communication;

    /** \deprecated Use Communication instead! Will be removed at some point in the future. */
    [[deprecated("Use Communication instead!")]]
    typedef Communication CollectiveCommunication;
  };

  // Definition of capabilities for the interface class
  namespace Capabilities
  {

    // capabilities for the interface class depend on the implementation
    template< int dim, int dimworld, typename ct, class GridFamily , int codim >
    struct hasEntity< Grid< dim, dimworld, ct, GridFamily >, codim >
    {
      static const bool v = hasEntity< typename GridFamily::Traits::Grid, codim >::v;
    };

    // capabilities for the interface class depend on the implementation
    template< int dim, int dimworld, typename ct, class GridFamily , int cdim >
    struct hasEntity< GridDefaultImplementation<dim,dimworld,ct,GridFamily>, cdim >
    {
      typedef GridDefaultImplementation<dim,dimworld,ct,GridFamily> GridType;
      typedef typename GridType::Traits::Grid GridImp;
      static const bool v = hasEntity<GridImp,cdim>::v;
    };

  } // end namespace Capabilities

  //! for creation of an engine interface object like Entity or Geometry
  //! one has to derive a class to create the object because the
  //! contructors of the interface object classes are protected
  //! therefore here a generic implementation for this object creation is
  //! provided
  template <class InterfaceType>
  struct MakeableInterfaceObject : public InterfaceType
  {
    typedef typename InterfaceType::Implementation ImplementationType;
    //! create interface object by calling the contructor of the base class
    explicit MakeableInterfaceObject ( const ImplementationType &realImp )
      : InterfaceType( realImp )
    {}
  };
}

#include "geometry.hh"
#include "entity.hh"
#include "intersection.hh"
#include "intersectioniterator.hh"
#include "entityiterator.hh"
#include "indexidset.hh"

#endif // #ifndef DUNE_GRID_COMMON_GRID_HH
