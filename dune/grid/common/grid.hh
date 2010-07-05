// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_HH
#define DUNE_GRID_HH

/** \file
    \brief Different resources needed by all grid implementations
 */
// system includes
#include <iostream>
#include <string>

// dune-common includes
#include <dune/common/deprecated.hh>
#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/geometrytype.hh>

// local includes
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/exceptions.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridview.hh>
#include <dune/grid/common/defaultgridview.hh>

// inlcude this file after all other, because other files might undef the
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

         - Dune::CollectiveCommunication defines an interface to global communication
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

         - %EntityPointer which is a model of Dune::EntityPointer. This is a dereferenceable
         type that delivers a reference to an entity. Moreover it is immutable, i.e. the
         referenced entity can not be modified.

         - %LevelIterator which is a model of Dune::LevelIterator is an immutable iterator
     that provides access to all entities of a given codimension and level of the
         grid. %EntityPointer is copy-constructible from a %LevelIterator.

         - %LeafIterator which is a model of Dune::LeafIterator is an immutable iterator
         that provides access to all entities of a given codimension of the leaf grid.
         %EntityPointer is copy-constructible from a %LeafIterator.

         - %HierarchicIterator which is a model of Dune::HierarchicIterator is an immutable
         iterator that provides access to all entities of codimension 0 that resulted from subdivision
         of a given entity of codimension 0. %EntityPointer is copy-constructible from a
         %HierarchicIterator.

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
         <TD>EntityPointer</TD>
         <TD>yes</TD>
         <TD>no</TD>
         <TD>yes</TD>
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
  // dim is necessary because Entity will be specialized for codim==0 _and_ codim==dim
  // EntityImp gets GridImp as 3rd template parameter to distinguish between const and mutable grid
  template<int codim, int dim, class GridImp,template<int,int,class> class EntityImp> class Entity;
  template<class GridImp, class EntityPointerImp> class EntityPointer;
  template<int codim, PartitionIteratorType pitype, class GridImp,
      template<int,PartitionIteratorType,class> class LevelIteratorImp> class LevelIterator;
  template<class GridImp, template<class> class IntersectionImp> class Intersection;
  template<class GridImp, template<class> class IntersectionIteratorImp, template<class> class IntersectionImp> class IntersectionIterator;
  template<class GridImp, template<class> class HierarchicIteratorImp> class HierarchicIterator;
  template<int codim, PartitionIteratorType pitype, class GridImp,
      template<int,PartitionIteratorType,class> class LeafIteratorImp> class LeafIterator;
  template<class GridImp> class GenericLeafIterator;
  template<class GridImp, class IndexSetImp, class IndexTypeImp=unsigned int> class IndexSet;
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

    //! A constant that exports the template parameter dim
    enum {
      //! \brief The dimension of the grid
      dimension=dim
    };

    //! A constant that exports the template parameter dimworld
    enum {
      //! \brief The dimension of the world the grid lives in.
      dimensionworld=dimworld
    };
    //@}

    //===========================================================
    /** @name Exported types
     */
    //@{
    //===========================================================

    /** \brief Types for GridView */
    template <PartitionIteratorType pitype>
    struct Partition
    {
      typedef typename GridFamily::Traits::template Partition<pitype>::LevelGridView
      LevelGridView;
      typedef typename GridFamily::Traits::template Partition<pitype>::LeafGridView
      LeafGridView;
    };
    /** \brief View types for All_Partition */
    typedef typename Partition< All_Partition > :: LevelGridView LevelGridView;
    typedef typename Partition< All_Partition > :: LeafGridView LeafGridView;


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

      //! A type that is a model of Dune::EntityPointer<cd,dim,...>.
      typedef typename GridFamily::Traits::template Codim<cd>::EntityPointer EntityPointer;

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

    /*! \brief A type that is a model of Dune::LeafIntersection, an
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

    /*! \brief A type that is a model of Dune::CollectiveCommunication.
       It provides a portable way for collective communication on the set
       of processes used by the grid.
     */
    typedef typename GridFamily::Traits::CollectiveCommunication CollectiveCommunication;

    //! Define type used for coordinates in grid module
    typedef ct ctype;
    //@}


    //===========================================================
    /** @name Grid id
     */
    //@{
    //===========================================================

    //! Return the id of the grid
    std::string name() const DUNE_DEPRECATED
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().name());
      return asImp().name();
    }

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

    //! View for a grid level
    template<PartitionIteratorType pitype>
    typename Partition<pitype>::LevelGridView levelView(int level) const {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().template levelView<pitype>(level)));
      return asImp().template levelView<pitype>(level);
    }

    //! View for the leaf grid
    template<PartitionIteratorType pitype>
    typename Partition<pitype>::LeafGridView leafView() const {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().template leafView<pitype>()));
      return asImp().template leafView<pitype>();
    }

    //! View for a grid level for All_Partition
    LevelGridView levelView(int level) const {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().levelView(level)));
      return asImp().levelView(level);
    }

    //! View for the leaf grid for All_Partition
    LeafGridView leafView() const {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().leafView()));
      return asImp().leafView();
    }

    //@}


    //===========================================================
    /** @name Iterators
     */
    //@{
    //===========================================================

    //! Iterator to first entity of given codim on level
    template<int cd, PartitionIteratorType pitype>
    typename Codim<cd>::template Partition<pitype>::LevelIterator lbegin (int level) const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().template lbegin<cd,pitype>(level)));
      return asImp().template lbegin<cd,pitype>(level);
    }

    //! one past the end on this level
    template<int cd, PartitionIteratorType pitype>
    typename Codim<cd>::template Partition<pitype>::LevelIterator lend (int level) const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().template lend<cd,pitype>(level)));
      return asImp().template lend<cd,pitype>(level);
    }

    //! Iterator to first entity of given codim on level for PartitionType All_Partition
    template<int cd>
    typename Codim<cd>::template Partition<All_Partition>::LevelIterator lbegin (int level) const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().template lbegin<cd>(level)));
      return asImp().template lbegin<cd>(level);
    }

    //! one past the end on this level for PartitionType All_Partition
    template<int cd>
    typename Codim<cd>::template Partition<All_Partition>::LevelIterator lend (int level) const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().template lend<cd>(level)));
      return asImp().template lend<cd>(level);
    }

    //! Iterator to first entity of given codim on leaf grid
    template<int cd, PartitionIteratorType pitype>
    typename Codim<cd>::template Partition<pitype>::LeafIterator leafbegin () const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().template leafbegin<cd,pitype>()));
      return asImp().template leafbegin<cd,pitype>();
    }

    //! one past the end on the leaf level grid
    template<int cd, PartitionIteratorType pitype>
    typename Codim<cd>::template Partition<pitype>::LeafIterator leafend () const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().template leafend<cd,pitype>()));
      return asImp().template leafend<cd,pitype>();
    }

    //! Iterator to first entity of given codim on leaf grid for PartitionType All_Partition
    template<int cd>
    typename Codim<cd>::template Partition<All_Partition>::LeafIterator leafbegin () const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().template leafbegin<cd,All_Partition>()));
      return asImp().template leafbegin<cd,All_Partition>();
    }

    //! one past the end on the leaf grid for PartitionType All_Partition
    template<int cd>
    typename Codim<cd>::template Partition<All_Partition>::LeafIterator leafend () const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().template leafend<cd,All_Partition>()));
      return asImp().template leafend<cd,All_Partition>();
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

    //! Return size of overlap for a given codim on a given level
    int overlapSize (int level, int codim) const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().overlapSize(level,codim));
      return asImp().overlapSize(level,codim);
    }

    //! Return size of overlap region for a given codim on the leaf grid
    int overlapSize (int codim) const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().overlapSize(codim));
      return asImp().overlapSize(codim);
    }

    //! Return size of ghost region for a given codim on a given level
    int ghostSize (int level, int codim) const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().ghostSize(level,codim));
      return asImp().ghostSize(level,codim);
    }

    //! Return size of ghost region for a given codim on the leaf grid
    int ghostSize (int codim) const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().ghostSize(codim));
      return asImp().ghostSize(codim);
    }

    /*! \brief Communicate information on distributed entities on a given level
          Template parameter is a model of Dune::CommDataHandleIF
     */
    template<class DataHandleImp, class DataTypeImp>
    void communicate (CommDataHandleIF<DataHandleImp,DataTypeImp> & data, InterfaceType iftype, CommunicationDirection dir, int level) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION((asImp().template communicate<DataHandleImp,DataTypeImp>(data,iftype,dir,level)));
      return;
    }

    /*! \brief Communicate information on distributed entities on the leaf grid
          Template parameter is a model of Dune::CommDataHandleIF
     */
    template<class DataHandleImp, class DataTypeImp>
    void communicate (CommDataHandleIF<DataHandleImp,DataTypeImp> & data, InterfaceType iftype, CommunicationDirection dir) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION((asImp().template communicate<DataHandleImp,DataTypeImp>(data,iftype,dir)));
      return;
    }

    //! return const reference to a collective communication object. The return type is a model of Dune::CollectiveCommunication.
    const CollectiveCommunication &comm () const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().comm());
      return asImp().comm();
    }
    //@}

    /*! \brief Re-balances the load each process has to handle for a parallel grid,
     *  if grid has changed , true is returned
     */
    bool loadBalance()
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().loadBalance());
      return asImp().loadBalance();
    }

    /*! \brief Re-balances the load each process has to handle for a parallel grid,
     * the DataHandle data works like the data handle for the communicate
     * methods. If grid has changed , true is returned.
     */
    template<class DataHandle>
    bool loadBalance (DataHandle& data)
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().loadBalance(data));
      return asImp().loadBalance(data);
    }

  protected:
    //!  Barton-Nackman trick
    GridImp& asImp () {return static_cast<GridImp &> (*this);}
    //!  Barton-Nackman trick
    const GridImp& asImp () const {return static_cast<const GridImp &>(*this);}
  };

#undef CHECK_INTERFACE_IMPLEMENTATION
#undef CHECK_AND_CALL_INTERFACE_IMPLEMENTATION


  //************************************************************************
  //
  //  Default Methods of Grid
  //
  //************************************************************************
  //

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
     \li Dune::ALUSimplexGrid and Dune::ALUCubeGrid and ALUConformGrid <br>
         <i> 2d/3D grid with support for non-conform adaptation and dynamic load balancing </i>
     \li Dune::OneDGrid <br>
         <i> Onedimensional adaptive grid</i>
     \li Dune::SGrid <br>
         <i> A structured mesh in d dimensions consisting of "cubes".</i>
     \li Dune::UGGrid <br>
         <i> Provides the meshes of the finite element toolbox UG.
             (http://sit.iwr.uni-heidelberg.de/~ug).</i>
     \li Dune::YaspGrid (Yet Another Structured Parallel Grid) <br>
         <i> Provides a distributed structured cube mesh.</i>

     For installation instructions for external grid managers see http://www.dune-project.org/external_libraries/index.html .

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

    //! View for a grid level
    template<PartitionIteratorType pitype>
    typename Traits::template Partition<pitype>::LevelGridView
    levelView(int level) const {
      typedef typename Traits::template Partition<pitype>::LevelGridView View;
      typedef typename View::GridViewImp ViewImp;
      return View(ViewImp(asImp(),level));
    }

    //! View for the leaf grid
    template<PartitionIteratorType pitype>
    typename Traits::template Partition<pitype>::LeafGridView leafView() const {
      typedef typename Traits::template Partition<pitype>::LeafGridView View;
      typedef typename View::GridViewImp ViewImp;
      return View(ViewImp(asImp()));
    }

    //! View for a grid level for All_Partition
    typename Traits::template Partition<All_Partition>::LevelGridView
    levelView(int level) const {
      typedef typename Traits::template Partition<All_Partition>::LevelGridView View;
      typedef typename View::GridViewImp ViewImp;
      return View(ViewImp(asImp(),level));
    }

    //! View for the leaf grid for All_Partition
    typename Traits::template Partition<All_Partition>::LeafGridView
    leafView() const {
      typedef typename Traits::template Partition<All_Partition>::LeafGridView View;
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

        - Note: this default implementation always returns false
          so grid with no adaptation doesn't need to implement these methods
     */
    bool adapt ()    { return false; }

    //! returns true, if at least one entity is marked for adaption
    bool preAdapt () { return false; }

    //! clean up some markers
    void postAdapt() {}

    /** \brief ghostSize is zero by default */
    int ghostSize (int level, int codim) const { return 0; }

    /** \brief overlapSize is zero by default */
    int overlapSize (int level, int codim) const { return 0; }

    /** \brief ghostSize is zero by default */
    int ghostSize (int codim) const { return 0; }

    /** \brief overlapSize is zero by default */
    int overlapSize (int codim) const { return 0; }

    /** dummy communicate, doing nothing  */
    template<class DataHandleImp, class DataTypeImp>
    void communicate (CommDataHandleIF<DataHandleImp,DataTypeImp> & data,
                      InterfaceType iftype, CommunicationDirection dir, int level) const
    {}

    /** dummy communicate, doing nothing  */
    template<class DataHandleImp, class DataTypeImp>
    void communicate (CommDataHandleIF<DataHandleImp,DataTypeImp> & data,
                      InterfaceType iftype, CommunicationDirection dir) const
    {}

    /*! \brief default implementation of load balance does nothing and returns false */
    bool loadBalance()
    {
      return false;
    }

    /*! \brief default implementation of load balance does nothing and returns false */
    template<class DataHandle>
    bool loadBalance (DataHandle& data)
    {
      return false;
    }
  protected:
    /**
     * @brief Helper class to choose correct implementation return type for getRealImplementation
     *
     * If the template parameter is const, const typename T::ImplementationType is returned otherwise
     * just typename ::ImplementationType.
     */
    template<class T>
    class ReturnImplementationType : public T // implement friendship via subclassing
    {
    public:
      /** @brief The correct type of the implementation to return. */
      typedef typename T::ImplementationType ImplementationType;
    private:
      // constructor in only need to compile
      ReturnImplementationType(const T& t) : T(t) {}
    };

    template<class T>
    class ReturnImplementationType<const T> : public T // implement friendship via subclassing
    {
    public:
      typedef const typename T::ImplementationType ImplementationType;
    private:
      // constructor in only need to compile
      ReturnImplementationType(const T& t) : T(t) {}
    };

    //! return real implementation of interface class
    template <class InterfaceType>
    static typename ReturnImplementationType<InterfaceType>::ImplementationType &
    getRealImplementation (InterfaceType &i) { return i.getRealImp(); }

  protected:
    using Grid< dim, dimworld, ct, GridFamily >::asImp;
  };

  /** @} */


  /**
     \brief A traits struct that collects all associated types of one grid model
     @ingroup GIMiscellaneous


     Template parameters:

     - <tt>dim</tt>
   */
  template <int dim, int dimw, class GridImp,
      template<int,int,class> class GeometryImp,
      template<int,int,class> class EntityImp,
      template<int,class> class EntityPointerImp,
      template<int,PartitionIteratorType,class> class LevelIteratorImp,
      template<class> class LeafIntersectionImp,
      template<class> class LevelIntersectionImp,
      template<class> class LeafIntersectionIteratorImp,
      template<class> class LevelIntersectionIteratorImp,
      template<class> class HierarchicIteratorImp,
      template<int,PartitionIteratorType,class> class LeafIteratorImp,
      class LevelIndexSetImp, class LeafIndexSetImp,
      class GlobalIdSetImp, class GIDType, class LocalIdSetImp, class LIDType, class CCType,
      template<class,PartitionIteratorType> class LevelGridViewTraits = DefaultLevelGridViewTraits,
      template<class,PartitionIteratorType> class LeafGridViewTraits = DefaultLeafGridViewTraits
      >
  struct GridTraits
  {
    /** \brief The type that implementing the grid. */
    typedef GridImp Grid;

    /** \brief The type of the intersection at the leafs of the grid. */
    typedef Dune::Intersection<const GridImp, LeafIntersectionImp>  LeafIntersection;
    /** \brief The type of the intersection at the levels of the grid. */
    typedef Dune::Intersection<const GridImp, LevelIntersectionImp> LevelIntersection;
    /** \brief The type of the intersection iterator at the leafs of the grid. */
    typedef Dune::IntersectionIterator<const GridImp, LeafIntersectionIteratorImp, LeafIntersectionImp>   LeafIntersectionIterator;
    /** \brief The type of the intersection iterator at the levels of the grid. */
    typedef Dune::IntersectionIterator<const GridImp, LevelIntersectionIteratorImp, LevelIntersectionImp> LevelIntersectionIterator;

    /** \brief The type of the  hierarchic iterator. */
    typedef Dune::HierarchicIterator<const GridImp, HierarchicIteratorImp> HierarchicIterator;

    /**
     * \brief Traits associated with a specific codim.
     * \tparam cd The codimension.
     */
    template <int cd>
    struct Codim
    {
      //! IMPORTANT: Codim<codim>::Geometry == Geometry<dim-codim,dimw>
      /** \brief The type of the geometry associated with the entity.*/
      typedef Dune::Geometry<dim-cd, dimw, const GridImp, GeometryImp> Geometry;
      /** \brief The type of the local geometry associated with the entity.*/
      typedef Dune::Geometry<dim-cd, dim, const GridImp, GeometryImp> LocalGeometry;
      /** \brief The type of the entity. */
      // we could - if needed - introduce another struct for dimglobal of Geometry
      typedef Dune::Entity<cd, dim, const GridImp, EntityImp> Entity;

      /** \brief The type of the iterator over all level entities of this codim. */
      typedef Dune::LevelIterator<cd,All_Partition,const GridImp,LevelIteratorImp> LevelIterator;

      /** \brief The type of the iterator over all leaf entities of this codim. */
      typedef Dune::LeafIterator<cd,All_Partition,const GridImp,LeafIteratorImp> LeafIterator;

      /** \brief The type of the entity pointer for entities of this codim.*/
      typedef Dune::EntityPointer<const GridImp,EntityPointerImp<cd,const GridImp> > EntityPointer;

      /**
       * \brief Traits associated with a specific grid partition type.
       * \tparam pitype The type of the grid partition.
       */
      template <PartitionIteratorType pitype>
      struct Partition
      {
        /** \brief The type of the iterator over the level entities of this codim on this partition. */
        typedef Dune::LevelIterator<cd,pitype,const GridImp,LevelIteratorImp> LevelIterator;
        /** \brief The type of the iterator over the leaf entities of this codim on this partition. */
        typedef Dune::LeafIterator<cd,pitype,const GridImp,LeafIteratorImp> LeafIterator;
      };
    private:
      friend class Dune::Entity<cd, dim, const GridImp, EntityImp>;
      typedef EntityPointerImp<cd,const GridImp> EntityPointerImpl;
    };

    /**
     * \brief Traits associated with a specific grid partition type.
     * \tparam pitype The type of the grid partition.
     */
    template <PartitionIteratorType pitype>
    struct Partition
    {
      /** \brief The type of the level grid view associated with this partition type. */
      typedef Dune::GridView<LevelGridViewTraits<const GridImp,pitype> >
      LevelGridView;

      /** \brief The type of the leaf grid view associated with this partition type. */
      typedef Dune::GridView<LeafGridViewTraits<const GridImp,pitype> >
      LeafGridView;
    };

    /** \brief The type of the level index set. */
    typedef IndexSet<const GridImp,LevelIndexSetImp> LevelIndexSet;
    /** \brief The type of the leaf index set. */
    typedef IndexSet<const GridImp,LeafIndexSetImp> LeafIndexSet;
    /** \brief The type of the global id set. */
    typedef IdSet<const GridImp,GlobalIdSetImp,GIDType> GlobalIdSet;
    /** \brief The type of the local id set. */
    typedef IdSet<const GridImp,LocalIdSetImp,LIDType> LocalIdSet;

    /** \brief The type of the collective communication. */
    typedef CCType CollectiveCommunication;
  };

  // define of capabilties for the interface class
  namespace Capabilities
  {
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
    typedef typename InterfaceType::ImplementationType ImplementationType;
    //! create interface object by calling the contructor of the base class
    explicit MakeableInterfaceObject ( const ImplementationType &realImp )
      : InterfaceType( realImp )
    {}
  };
}

#include "geometry.hh"
#include "entity.hh"
#include "entitypointer.hh"
#include "leveliterator.hh"
#include "intersection.hh"
#include "intersectioniterator.hh"
#include "hierarchiciterator.hh"
#include "leafiterator.hh"
#include "indexidset.hh"

template<int d, int dw, class ct, class gf>
inline std::ostream& operator<< (std::ostream& s,
                                 const Dune::Grid<d,dw,ct,gf> & g)
{
  s << g.name();
  return s;
}

#endif
