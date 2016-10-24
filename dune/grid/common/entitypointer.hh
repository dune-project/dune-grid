// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_ENTITYPOINTER_HH
#define DUNE_GRID_ENTITYPOINTER_HH

#include <utility>

#include <dune/common/proxymemberaccess.hh>
#include <dune/common/iteratorfacades.hh>
#include <dune/common/deprecated.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/geometry.hh>
#include <dune/grid/common/gridenums.hh>

/** \file
    \brief Wrapper and interface class for a static iterator (EntityPointer)
 */

#define DUNE_ENTITYPOINTER_DEPRECATED_MSG  DUNE_DEPRECATED_MSG("EntityPointer is deprecated and will be removed after the release of dune-grid-2.4. Instead, you can copy and store entities directly now. Note, this might lead to a decreased performance until all grid implementations properly addressed this interface change.")
namespace Dune
{

  // External forward declaration
  // ----------------------------

  template< int, int, class, template< int, int, class > class >
  class Entity;



  /**
     @brief Wrapper class for pointers to entities

     Template parameters are:

     - <tt>GridImp</tt> Type that is a model of Dune::Grid
     - <tt>IteratorImp</tt> Class template that is a model of Dune::EntityPointer

     <H3>Engine Concept</H3>

     The EntityPointer class template wraps an object of type IteratorImp and forwards all member
     function calls to corresponding members of this class. In that sense EntityPointer
     defines the interface and IteratorImp supplies the implementation.

     <H3>Relation of EntityPointer and Iterators</H3>

      The EntityPointer can be used like a static iterator. It points to a
      Dune::Entity and can be dereferenced, compared and it knows the
      Entity's level.

      You should be able to initialize and interpret every Dune::EntityIterator
      as a Dune::EntityPointer. Therefore we need an inheritance hierarchy of
      the Iterator wrappers:
      \code
      class Dune::EntityPointer< Grid, IteratorImp >;

      class Dune::EntityIterator< codim, Grid, IteratorImp >
      : public Dune::EntityPointer< Grid, IteratorImp >;
      \endcode

      This hierarchy must be mimicked in the implementation (e.g. SGrid):
      \code
      class SEntityPointer<...>;

      class SLevelIterator<...> :
         public SEntityPointer <...>;

      class SHierarchicIterator<...> :
         public SEntityPointer <...>;

      ...
      \endcode
      Please note that dereference(...), equals(...) and level() are only
      implemented in SEntityPointer -- SLevelIterator inherits these methods.
      And it is not possible to specialize these, because EntityPointer always
      uses the base class.

      This leads to a hierarchy where
      Dune::LevelIterator<..., SLevelIterator> inherits
      Dune::EntityPointer<..., SLevelIterator> and
      Dune::HierarchicIterator<..., SHierarchicIterator> inherits
      Dune::EntityPointer<..., SHierarchicIterator>.
      And virtually all Dune::EntityPointer<..., SXxxIterator> are descendents
      of Dune::EntityPointer<..., SEntityPointer>.

      Now you can compare Dune::LevelIterator with Dune::EntityPointer and
      Dune::LeafIterator with Dune::HierarchicIterator. And you can assign
      Dune::EntityPointer from any Dune::XxxIterator class.

      The compiler takes care that you only assign/compare Iterators from the same
      Grid.

      The downside (or advantage) of this inheritance is that you cannot
      use different comparison operators and different dereference
      operators for the different Iterators in one Grid. On the first
      sight it is a downside because one might consider it a good idea
      to have special treatment for different iterators. On the other
      hand it's very confusing for the user if different Iterators show
      different behavior in the same situation. So now they are forced to
      show the same behavior.

      \deprecated The EntityPointer is deprecated and will be removed after the release
                  of dune-grid-2.4. It is not needed anymore because starting with
                  dune-grid-2.4, you can now simply copy and store entities directly.
                  If you need to store many entities for an extended time, use EntitySeed
                  instead. Please note that due to the effort required by this change,
                  those grids that are deprecated in dune-grid-2.4 will not have copyable
                  entities, so if you are forced to use one of those grids, you will have
                  to continue using EntityPointer as well.

      \ingroup GIEntityPointer
   */
  template<class GridImp, class IteratorImp>
  class EntityPointer
  {
    // need to make copy constructor of EntityPointer work for any iterator
    //friend class EntityPointer<GridImp,typename IteratorImp::EntityPointerImp>;
    template< class, class > friend class EntityPointer;

#if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
  public:
#else
  protected:
    // give the GridDefaultImplementation class access to the realImp
    friend class GridDefaultImplementation<
        GridImp::dimension, GridImp::dimensionworld,
        typename GridImp::ctype,
        typename GridImp::GridFamily> ;
#endif
    // type of underlying implementation, for internal use only
    typedef IteratorImp Implementation;

    //! return reference to the real implementation
    Implementation &impl () { return realIterator; }
    //! return reference to the real implementation
    const Implementation &impl () const { return realIterator; }

  protected:
    Implementation realIterator;

  public:
    //! codimension of entity pointer
    enum { codimension = IteratorImp::codimension };

    /** \brief The Entity that this EntityPointer can point to */
    typedef typename IteratorImp::Entity Entity;

    /** \brief Tpy of the reference used when derefencing the Ptr */
    typedef typename std::conditional<
      std::is_lvalue_reference<
        decltype(realIterator.dereference())
        >::value,
      const Entity&,
      Entity
      >::type Reference;

    //===========================================================
    /** @name Constructor & conversion
     */
    //@{
    //===========================================================

    /** \brief Templatized copy constructor from arbitrary IteratorImp.
            This enables that an EntityPointer can be copy-constructed from
            LevelIterator, LeafIterator and HierarchicIterator (because
            these are derived from EntityPointer<...> with their
            corresponding implementation.
     */
    template< class ItImp >
    explicit EntityPointer ( const EntityPointer< GridImp, ItImp > &ep )
      : realIterator( ep.realIterator )
    {}

    /** \brief Default constructor of an empty (undefined) EntityPointer */
    EntityPointer()
    {}

    /** \brief Templatized constructor from type of entity that
        this entity pointer points to. This constructor can be used to
        create an entity pointer from an entity in order to store an
        entity. The implementation of EntityPointer has to have a
        constructor taking a Dune::Entity.
     */
    DUNE_ENTITYPOINTER_DEPRECATED_MSG
    EntityPointer(const Entity& entity)
      : realIterator( entity.impl() )
    {}

    /** \brief Constructor from type of entity implementation that
        this entity pointer points to. This constructor is only
        used in the EntityDefaultImplementation to implement the method
        seed() by default when the type of EntitySeed and EntityPointer coniside.
     */
    EntityPointer ( const typename Entity::Implementation &entityImp )
      : realIterator( entityImp )
    {}

    template< class ItImp >
    DUNE_ENTITYPOINTER_DEPRECATED_MSG
    inline EntityPointer & operator= ( const EntityPointer< GridImp, ItImp > &ep )
    {
      realIterator = ep.realIterator;
      return *this;
    }

    //@}

    //===========================================================
    /** @name Dereferencing
     */
    //@{
    //===========================================================

    // The behavior when dereferencing the EntityPointer facade depends on
    // the way the grid implementation handles returning entities. The implementation
    // may either return a reference to an entity stored inside the EntityPointer
    // implementation or a temporary Entity object. This object has to be forwarded through
    // the facade to the user, which requires a little trickery, especially for operator->().
    //
    // In order to avoid confusing users reading the Doxygen documentation, we provide "clean"
    // function signatures to Doxygen and hide the actual implementations.

#ifdef DOXYGEN

    /** \brief Dereferencing operator. */
    Entity operator*() const
    DUNE_ENTITYPOINTER_DEPRECATED_MSG;


    /** \brief Pointer operator. */
    const Entity* operator->() const
    DUNE_ENTITYPOINTER_DEPRECATED_MSG;

#else // DOXYGEN

    /** \brief Dereferencing operator. */
    Reference
    operator*() const
    DUNE_ENTITYPOINTER_DEPRECATED_MSG
    {
      return realIterator.dereference();
    }

    /** \brief Pointer operator. */
    decltype(handle_proxy_member_access(realIterator.dereference()))
    operator->() const
    DUNE_ENTITYPOINTER_DEPRECATED_MSG
    {
      return handle_proxy_member_access(realIterator.dereference());
    }

    template<typename T>
    // this construction, where the deprecation warning is triggered by a separate function,
    // is slightly convoluted, but I could not get the warning to trigger reliably when attached
    // directly to the cast operator.
    DUNE_DEPRECATED_MSG("The implicit cast from EntityPointer to an Entity reference is DANGEROUS. It's mainly there for writing backwards compatible code that doesn't trigger a deprecation warning for ported grids and must ONLY be used if the returned reference is used in an rvalue-like setting!")
    void trigger_entity_cast_warning() const
    {}

    template<typename T, typename std::enable_if<std::is_same<T,Entity>::value,int>::type = 0>
    operator const T&() const
    {
      static_assert(std::is_same<T,Entity>::value,"invalid cast");
      trigger_entity_cast_warning<T>();
      return realIterator.dereference();
    }

#endif // DOXYGEN

    //@}

    //===========================================================
    /** @name Compare methods
     */
    //@{
    //===========================================================

    /** \brief Checks for equality.
            Only works for EntityPointers and iterators on the same grid.
            Due to the conversion operators one can compare
            all kinds of iterators and EntityPointer.
     */
    template< class ItImp >
    bool operator== ( const EntityPointer< GridImp, ItImp > &rhs ) const
    {
      return equals( rhs );
    }

    /** \brief Checks for inequality.
            Only works for EntityPointers and iterators on the same grid.
            Due to the conversion operators one can compare
            all kinds of iterators and EntityPointer.
     */
    template< class ItImp >
    bool operator!= ( const EntityPointer< GridImp, ItImp > &rhs ) const
    {
      return !equals( rhs );
    }
    //@}

    /** \brief Compares an EntityPointer with an Entity for equality.
     *
     * \deprecated This method only exists for backwards compatibility during the 2.4
     *             release cycle and will be removed after dune-grid-2.4 is released.
     */
    DUNE_ENTITYPOINTER_DEPRECATED_MSG
    bool operator==(const Entity& rhs) const
    {
      return (**this) == rhs;
    }

    /** \brief Compares an EntityPointer with an Entity for inequality.
     *
     * \deprecated This method only exists for backwards compatibility during the 2.4
     *             release cycle and will be removed after dune-grid-2.4 is released.
     */
    DUNE_ENTITYPOINTER_DEPRECATED_MSG
    bool operator!=(const Entity& rhs) const
    {
      return (**this) != rhs;
    }


    //===========================================================
    /** @name Query methods
     */
    //@{
    //===========================================================

    /** \brief Ask for level of entity.
     *
            This method is redundant and is only there for efficiency reasons.
            It allows an implementation to return the level without actually
            constructing the entity.

       \deprecated Will be removed after the release of dune-grid-2.4. Use the
                   method level() from the dereferenced Entity instead.
     */
    int level () const
    DUNE_ENTITYPOINTER_DEPRECATED_MSG
    {
      return realIterator.level();
    }

    //@}


    //===========================================================
    /** @name Implementor interface
     */
    //@{
    //===========================================================

    /** \brief Copy Constructor from an Iterator implementation.

       You can supply LeafIterator, LevelIterator, HierarchicIterator
       or EntityPointer.
     */
    EntityPointer(const IteratorImp & i) :
      realIterator(i) {}

    /** @brief Forward equality check to realIterator */
    template< class ItImp >
    bool equals ( const EntityPointer< GridImp, ItImp > &rhs ) const
    {
      return realIterator.equals( rhs.realIterator );
    }
    //@}

    //===========================================================
    /** @name Methods and Types of the Entity interface. These are here just for transition purpose
     */
    //@{
    //===========================================================

    /** \brief The geometry type of this entity */
    typedef typename GridImp::template Codim<codimension>::Geometry Geometry;

    //! \brief The corresponding entity seed (for storage of entities)
    typedef typename GridImp::template Codim<codimension>::EntitySeed EntitySeed;

    /** \brief The geometry type of this entity when the geometry is expressed
       embedded in the father element.

       This differs from Geometry in particular when dim != dimworld,
       but even when dim == dimworld the implementation may choose to use
       a different type here.
     */
    typedef typename GridImp::template Codim<codimension>::LocalGeometry LocalGeometry;

    /** \brief EntityPointer types of the different codimensions */
    template <int cd>
    struct Codim
    {
      typedef typename GridImp::template Codim<cd>::EntityPointer EntityPointer;
      typedef typename GridImp::template Codim<cd>::Entity Entity;
    };

    /** \brief The codim==0 EntityPointer type */
    // typedef typename Entity::EntityPointer EntityPointer;

    /** \brief The HierarchicIterator type*/
    typedef typename GridImp::HierarchicIterator HierarchicIterator;

    enum {
      //! Know the grid's dimension
      dimension=Entity::dimension
    };
    enum {
      /** \brief Know dimension of the entity */
      mydimension=Entity::dimension
    };

    //! Partition type of this entity
    PartitionType partitionType () const { return realIterator.dereference().partitionType(); }

    /** \brief obtain geometric realization of the entity
     *
     *  Each entity provides an object of type
     *  Dune::Geometry< dimension-codimension, dimensionworld, ... > that
     *  represents the map from a reference element to world coordinates.
     *
     *  \note Previously, the geometry was encapsulated in the entity object and
     *        a const reference was returned.
     *
     *  \note The returned geometry object is guaranteed to remain valid until the
     *        grid is modified (or deleted).
     */
    Geometry geometry () const { return realIterator.dereference().geometry(); }

    /** \brief Return the name of the reference element. The type can
       be used to access the Dune::ReferenceElement.
     */
    GeometryType type () const { return realIterator.dereference().type(); }

    /** \brief Return the entity seed which contains sufficient information
     *  to generate the entity again and uses as little memory as possible
     */
    EntitySeed seed () const { return realIterator.dereference().seed(); }

#define CHECK_CODIM0 int ecodim = codimension, typename std::enable_if<ecodim == 0,int>::type = 0
#define ONLY_CODIM0 template<int ecodim = codimension, typename std::enable_if<ecodim == 0,int>::type = 0>

    template< int codim, CHECK_CODIM0 >
    typename Codim<codim>::Entity
    subEntity ( int i ) const
    {
      return realIterator.dereference().template subEntity< codim >(i);
    }

    /**\brief Return true if entity has a father entity which can be accessed
       using the father() method.
     */
    ONLY_CODIM0
    bool hasFather () const
    {
      return realIterator.dereference().hasFather();
    }

    //! Returns true if the entity is contained in the leaf grid
    ONLY_CODIM0
    bool isLeaf () const
    {
      return realIterator.dereference().isLeaf();
    }

    /** @brief Returns true if element is of regular type in red/green type refinement.
       In bisection or hanging node refinement this is always true.
     */
    ONLY_CODIM0
    bool isRegular() const { return realIterator.dereference().isRegular(); }

    /** \brief Provides information how this element has been subdivided from its
     *         father element.
     *
     *  The returned LocalGeometry is a model of
     *  Dune::Geometry<dimension,dimension,...>, mapping the reference element of
     *  the given entity to the reference element of its father.
     *
     *  This information is sufficient to interpolate all degrees of freedom in
     *  the conforming case.
     *  Nonconforming may require access to neighbors of the father and
     *  calculations with local coordinates.
     *  The on-the-fly case is somewhat inefficient since degrees of freedom may be
     *  visited several times.
     *  If we store interpolation matrices, this is tolerable.
     *  We assume that on-the-fly implementation of interpolation is only done for
     *  simple discretizations.
     *
     *  \note For ghost entities, this method is not guaranteed to be implemented.
     *
     *  \note Previously, the geometry was encapsulated in the entity object and
     *        a const reference was returned.
     *
     *  \note The returned geometry object is guaranteed to remain valid until the
     *        grid is modified (or deleted).
     */
    ONLY_CODIM0
    LocalGeometry geometryInFather () const { return realIterator.dereference().geometryInFather(); }

    /**\brief Inter-level access to elements that resulted from (recursive)
       subdivision of this element.

       \param[in] maxlevel Iterator does not stop at elements with level greater than maxlevel.
       \return Iterator to the first son (level is not greater than maxlevel)

       \note If the partitionType of the Entity is GhostEntity,
           it is not guaranteed that this method is working
           or implemented in general.
           For some grids it might be available, though.
     */
    ONLY_CODIM0
    HierarchicIterator hbegin (int maxLevel) const
    {
      return realIterator.dereference().hbegin(maxLevel);
    }

    /** \brief Returns iterator to one past the last son element

       \note If the partitionType of the Entity is GhostEntity,
             it is not guaranteed that this method is working
             or implemented in general.
             For some grids it might be available, though.
     */
    ONLY_CODIM0
    HierarchicIterator hend (int maxLevel) const
    {
      return realIterator.dereference().hend(maxLevel);
    }

    /**\brief Returns true, if the entity has been created during the last call to adapt()
     */
    ONLY_CODIM0
    bool isNew () const { return realIterator.dereference().isNew(); }

    /**\brief Returns true, if entity might disappear during the next call to adapt().
     * If the method returns false, the entity is guaranteed to still be present after
     * adaptation.
     */
    ONLY_CODIM0
    bool mightVanish () const { return realIterator.dereference().mightVanish(); }

    /**\brief Returns true, if entity has intersections with boundary
     */
    ONLY_CODIM0
    bool hasBoundaryIntersections () const { return realIterator.dereference().hasBoundaryIntersections(); }
    //@}

  };



#ifndef DOXYEN

  // DefaultEntityPointer
  // --------------------

  /* The EntityPointer class defined above has been deprecated. Unitil its
     final removal valid Dune grids still have to provide at least a suitable
     EntityPointer typedef. This class provides a default implementation of an
     entity pointer from a given Dune::Entity type:
     \code
     struct GridFamily
     {
       ...
       typedef ImplementationDefined Entity;
       typedef DefaultEntityPointer<Entity> EntityPointer;
       ...
     };
     \endcode

     This class will retain a possible compatibility support with the
     deprecated interface behavior if the iterator classes in the grid
     implementation provide the following two additional methods:
     \code
     class Iterator
     {
       // dereference() method required by Dune::EntityIterator
       typedef ImplemenatationDefined Entity;
       Entity dereference () const;

       // allow for (deprecated) construction/assignment of EntityPointer from a given iterator
       operator Dune::DefaultEntityPointer<Entity>() const
       {
         return Dune::DefaultEntityPointer<Entity>(dereference());
       }

       // allow for (deprecated) comparison of an iterator with an entity pointer
       bool equals(const Dune::DefaultEntityPointer<Entity> &rhs) const
       {
         return dereference() == rhs.dereference();
       }
     };
     \endcode
  */
  template< class E >
  class DefaultEntityPointer;

  template< int codim, int dim, class Grid, template< int, int, class > class EntityImp >
  class DefaultEntityPointer< Dune::Entity< codim, dim, Grid, EntityImp > >
  {
  public:
    static const int codimension = codim;

    typedef Dune::Entity< codim, dim, Grid, EntityImp > Entity;

    DefaultEntityPointer () {}

    explicit DefaultEntityPointer ( Entity entity )
      : entity_( std::move( entity ) )
    {}

    explicit DefaultEntityPointer ( EntityImp< codim, dim, Grid > entity )
      : entity_( std::move( entity ) )
    {}

    const Entity &dereference () const { return entity_; }

    bool equals ( const DefaultEntityPointer &rhs ) const
    {
      return entity_ == rhs.entity_;
    }

    int level () const { return entity_.level(); }

  private:
    Entity entity_;
  };

#endif // #ifndef DOXYEN

}
#undef DUNE_ENTITYPOINTER_DEPRECATED_MSG

#endif // DUNE_GRID_ENTITYPOINTER_HH
