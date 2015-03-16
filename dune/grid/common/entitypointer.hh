// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_ENTITYPOINTER_HH
#define DUNE_GRID_ENTITYPOINTER_HH

#include <dune/common/proxymemberaccess.hh>
#include <dune/common/iteratorfacades.hh>
#include <dune/common/deprecated.hh>
#include <dune/grid/common/grid.hh>

/** \file
    \brief Wrapper and interface class for a static iterator (EntityPointer)
 */

#define DUNE_ENTITYPOINTER_DEPRECATED_MSG  DUNE_DEPRECATED_MSG("EntityPointer is deprecated and will be removed after the release of dune-grid-2.4. Instead, you can copy and store entities directly now. Note, this might lead to a decreased performance until all grid implementations properly addressed this interface change.")
namespace Dune
{

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
      And virtualy all Dune::EntityPointer<..., SXxxIterator> are descendents
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
  };

}
#undef DUNE_ENTITYPOINTER_DEPRECATED_MSG

#endif // DUNE_GRID_ENTITYPOINTER_HH
