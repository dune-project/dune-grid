// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_ENTITY_HH
#define DUNE_GRID_ENTITY_HH

#include <dune/common/typetraits.hh>
#include "grid.hh"
#include "entitypointer.hh"

namespace Dune
{

  /**
     @brief Wrapper class for entities


     Template parameters are:

     - <tt>cd</tt> Codimension of the entity
     - <tt>dim</tt> Dimension of the grid
     - <tt>GridImp</tt> Type that is a model of Dune::Grid
     - <tt>EntityImp</tt> Class template that is a model of Dune::Entity


     <H3>Engine Concept</H3>

     This class wraps a object of type EntityImp and forwards all member
     function calls to corresponding members of this class. In that sense Entity
     defines the interface and EntityImp supplies the implementation.
     For various reasons we do not use an inheritance hierarchy and the
     Barton-Nackman trick here.


     <H3>Specialization</H3>

     The Entity class template is specialized for <tt>cd=0</tt> (elements,
     Dune::Entity<0,dim,GridImp,EntityImp>).
     This case has an extended interface.
     The methods defined in the general template
     are provided by the specialization as well. We did not use inheritance
     because different implementations for different codimensions may be required
     and virtual functions had to be avoided.

     <H3>View concept</H3>

     Entities can not be created, assigned or otherwise modified outside
     the interface in the user code. They are only accessible by immutable
     iterators provided on the corresponding grid class.

     The only way to modify the entities of a grid is through grid adaptation which
     consists of tagging entities (of codimension 0) for refinement and then
     calling the adapt() method on the grid.


     \ingroup GIEntity
     \nosubgrouping
   */
  template<int cd, int dim, class GridImp, template<int,int,class> class EntityImp>
  class Entity
  {
#if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
  public:
#else
  protected:
    // give the GridDefaultImplementation class access to the realImp
    friend class GridDefaultImplementation<
        GridImp::dimension, GridImp::dimensionworld,
        typename GridImp::ctype,
        typename GridImp::GridFamily> ;

    // Default*GridView classes need access to intersection iterators
    template< class, PartitionIteratorType > friend class DefaultLevelGridView;
    template< class, PartitionIteratorType > friend class DefaultLeafGridView;
#endif
    // type of underlying implementation, for internal use only
    typedef EntityImp< cd, dim, GridImp > Implementation;

    //! Return reference to the real implementation
    Implementation &impl () { return realEntity; }
    //! Return const reference to the real implementation
    const Implementation &impl () const { return realEntity; }

  protected:
    Implementation realEntity;

  public:

    //===========================================================
    /** @name Exported types and constants
     */
    //@{
    //===========================================================

    //! \brief The corresponding geometry type
    typedef typename GridImp::template Codim<cd>::Geometry Geometry;

    //! \brief The corresponding entity seed (for storage of entities)
    typedef typename GridImp::template Codim<cd>::EntitySeed EntitySeed;

    enum {
      //! \brief Know your own codimension.
      codimension=cd
    };
    enum {
      //! \brief Know the grid dimension.
      dimension=dim
    };
    enum {
      //! \brief Dimensionality of the reference element of the entity.
      mydimension=dim-cd
    };
    enum {
      //! \brief Know the dimension of world.
      dimensionworld=GridImp::dimensionworld
    };

    //! @brief coordinate type of the Grid
    typedef typename GridImp::ctype ctype;
    //@}



    //===========================================================
    /** @name Methods shared by entities of all codimensions
     */
    //@{
    //===========================================================

    //! The level of this entity
    int level () const { return realEntity.level(); }

    //! Partition type of this entity
    PartitionType partitionType () const { return realEntity.partitionType(); }

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
    Geometry geometry () const { return realEntity.geometry(); }
    //@}

    /** \brief Return the name of the reference element. The type can
       be used to access the Dune::ReferenceElement.
     */
    GeometryType type () const { return realEntity.type(); }

    /** \brief Return the entity seed which contains sufficient information
     *  to generate the entity again and uses as less memory as possible
     */
    EntitySeed seed () const { return realEntity.seed(); }

    //===========================================================
    /** @name Interface for the implementor
     */
    //@{
    //===========================================================

    //! Copy constructor from EntityImp
    explicit Entity(const EntityImp<cd,dim,GridImp> & e) : realEntity(e) {}

    //@}

  protected:
    typedef typename remove_const<GridImp>::type mutableGridImp;

    //===========================================================
    /** @name Protected methods
     */
    //@{
    //===========================================================

    // need to make copy constructor of EntityPointer work for any iterator
    //friend class Dune::EntityPointer<GridImp,
    //                                 typename GridImp::GridFamily::Traits::template Codim<cd>::EntityPointerImpl>;
    template< class, class > friend class Dune::EntityPointer;

  protected:
    /** hide copy constructor */
    Entity(const Entity& rhs) : realEntity(rhs.realEntity) {}
    /** hide assignment operator */
    Entity & operator = (const Entity& rhs) {
      realEntity = rhs.realEntity;
      return *this;
    }
    //@}
  };

  /**
     @brief Template specialization of Dune::Entity for Elements (codim==0)

     @see Dune::Entity (general version) for the full documentation

     \extends Entity<int cd, int dim, class GridImp, template<int,int,class> class EntityImp>

     \ingroup GIEntity
     \nosubgrouping
   */
  template<int dim, class GridImp, template<int,int,class> class EntityImp>
  class Entity <0,dim,GridImp,EntityImp>
  {
#if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
  public:
#else
  protected:
    // give the GridDefaultImplementation class access to the realImp
    friend class GridDefaultImplementation<
        GridImp::dimension, GridImp::dimensionworld,
        typename GridImp::ctype,
        typename GridImp::GridFamily> ;

    // Default*GridView classes need access to intersection iterators
    template< class, PartitionIteratorType > friend class DefaultLevelGridView;
    template< class, PartitionIteratorType > friend class DefaultLeafGridView;
#endif
    // type of underlying implementation, for internal use only
    typedef EntityImp< 0, dim, GridImp > Implementation;

    //! Return reference to the real implementation
    Implementation &impl () { return realEntity; }
    //! Return const reference to the real implementation
    const Implementation &impl () const { return realEntity; }

  protected:
    Implementation realEntity;

    typedef typename remove_const<GridImp>::type mutableGridImp;

  public:

    //===========================================================
    /** @name Exported types and constants
     */
    //@{
    //===========================================================

    /** \brief The geometry type of this entity */
    typedef typename GridImp::template Codim<0>::Geometry Geometry;

    //! \brief The corresponding entity seed (for storage of entities)
    typedef typename GridImp::template Codim<0>::EntitySeed EntitySeed;

    /** \brief The geometry type of this entity when the geometry is expressed
       embedded in the father element.

       This differs from Geometry in particular when dim != dimworld,
       but even when dim == dimworld the implementation may choose to use
       a different type here.
     */
    typedef typename GridImp::template Codim<0>::LocalGeometry LocalGeometry;

    /** \brief EntityPointer types of the different codimensions */
    template <int cd>
    struct Codim
    {
      typedef typename GridImp::template Codim<cd>::EntityPointer EntityPointer;
    };

    /** \brief The codim==0 EntityPointer type */
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;

    /** \brief The Dune::IntersectionIterator type for the LeafGridView */
    typedef typename GridImp::LeafIntersectionIterator LeafIntersectionIterator;

    /** \brief The Dune::IntersectionIterator type for the LevelGridView */
    typedef typename GridImp::LevelIntersectionIterator LevelIntersectionIterator;

    /** \brief The HierarchicIterator type*/
    typedef typename GridImp::HierarchicIterator HierarchicIterator;

    enum {
      //! Know your own codimension
      codimension=0
    };
    enum {
      //! Know the grid's dimension
      dimension=dim
    };
    enum {
      /** \brief Know dimension of the entity */
      mydimension=dim
    };
    enum {
      //! Know the world dimension
      dimensionworld=GridImp::dimensionworld
    };
    //! Type used for coordinates
    typedef typename GridImp::ctype ctype;
    //@}


    //===========================================================
    /** @name Methods shared by entities of all codimensions
     */
    //@{
    //===========================================================

    //! @copydoc Dune::Entity::level()
    int level () const { return realEntity.level(); }

    //! @copydoc Dune::Entity::partitionType()
    PartitionType partitionType () const { return realEntity.partitionType(); }

    //! @copydoc Dune::Entity::geometry()
    Geometry geometry () const { return realEntity.geometry(); }
    //@}

    /** \brief Return the name of the reference element. The type can
        be used to access the Dune::ReferenceElement.
     */
    GeometryType type () const { return realEntity.type(); }

    /** \brief Return the entity seed which contains sufficient information
     *  to generate the entity again and uses as little memory as possible
     */
    EntitySeed seed () const { return realEntity.seed(); }

    //===========================================================
    /** @name Extended interface of entities of codimension 0
     */
    //@{
    //===========================================================

    /**\brief Number of subentities with codimension <tt>cc</tt>. This method is in
       principle redundant because this information can be obtained via the reference
       element of the geometry. It is there for efficiency reasons and to make
       the interface self-contained.
     */
    template<int cc> int count () const { return realEntity.template count<cc>(); }

    /** \brief obtain a pointer to a subentity
     *
     *  \tparam  codim  codimension of the desired subentity
     *
     *  \param[in]  i  number of the subentity (in generic numbering)
     *
     *  \returns an EntityPointer to the specified subentity
     *
     *  \note The subentities are numbered 0, ..., count< codim >-1
     */
    template< int codim >
    typename Codim< codim >::EntityPointer subEntity ( int i ) const
    {
      return realEntity.template subEntity< codim >( i );
    }

    /**\brief Access to intersections with neighboring leaf elements.
       A neighbor is an entity of codimension 0
       which has an intersection of codimension 1 in common with this entity.
       Access to those neighbors is provided using the IntersectionIterator.
       This method returns an iterator refering to the first neighbor.

       \note If the partitionType of the Entity is GhostEntity,
             this method might give you only one neighbor, which is the
             interior Entity the GhostEntity is connected to.

       \deprecated This method is deprecated and will be removed after
                   Dune 2.3. Use LeafGridView.ibegin(Entity) instead.
     */
    LeafIntersectionIterator ileafbegin () const
    DUNE_DEPRECATED_MSG("Use LeafGridView.ibegin(Entity) instead.")
    {
      return realEntity.ileafbegin();
    }

    /**\brief Reference to an IntersectionIterator one past the last intersection

       \note If the partitionType of the Entity is GhostEntity,
             this method might give you only one neighbor, which is the
             interior Entity the GhostEntity is connected to.

       \deprecated This method is deprecated and will be removed after
                   Dune 2.3. Use LeafGridView.iend(Entity) instead.
     */
    LeafIntersectionIterator ileafend () const
    DUNE_DEPRECATED_MSG("Use LeafGridView.iend(Entity) instead.")
    {
      return realEntity.ileafend();
    }

    /**\brief Intra-level access to intersections with neighboring elements.
       A neighbor is an entity of codimension 0
       which has an intersection of codimension 1 in common with this entity.
       Access to those neighbors is provided using the IntersectionIterator.
       This method returns an iterator refering to the first neighbor.

       \note If the partitionType of the Entity is GhostEntity,
             this method might give you only one neighbor, which is the
             interior Entity the GhostEntity is connected to.

       \deprecated This method is deprecated and will be removed after
                   Dune 2.3. Use LevelGridView.ibegin(Entity) instead.
     */
    LevelIntersectionIterator ilevelbegin () const
    DUNE_DEPRECATED_MSG("Use LevelGridView.ibegin(Entity) instead.")
    {
      return realEntity.ilevelbegin();
    }

    /**\brief Reference to an IntersectionIterator one past the last intersection

       \note If the partitionType of the Entity is GhostEntity,
             this method might give you only one neighbor, which is the
             interior Entity the GhostEntity is connected to.

       \deprecated This method is deprecated and will be removed after
                   Dune 2.3. Use LevelGridView.iend(Entity) instead.
     */
    LevelIntersectionIterator ilevelend () const
    DUNE_DEPRECATED_MSG("Use LevelGridView.iend(Entity) instead.")
    {
      return realEntity.ilevelend();
    }

    /**\brief Inter-level access to father entity on the next-coarser grid.
       The given entity resulted directly from a subdivision of its father
       entity. For the macro elements dereferencing the EntityPointer is undefined.

       \note If the partitionType of the Entity is GhostEntity,
             it is not guaranteed that this method is working
             or implemented in general.
             For some grids it might be available, though.
     */
    EntityPointer father () const
    {
      return realEntity.father();
    }

    /**\brief Return true if entity has a father entity which can be accessed
       using the father() method.
     */
    bool hasFather () const
    {
      return realEntity.hasFather();
    }

    //! Returns true if the entity is contained in the leaf grid
    bool isLeaf () const
    {
      return realEntity.isLeaf();
    }

    /** @brief Returns true if element is of regular type in red/green type refinement.
       In bisection or hanging node refinement this is always true.
     */
    bool isRegular() const { return realEntity.isRegular(); }

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
    LocalGeometry geometryInFather () const { return realEntity.geometryInFather(); }

    /**\brief Inter-level access to elements that resulted from (recursive)
       subdivision of this element.

       \param[in] maxlevel Iterator does not stop at elements with level greater than maxlevel.
       \return Iterator to the first son (level is not greater than maxlevel)

       \note If the partitionType of the Entity is GhostEntity,
           it is not guaranteed that this method is working
           or implemented in general.
           For some grids it might be available, though.
     */
    HierarchicIterator hbegin (int maxlevel) const
    {
      return realEntity.hbegin(maxlevel);
    }

    /** \brief Returns iterator to one past the last son element

       \note If the partitionType of the Entity is GhostEntity,
             it is not guaranteed that this method is working
             or implemented in general.
             For some grids it might be available, though.
     */
    HierarchicIterator hend (int maxlevel) const
    {
      return realEntity.hend(maxlevel);
    }

    /**\brief Returns true, if the entity has been created during the last call to adapt()
     */
    bool isNew () const { return realEntity.isNew(); }

    /**\brief Returns true, if entity might disappear during the next call to adapt().
     * If the method returns false, the entity is guaranteed to still be present after
     * adaptation.
     */
    bool mightVanish () const { return realEntity.mightVanish(); }

    //===========================================================
    /** @name Interface for the implementor
     */
    //@{
    //===========================================================
    /**\brief Returns true, if entity has intersections with boundary, see
         default implementation
     */
    bool hasBoundaryIntersections () const { return realEntity.hasBoundaryIntersections(); }

    //! Copy constructor from EntityImp
    explicit Entity(const EntityImp<0,dim,GridImp> & e) : realEntity(e) {}

    //@}


  protected:
    //===========================================================
    /** @name Protected methods
     */
    //@{
    //===========================================================

    // need to make copy constructor of EntityPointer work for any iterator
    //friend class Dune::EntityPointer<GridImp,
    //                                 typename GridImp::GridFamily::Traits::template Codim<0>::EntityPointerImpl>;
    template< class, class > friend class Dune::EntityPointer;

  protected:
    /** hide copy constructor */
    Entity(const Entity& rhs) : realEntity(rhs.realEntity) {}
    /** hide assignement operator */
    Entity & operator = (const Entity& rhs) {
      realEntity = rhs.realEntity;
      return *this;
    }
    //@}
  };


  //********************************************************************
  /**
     @brief Default Implementations for EntityImp

     EntityDefaultImplementation provides default implementations for Entity which uses
     the implemented interface which has to be done by the user.

     @note this is the general version, but there is a specialization for cd=0

     @ingroup GridDevel
   */
  template<int cd, int dim, class GridImp, template<int,int,class> class EntityImp>
  class EntityDefaultImplementation
  {
  public:
    //! know your own codimension
    enum { codimension=cd };

    //! know your own dimension
    enum { dimension=dim };

    /** \brief Know dimension of the entity */
    enum { mydimension=dim-cd };

    //! know your own dimension of world
    enum { dimensionworld=GridImp::dimensionworld };

    //! define type used for coordinates in grid module
    typedef typename GridImp::ctype ctype;

    //! \brief The corresponding entity seed (for storage of entities)
    typedef typename GridImp::template Codim<cd>::EntitySeed EntitySeed;

    //! \brief The corresponding entity seed (for storage of entities)
    typedef typename GridImp::template Codim<cd>::EntityPointer EntityPointer;

    /** \brief Return the name of the reference element. The type can
        be used to access the Dune::ReferenceElement.
     */
    GeometryType type () const { return asImp().geometry().type(); }

  private:
    //!  Barton-Nackman trick
    EntityImp<cd,dim,GridImp>& asImp ()
    {
      return static_cast<EntityImp<cd,dim,GridImp>&>(*this);
    }
    const EntityImp<cd,dim,GridImp>& asImp () const
    {
      return static_cast<const EntityImp<cd,dim,GridImp>&>(*this);
    }
  }; // end EntityDefaultImplementation

  //********************************************************************
  /**
     @brief Default Implementations for EntityImp (Elements [cd=0])

     EntityDefaultImplementation provides default implementations for Entity which uses
     the implemented interface which has to be done by the user.

     \extends EntityDefaultImplementation<int cd, int dim, class GridImp, template<int,int,class> class EntityImp>

     @ingroup GridDevel
   */
  template<int dim, class GridImp, template<int,int,class> class EntityImp>
  class EntityDefaultImplementation <0,dim,GridImp,EntityImp>
  {
  public:
    //! know your own codimension
    enum { codimension=0 };

    //! know your own dimension
    enum { dimension=dim };

    /** \brief Know dimension of the entity */
    enum { mydimension=dim };

    //! know your own dimension of world
    enum { dimensionworld=GridImp::dimensionworld };

    //! define type used for coordinates in grid module
    typedef typename GridImp::ctype ctype;

    //! \brief The corresponding entity seed (for storage of entities)
    typedef typename GridImp::template Codim<0>::EntitySeed EntitySeed;

    //! \brief The corresponding entity seed (for storage of entities)
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;

    /** @brief Returns true if element is of regular type in red/green type refinement.
       In bisection or hanging node refinement this is always true.
     */
    bool isRegular() const { return true; }

    /** \brief Return the name of the reference element. The type can
        be used to access the Dune::ReferenceElement.
     */
    GeometryType type () const { return asImp().geometry().type(); }

    /**\brief Returns true, if the entity has been created during the last call to adapt()
     */
    bool isNew () const { return false; }

    /**\brief Returns true, if entity might disappear during the next call to adapt()
     */
    bool mightVanish () const { return false; }

    /**\brief Returns true, if entity has intersections with boundary,
       this implementation uses the Level- and LeafIntersectionIterator to
       check for boundary intersections
     */
    bool hasBoundaryIntersections () const
    {
      {
        typedef typename GridImp::LevelIntersectionIterator IntersectionIterator;
        IntersectionIterator end = asImp().ilevelend();
        for(IntersectionIterator it = asImp().ilevelbegin(); it != end; ++it)
        {
          if( it->boundary() ) return true;
        }
      }

      {
        typedef typename GridImp::LeafIntersectionIterator IntersectionIterator;
        IntersectionIterator end = asImp().ileafend();
        for(IntersectionIterator it = asImp().ileafbegin(); it != end; ++it)
        {
          if( it->boundary() ) return true;
        }
      }

      return false;
    }

  private:
    //  Barton-Nackman trick
    EntityImp<0,dim,GridImp>& asImp () { return static_cast<EntityImp<0,dim,GridImp>&>(*this); }
    const EntityImp<0,dim,GridImp>& asImp () const { return static_cast<const EntityImp<0,dim,GridImp>&>(*this); }
  };

}

#endif // DUNE_GRID_ENTITY_HH
