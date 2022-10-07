// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_ENTITY_HH
#define DUNE_GRID_COMMON_ENTITY_HH

#include <type_traits>

#include <dune/common/iteratorrange.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/dimension.hh>
#include <dune/geometry/referenceelements.hh>

#include "grid.hh"
#include "rangegenerators.hh"

namespace Dune
{

  /**
     @brief Wrapper class for entities

     \tparam cd Codimension of the entity
     \tparam dim Dimension of the grid
     \tparam GridImp Type that is a model of Dune::Grid
     \tparam EntityImp Class template that is a model of Dune::Entity


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
  public:
    /**
     * \brief type of underlying implementation
     *
     * \warning Implementation details may change without prior notification.
     **/
    typedef EntityImp< cd, dim, GridImp > Implementation;

    /**
     * \brief access to the underlying implementation
     *
     * \warning Implementation details may change without prior notification.
     **/
    Implementation &impl () { return realEntity; }
    /**
     * \brief access to the underlying implementation
     *
     * \warning Implementation details may change without prior notification.
     **/
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

    //! \brief Know your own codimension.
    constexpr static int codimension = cd;

    //! \brief Know the grid dimension.
    constexpr static int dimension = dim;

    //! \brief Dimensionality of the reference element of the entity.
    constexpr static int mydimension = dim - cd;
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

    /** \brief Return the name of the reference element. The type can
       be used to access the Dune::ReferenceElement.
     */
    GeometryType type () const { return realEntity.type(); }

    /**
     * \brief Number of subentities for a given codimension
     *
     * \param  codim  codimension to obtain number of subentities for
     *
     * \note The codimension is specified with respect to the grid dimension.
     *
     * \note Unless the geometry type is None, this method is redundant and
     *       the same information can be obtained from the corresponding
     *       reference element.
     **/
    unsigned int subEntities ( unsigned int codim ) const
    {
      return realEntity.subEntities(codim);
    }

    /** \brief Return the entity seed which contains sufficient information
     *  to generate the entity again and uses as little memory as possible
     */
    EntitySeed seed () const { return realEntity.seed(); }

    //! Compares two entities for equality.
    bool operator==(const Entity& other) const
    {
      return realEntity.equals(other.realEntity);
    }

    //! Compares two entities for inequality.
    bool operator!=(const Entity& other) const
    {
      return !realEntity.equals(other.realEntity);
    }

    Entity()
    {}

    //! Copy constructor from an existing entity.
    Entity(const Entity& other)
      : realEntity(other.realEntity)
    {}

    //! Move constructor from an existing entity.
    Entity(Entity&& other)
      : realEntity(std::move(other.realEntity))
    {}

    //! Copy assignment operator from an existing entity.
    Entity& operator=(const Entity& other)
    {
      realEntity = other.realEntity;
      return *this;
    }

    //! Move assignment operator from an existing entity.
    Entity& operator=(Entity&& other)
    {
      realEntity = std::move(other.realEntity);
      return *this;
    }

    //@}

    //===========================================================
    /** @name Interface for the implementor
     */
    //@{
    //===========================================================

    //! Copy constructor from EntityImp
    Entity(const EntityImp<cd,dim,GridImp> & e) : realEntity(e) {}

    //! Move constructor from EntityImp
    Entity(EntityImp<cd,dim,GridImp> && e) : realEntity(std::move(e)) {}

    //@}
  };

  /**
     @brief Template specialization of Dune::Entity for Elements (codim==0)

     \tparam dim Dimension of the grid
     \tparam GridImp Type that is a model of Dune::Grid
     \tparam EntityImp Class template that is a model of Dune::Entity

     @see Dune::Entity (general version) for the full documentation

     \extends Entity<int cd, int dim, class GridImp, template<int,int,class> class EntityImp>

     \ingroup GIEntity
     \nosubgrouping
   */
  template<int dim, class GridImp, template<int,int,class> class EntityImp>
  class Entity <0,dim,GridImp,EntityImp>
  {
  public:
    /**
     * \brief Type of underlying implementation
     *
     * \note This code may change without prior warning
     *
     **/
    typedef EntityImp< 0, dim, GridImp > Implementation;

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

    /** \brief Entity types of the different codimensions */
    template <int cd>
    struct Codim
    {
      typedef typename GridImp::template Codim<cd>::Entity Entity;
    };

    /** \brief The HierarchicIterator type*/
    typedef typename GridImp::HierarchicIterator HierarchicIterator;

    //! Know your own codimension
    constexpr static int codimension = 0;

    //! Know the grid's dimension
    constexpr static int dimension = dim;

    /** \brief Know dimension of the entity */
    constexpr static int mydimension = dim;
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

    /**
     * \brief Number of subentities for a given codimension
     *
     * \param  codim  codimension to obtain number of subentities for
     *
     * \note The codimension is specified with respect to the grid dimension.
     *
     * \note Unless the geometry type is None, this method is redundant and
     *       the same information can be obtained from the corresponding
     *       reference element.
     **/
    unsigned int subEntities ( unsigned int codim ) const
    {
      return realEntity.subEntities(codim);
    }

    /** \brief Return the name of the reference element. The type can
        be used to access the Dune::ReferenceElement.
     */
    GeometryType type () const { return realEntity.type(); }

    /** \brief Return the entity seed which contains sufficient information
     *  to generate the entity again and uses as little memory as possible
     */
    EntitySeed seed () const { return realEntity.seed(); }

    //! Compares two entities for equality.
    bool operator==(const Entity& other) const
    {
      return realEntity.equals(other.realEntity);
    }

    //! Compares two entities for inequality.
    bool operator!=(const Entity& other) const
    {
      return !realEntity.equals(other.realEntity);
    }

    Entity()
    {}

    //! Copy constructor from an existing entity.
    Entity(const Entity& other)
      : realEntity(other.realEntity)
    {}

    //! Move constructor from an existing entity.
    Entity(Entity&& other)
      : realEntity(std::move(other.realEntity))
    {}

    //! Copy assignment operator from an existing entity.
    Entity& operator=(const Entity& other)
    {
      realEntity = other.realEntity;
      return *this;
    }

    //! Move assignment operator from an existing entity.
    Entity& operator=(Entity&& other)
    {
      realEntity = std::move(other.realEntity);
      return *this;
    }

    //@}

    //===========================================================
    /** @name Extended interface of entities of codimension 0
     */
    //@{
    //===========================================================

    /** \brief Obtain a subentity
     *
     *  \tparam  codim  codimension of the desired subentity
     *
     *  \param[in]  i  number of the subentity (in generic numbering)
     *
     *  \returns the specified subentity
     *
     *  \note The subentities are numbered 0, ..., subEntities( codim )-1
     */
    template< int codim >
    typename Codim< codim >::Entity
    subEntity ( int i ) const
    {
      return realEntity.template subEntity< codim >( i );
    }

    /**\brief Inter-level access to father entity on the next-coarser grid.
       The given entity resulted directly from a subdivision of its father
       entity. The behaviour for elements on the macro grid, that is when
       \ref hasFather() is false, is undefined.

       \note If the partitionType of the Entity is GhostEntity,
             it is not guaranteed that this method is working
             or implemented in general.
             For some grids it might be available, though.
     */
    Entity father () const
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

       \param[in] maxLevel Iterator does not stop at elements with level greater than maxlevel.
       \return Iterator to the first son (level is not greater than maxlevel)

       \note If the partitionType of the Entity is GhostEntity,
           it is not guaranteed that this method is working
           or implemented in general.
           For some grids it might be available, though.
     */
    HierarchicIterator hbegin (int maxLevel) const
    {
      return realEntity.hbegin(maxLevel);
    }

    /** \brief Returns iterator to one past the last son element

       \note If the partitionType of the Entity is GhostEntity,
             it is not guaranteed that this method is working
             or implemented in general.
             For some grids it might be available, though.
     */
    HierarchicIterator hend (int maxLevel) const
    {
      return realEntity.hend(maxLevel);
    }

    /**\brief Returns true, if the entity has been created during the last call to adapt()
     */
    bool isNew () const { return realEntity.isNew(); }

    /**\brief Returns true, if entity might disappear during the next call to adapt().
     * If the method returns false, the entity is guaranteed to still be present after
     * adaptation.
     */
    bool mightVanish () const { return realEntity.mightVanish(); }

    /**\brief Returns true, if entity has intersections with boundary
     */
    bool hasBoundaryIntersections () const { return realEntity.hasBoundaryIntersections(); }


    //===========================================================
    /** @name Interface for the implementor
     */
    //@{
    //===========================================================

    //! Copy constructor from EntityImp
    Entity(const EntityImp<0,dim,GridImp> & e) : realEntity(e) {}

    //! Move constructor from EntityImp
    Entity(EntityImp<0,dim,GridImp> && e) : realEntity(std::move(e)) {}

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
    constexpr static int codimension = cd;

    //! Dimension of the grid
    constexpr static int dimension = dim;

    /** \brief Know dimension of the entity */
    constexpr static int mydimension = dim - cd;

    //! \brief The corresponding entity seed (for storage of entities)
    typedef typename GridImp::template Codim<cd>::EntitySeed EntitySeed;

    /**
     * \brief Number of subentities for a given codimension
     *
     * \param  codim  codimension to obtain number of subentities for
     *
     * \note The codimension is specified with respect to the grid dimension.
     *
     * \note Unless the geometry type is None, this method is redundant and
     *       the same information can be obtained from the corresponding
     *       reference element.
     **/
    unsigned int subEntities ( unsigned int codim ) const
    {
      typedef typename std::remove_const< GridImp >::type::ctype ctype;
      return ReferenceElements< ctype, mydimension >::general( asImp().type() ).size( codim - codimension );
    }

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
    constexpr static int codimension = 0;

    //! Dimension of the grid
    constexpr static int dimension = dim;

    /** \brief Know dimension of the entity */
    constexpr static int mydimension = dim;

    //! \brief The corresponding entity seed (for storage of entities)
    typedef typename GridImp::template Codim<0>::EntitySeed EntitySeed;

    /** @brief Returns true if element is of regular type in red/green type refinement.
       In bisection or hanging node refinement this is always true.
     */
    bool isRegular() const { return true; }

    /**
     * \brief Number of subentities for a given codimension
     *
     * \param  codim  codimension to obtain number of subentities for
     *
     * \note The codimension is specified with respect to the grid dimension.
     *
     * \note Unless the geometry type is None, this method is redundant and
     *       the same information can be obtained from the corresponding
     *       reference element.
     **/
    unsigned int subEntities ( unsigned int codim ) const
    {
      typedef typename std::remove_const< GridImp >::type::ctype ctype;
      return ReferenceElements< ctype, mydimension >::general( asImp().type() ).size( codim - codimension );
    }

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
      typedef typename GridImp::LevelIntersectionIterator IntersectionIterator;
      IntersectionIterator end = asImp().ilevelend();
      for (IntersectionIterator it = asImp().ilevelbegin(); it != end; ++it)
        if( it->boundary() )
          return true;

      return false;
    }

  private:
    //  Barton-Nackman trick
    EntityImp<0,dim,GridImp>& asImp () { return static_cast<EntityImp<0,dim,GridImp>&>(*this); }
    const EntityImp<0,dim,GridImp>& asImp () const { return static_cast<const EntityImp<0,dim,GridImp>&>(*this); }
  };

  //! Second-level dispatch to select the correct reference element for a grid entity.
  /**
   * This function is the default implementation of the second-level reference element dispatch
   * performed by Entity.
   *
   * When referenceElement() is called with an Entity, it will forward the call to
   * `referenceElement<ctype, mydim>(const GeometryType&)`. This default implementation
   * will do the right thing as long as your geometry is based on a standard Dune ReferenceElement. If
   * it is not and you want to supply your own reference element implementation, provide an override of
   * this function for your specific geometry implementation.
   *
   * \related Entity
   */
  template< int cd, int dim, class GridImp, template<int,int,class> class EntityImp >
  auto referenceElement(const Entity< cd, dim, GridImp, EntityImp >& entity )
    -> decltype(referenceElement<typename GridImp::ctype,GridImp::template Codim<cd>::Geometry::mydimension>(entity.type()))
  {
    typedef typename GridImp::template Codim<cd>::Geometry Geo;
    return referenceElement< typename Geo::ctype, Geo::mydimension >(entity.type());
  }
}

#endif // DUNE_GRID_COMMON_ENTITY_HH
