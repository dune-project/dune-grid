// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGGRIDENTITY_HH
#define DUNE_UGGRIDENTITY_HH

/** \file
 * \brief The UGGridEntity class and its specializations
 */

#include <dune/grid/common/gridenums.hh>

#include "uggridrenumberer.hh"


namespace Dune {

  // Forward declarations
  template<int codim, int dim, class GridImp>
  class UGGridEntity;
  template<int dim>
  class UGGrid;
  template<int codim, class GridImp>
  class UGGridEntityPointer;
  template<int codim, class GridImp>
  class UGGridEntitySeed;
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class UGGridLevelIterator;
  template<class GridImp>
  class UGGridLevelIntersectionIterator;
  template<class GridImp>
  class UGGridLeafIntersectionIterator;
  template<class GridImp>
  class UGGridHierarchicIterator;

  template<int codim, int dim, class GridImp>
  class UGMakeableEntity :
    public GridImp::template Codim<codim>::Entity
  {
  public:

    UGMakeableEntity(typename UG_NS<dim>::template Entity<codim>::T* target, const GridImp* gridImp) :
      GridImp::template Codim<codim>::Entity (UGGridEntity<codim, dim, const GridImp>())
    {
      this->realEntity.setToTarget(target,gridImp);
    }

    UGMakeableEntity() :
      GridImp::template Codim<codim>::Entity (UGGridEntity<codim, dim, const GridImp>())
    {}

    void setToTarget(typename UG_NS<dim>::template Entity<codim>::T* target, const GridImp* gridImp) {
      this->realEntity.setToTarget(target,gridImp);
    }

    typename UG_NS<dim>::template Entity<codim>::T* getTarget() {
      return this->realEntity.target_;
    }

  };

  //**********************************************************************
  //
  // --UGGridEntity
  // --Entity
  //
  /** \brief The implementation of entities in a UGGrid
     \ingroup UGGrid

     A Grid is a container of grid entities. An entity is parametrized by the codimension.
     An entity of codimension c in dimension d is a d-c dimensional object.

   */
  template<int codim, int dim, class GridImp>
  class UGGridEntity
  {

    template <int codim_, PartitionIteratorType PiType_, class GridImp_>
    friend class UGGridLevelIterator;

    friend class UGMakeableEntity<codim,dim,const UGGrid<dim> >;
    friend class UGMakeableEntity<codim,dim,UGGrid<dim> >;

    friend class UGGrid<dim>;

    template <class GridImp_>
    friend class UGGridLevelIndexSet;

    template <class GridImp_>
    friend class UGGridLeafIndexSet;

    template <class GridImp_, bool Local>
    friend class UGGridIdSet;

    friend class UGGridEntityPointer<codim, GridImp>;
    friend class UGGridEntitySeed<codim, GridImp>;

    typedef typename GridImp::ctype UGCtype;

  public:

    UGGridEntity()
      : geo_(UGGridGeometry<dim-codim,dim,GridImp>())
    {}

    typedef typename GridImp::template Codim<codim>::Geometry Geometry;

    /** \brief The type of UGGrid Entity seeds */
    typedef typename GridImp::Traits::template Codim<codim>::EntitySeed EntitySeed;

    //! level of this element
    int level () const {
      return UG_NS<dim>::myLevel(target_);
    }

    /** \brief Return the entity type identifier */
    GeometryType type() const;

    /** \brief The partition type for parallel computing */
    PartitionType partitionType () const
    {
#ifndef ModelP
      return InteriorEntity;
#else
#define PARHDRE(p) (&((p)->ddd))
#define EPRIO(e) (PARHDRE(e)->prio)
      if (codim != dim) {
        // TODO: faces (elements and edges are done below)
        return InteriorEntity;
      }

      typename UG_NS<dim>::Node *node =
        static_cast<typename UG_NS<dim>::Node *>(target_);

      if (EPRIO(node)    == UG_NS<dim>::PrioHGhost
          || EPRIO(node) == UG_NS<dim>::PrioVGhost
          || EPRIO(node) == UG_NS<dim>::PrioVHGhost)
        return GhostEntity;
      else if (EPRIO(node) == UG_NS<dim>::PrioBorder || hasBorderCopy_(node))
        return BorderEntity;
      else if (EPRIO(node) == UG_NS<dim>::PrioMaster || EPRIO(node) == UG_NS<dim>::PrioNone)
        return InteriorEntity;
      else
        DUNE_THROW(GridError, "Unknown priority " << EPRIO(node));

#undef EPRIO
#undef PARHDRE
#endif
    }

  protected:
#ifdef ModelP
    bool hasBorderCopy_(typename UG_NS<dim>::Node *node) const {
#define PARHDR(p)         (&((p)->ddd))
#define PRIOLIST(n)        UG_NS<dim>::DDD_InfoProcList(PARHDR(n))

      int  *plist = PRIOLIST(node);
      for (int i = 0; plist[i] >= 0; i += 2)
        if (plist[i + 1] == UG_NS<dim>::PrioBorder)
          return true;

      return false;
#undef PARHDR
#undef PRIOLIST
    }
#endif

  public:


    /*! Intra-element access to entities of codimension cc > codim. Return number of entities
       with codimension cc.
     */
    //!< Default codim 1 Faces and codim == dim Vertices
    template<int cc> int count () const;

    //! geometry of this entity
    const Geometry& geometry () const {return geo_;}

    /** \brief Get the seed corresponding to this entity */
    EntitySeed seed () const { return EntitySeed( *this ); }

  private:
    /** \brief Set this entity to a particular UG entity */
    void setToTarget(typename UG_NS<dim>::template Entity<codim>::T* target,const GridImp* gridImp) {
      target_ = target;
      GridImp::getRealImplementation(geo_).setToTarget(target);
      gridImp_ = gridImp;
    }

    //! the current geometry
    MakeableInterfaceObject<Geometry> geo_;

    typename UG_NS<dim>::template Entity<codim>::T* target_;

    /** \brief gridImp Not actually used, only the codim-0 specialization needs it
     * But code is simpler if we just keep it everywhere.
     */
    const GridImp* gridImp_;
  };

  /*! \brief Edge entity
   * \ingroup UGGrid
   */
  template<int dim, class GridImp>
  class UGEdgeEntity
  {
    enum {codim = dim - 1};
    template <int codim_, PartitionIteratorType PiType_, class GridImp_>
    friend class UGGridLevelIterator;

    friend class UGMakeableEntity<codim,dim,const UGGrid<dim> >;
    friend class UGMakeableEntity<codim,dim,UGGrid<dim> >;

    friend class UGGrid<dim>;

    template <class GridImp_>
    friend class UGGridLevelIndexSet;

    template <class GridImp_>
    friend class UGGridLeafIndexSet;

    template <class GridImp_, bool Local>
    friend class UGGridIdSet;

    friend class UGGridEntityPointer<codim, GridImp>;
    friend class UGGridEntitySeed<codim, GridImp>;
    typedef typename GridImp::ctype UGCtype;

  public:

    typedef typename GridImp::template Codim<codim>::Geometry Geometry;

    /** \brief The type of UGGrid Entity seeds */
    typedef typename GridImp::Traits::template Codim<codim>::EntitySeed EntitySeed;

    //! level of this element
    int level () const {
      return UG_NS<dim>::myLevel(target_);
    }

    /** \brief Return the entity type identifier */
    GeometryType type() const
    {
      return GeometryType(1);
    }

    /** \brief The partition type for parallel computing
     * \todo So far it always returns InteriorEntity */
    PartitionType partitionType () const
    {
#ifndef ModelP
      return InteriorEntity;
#else

#define PARHDRE(p) (&((p)->ddd))
#define EPRIO(e) (PARHDRE(e)->prio)
      typename UG_NS<dim>::Edge *edge =
        static_cast<typename UG_NS<dim>::Edge *>(target_);

      if (EPRIO(edge)    == UG_NS<dim>::PrioHGhost
          || EPRIO(edge) == UG_NS<dim>::PrioVGhost
          || EPRIO(edge) == UG_NS<dim>::PrioVHGhost)
        return GhostEntity;
      else if (EPRIO(edge) == UG_NS<dim>::PrioBorder || hasBorderCopy_(edge))
        return BorderEntity;
      else if (EPRIO(edge) == UG_NS<dim>::PrioMaster || EPRIO(edge) == UG_NS<dim>::PrioNone)
        return InteriorEntity;
      else
        DUNE_THROW(GridError, "Unknown priority " << EPRIO(edge));
#undef EPRIO
#undef PARHDRE
#endif
    }

    /*! Intra-element access to entities of codimension cc > codim. Return number of entities
       with codimension cc.
     */
    //!< Default codim 1 Faces and codim == dim Vertices
    template<int cc> int count () const;

    /** \brief Get the seed corresponding to this entity */
    EntitySeed seed () const { return EntitySeed( *this ); }

  protected:
#ifdef ModelP
    bool hasBorderCopy_(typename UG_NS<dim>::Edge *edge) const {
#define PARHDR(p)         (&((p)->ddd))
#define PRIOLIST(n)        UG_NS<dim>::DDD_InfoProcList(PARHDR(n))

      int  *plist = PRIOLIST(edge);
      for (int i = 0; plist[i] >= 0; i += 2)
        if (plist[i + 1] == UG_NS<dim>::PrioBorder)
          return true;

      return false;
#undef PARHDR
#undef PRIOLIST
    }
#endif

    void setToTarget(typename UG_NS<dim>::template Entity<codim>::T* target, const GridImp* gridImp) {
      target_ = target;
      gridImp_ = gridImp;
    }

    typename UG_NS<dim>::template Entity<codim>::T* target_;

    /** \brief gridImp Not actually used, only the codim-0 specialization needs it
     * But code is simpler if we just keep it everywhere.
     */
    const GridImp* gridImp_;
  };

  /*! \brief Specialization for edge in 2D
   * \ingroup UGGrid
   */
  template<class GridImp>
  class UGGridEntity<1,2,GridImp>
    : public UGEdgeEntity<2,GridImp>
  {};

  /*! \brief Specialization for edge in 3D
   * \ingroup UGGrid
   */
  template<class GridImp>
  class UGGridEntity<2,3,GridImp>
    : public UGEdgeEntity<3,GridImp>
  {};

  //***********************
  //
  //  --UGGridEntity
  //  --0Entity
  //
  //***********************

  /** \brief Specialization for codim-0-entities.
   * \ingroup UGGrid
   *
   * This class embodies the topological parts of elements of the grid.
   * It has an extended interface compared to the general entity class.
   * For example, Entities of codimension 0  allow to visit all neighbors.
   *
   * UGGrid only implements the cases dim==dimworld==2 and dim=dimworld==3.
   */
  template<int dim, class GridImp>
  class UGGridEntity<0,dim,GridImp>
  {
    friend class UGGrid<dim>;
    friend class UGGridLeafIntersectionIterator <GridImp>;
    friend class UGGridHierarchicIterator <GridImp>;

    template <int codim_, PartitionIteratorType PiType_, class GridImp_>
    friend class UGGridLevelIterator;

    friend class UGMakeableEntity<0,dim,GridImp>;

    friend class UGGridEntityPointer<0, GridImp>;

    typedef typename GridImp::ctype UGCtype;

  public:

    typedef typename GridImp::template Codim<0>::Geometry Geometry;

    typedef typename GridImp::template Codim<0>::LocalGeometry LocalGeometry;

    //! The Iterator over neighbors on this level
    typedef UGGridLevelIntersectionIterator<GridImp> LevelIntersectionIterator;

    //! The Iterator over neighbors on the leaf level
    typedef UGGridLeafIntersectionIterator<GridImp> LeafIntersectionIterator;

    //! Iterator over descendants of the entity
    typedef UGGridHierarchicIterator<GridImp> HierarchicIterator;

    /** \brief The type of UGGrid Entity seeds */
    typedef typename GridImp::Traits::template Codim<0>::EntitySeed EntitySeed;

    UGGridEntity()
      : geo_(UGGridGeometry<dim,dim,GridImp>()),
        geometryInFather_(UGGridLocalGeometry<dim,dim,GridImp>()),
        gridImp_(NULL)
    {}

    //! Level of this element
    int level () const {
      return UG_NS<dim>::myLevel(target_);
    }

    /** \brief Return the entity type identifier */
    GeometryType type() const;

    /** \brief The partition type for parallel computing
        \todo Do not copy macro from UG */
    PartitionType partitionType () const {
#ifndef ModelP
      return InteriorEntity;
#else
#define PARHDRE(p) (&((p)->ge.ddd))
#define EPRIO(e) DDD_InfoPriority(PARHDRE(e))
      if (EPRIO(target_) == UG_NS<dim>::PrioHGhost
          || EPRIO(target_) == UG_NS<dim>::PrioVGhost
          || EPRIO(target_) == UG_NS<dim>::PrioVHGhost)
        return GhostEntity;
      else
        return InteriorEntity;
#undef EPRIO
#undef PARHDRE
#endif
    }

    //! Geometry of this entity
    const Geometry& geometry () const {
      return geo_;
    }

    /** \brief Get the seed corresponding to this entity */
    EntitySeed seed () const { return EntitySeed( *this ); }

    /** \brief Return the number of subEntities of codimension cc.
     */
    template<int cc>
    int count () const;

    /** \brief Provide access to sub entity i of given codimension. Entities
     *  are numbered 0 ... count<cc>()-1
     */
    template<int cc>
    typename GridImp::template Codim<cc>::EntityPointer subEntity (int i) const;

    /** \todo It would be faster to not use -1 as the end marker but
        number of sides instead */
    UGGridLeafIntersectionIterator<GridImp> ileafbegin () const {
      return UGGridLeafIntersectionIterator<GridImp>(target_, (isLeaf()) ? 0 : UG_NS<dim>::Sides_Of_Elem(target_),gridImp_);
    }

    UGGridLevelIntersectionIterator<GridImp> ilevelbegin () const {
      return UGGridLevelIntersectionIterator<GridImp>(target_, 0, gridImp_);
    }

    //! Reference to one past the last leaf neighbor
    UGGridLeafIntersectionIterator<GridImp> ileafend () const {
      return UGGridLeafIntersectionIterator<GridImp>(target_, UG_NS<dim>::Sides_Of_Elem(target_), gridImp_);
    }

    //! Reference to one past the last level neighbor
    UGGridLevelIntersectionIterator<GridImp> ilevelend () const {
      return UGGridLevelIntersectionIterator<GridImp>(target_, UG_NS<dim>::Sides_Of_Elem(target_),gridImp_);
    }

    //! returns true if Entity has NO children
    bool isLeaf() const {
      return UG_NS<dim>::isLeaf(target_);
    }

    //! returns true if element is of regular type
    bool isRegular() const {
      return UG_NS<dim>::isRegular(target_);
    }

    /** \brief Returns true if the entity has intersections with the boundary
     */
    bool hasBoundaryIntersections() const {
      return UG_NS<dim>::isBoundaryElement(target_);
    }

    //! Inter-level access to father element on coarser grid.
    //! Assumes that meshes are nested.
    typename GridImp::template Codim<0>::EntityPointer father () const {
      return typename GridImp::template Codim<0>::EntityPointer (UGGridEntityPointer<0,GridImp>(UG_NS<dim>::EFather(target_),gridImp_));
    }

    //! returns true if father entity exists
    bool hasFather () const
    {
      return UG_NS<dim>::EFather(target_) != NULL;
    }

    /*! Location of this element relative to the reference element element of the father.
     */
    const LocalGeometry& geometryInFather () const;

    /*! Inter-level access to son elements on higher levels<=maxlevel.
       This is provided for sparsely stored nested unstructured meshes.
       Returns iterator to first son.
     */
    UGGridHierarchicIterator<GridImp> hbegin (int maxlevel) const;

    //! Returns iterator to one past the last son
    UGGridHierarchicIterator<GridImp> hend (int maxlevel) const;

    //***************************************************************
    //  Interface for Adaptation
    //***************************************************************

    //! returns true, if entity was refined during last adaptation cycle
    bool isNew() const;

    /*! \brief
       returns true, if entity might be coarsened during next adaptation
       cycle, which is true for entities that have been marked for
       coarsening or for entities that are not regular (i.e. isRegular
       returns false) */
    bool mightVanish() const;

    //!
    void setToTarget(typename UG_NS<dim>::Element* target, const GridImp* gridImp);

    //! the current geometry
    MakeableInterfaceObject<Geometry> geo_;

    //! geometry for mapping into father's reference element
    mutable MakeableInterfaceObject<LocalGeometry> geometryInFather_;

    typename UG_NS<dim>::Element* target_;

    /** \brief Pointer to the grid that we are part of.
     *
     * We need that only to hand it over to the intersections,
     * which need it.
     */
    const GridImp* gridImp_;

  }; // end of UGGridEntity codim = 0

} // namespace Dune

#endif
