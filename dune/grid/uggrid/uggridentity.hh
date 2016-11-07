// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGGRIDENTITY_HH
#define DUNE_UGGRIDENTITY_HH

/** \file
 * \brief The UGGridEntity class and its specializations
 */

#include <memory>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/gridenums.hh>


namespace Dune {

  // Forward declarations
  template<int codim, int dim, class GridImp>
  class UGGridEntity;
  template<int dim>
  class UGGrid;
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
  template <class GridType>
  class GridFactory;

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

    typename UG_NS<dim>::template Entity<codim>::T* getTarget() const {
      return this->realEntity.getTarget();
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

    template <class GridImp_>
    friend class UGGridIdSet;

    friend class UGGridEntitySeed<codim, GridImp>;

    typedef typename GridImp::ctype UGCtype;

    typedef typename GridImp::Traits::template Codim<codim>::GeometryImpl GeometryImpl;

    template <class GridType>
    friend class GridFactory;

  public:
    UGGridEntity()
      : target_(nullptr)
      , gridImp_(nullptr)
    {}

    UGGridEntity(typename UG_NS<dim>::template Entity<codim>::T* target, const GridImp* gridImp)
    {
      setToTarget(target,gridImp);
    }

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
      if (codim != dim) {   // other codims are done below
        return InteriorEntity;
      }

      typename UG_NS<dim>::Node *node =
        static_cast<typename UG_NS<dim>::Node *>(target_);

      if (UG_NS<dim>::Priority(node)    == UG_NS<dim>::PrioHGhost
          || UG_NS<dim>::Priority(node) == UG_NS<dim>::PrioVGhost
          || UG_NS<dim>::Priority(node) == UG_NS<dim>::PrioVHGhost)
        return GhostEntity;
      else if (hasBorderCopy_(node))
        return BorderEntity;
      else if (UG_NS<dim>::Priority(node) == UG_NS<dim>::PrioMaster || UG_NS<dim>::Priority(node) == UG_NS<dim>::PrioNone)
        return InteriorEntity;
      else
        DUNE_THROW(GridError, "Unknown priority " << UG_NS<dim>::Priority(node));
#endif
    }

  protected:
#ifdef ModelP
    bool hasBorderCopy_(typename UG_NS<dim>::Node *node) const {
      int  *plist = UG_NS<dim>::DDD_InfoProcList(UG_NS<dim>::ParHdr(node));
      for (int i = 0; plist[i] >= 0; i += 2)
        if (plist[i + 1] == UG_NS<dim>::PrioBorder)
          return true;

      return false;
    }
#endif

  public:


    /*! Intra-element access to entities of codimension cc > codim. Return number of entities
       with codimension cc.
     */
    //!< Default codim 1 Faces and codim == dim Vertices
    template<int cc> int count () const;

    //! geometry of this entity
    Geometry geometry () const { return Geometry( geo_ ); }

    /** \brief Get the seed corresponding to this entity */
    EntitySeed seed () const { return EntitySeed( *this ); }

    /** \brief Return the number of subEntities of codimension codim.
     */
    unsigned int subEntities (unsigned int cd) const
    {
      return ReferenceElements<UGCtype, dim-codim>::general(type()).size(cd-codim);
    }

    typename UG_NS<dim>::template Entity<codim>::T* getTarget() const
    {
      return target_;
    }

    //! equality
    bool equals(const UGGridEntity& other) const
    {
      return getTarget() == other.getTarget();
    }

  private:
    /** \brief Set this entity to a particular UG entity */
    void setToTarget(typename UG_NS<dim>::template Entity<codim>::T* target,const GridImp* gridImp) {
      target_ = target;
      geo_.setToTarget(target);
      gridImp_ = gridImp;
    }

    //! the current geometry
    GeometryImpl geo_;

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

    template <class GridImp_>
    friend class UGGridIdSet;

    friend class UGGridEntitySeed<codim, GridImp>;
    typedef typename GridImp::ctype UGCtype;

    typedef typename GridImp::Traits::template Codim<codim>::GeometryImpl GeometryImpl;

  public:

    typedef typename GridImp::template Codim<codim>::Geometry Geometry;

    UGEdgeEntity()
      : target_(nullptr)
      , gridImp_(nullptr)
    {}

    UGEdgeEntity(typename UG_NS<dim>::template Entity<codim>::T* target, const GridImp* gridImp)
    {
      setToTarget(target,gridImp);
    }

    //! level of the edge
    int level () const {
      return UG_NS<dim>::myLevel(target_);
    }

    /** \brief Return the entity type identifier */
    GeometryType type() const
    {
      return GeometryType(1);
    }

    /** \brief Return the number of subEntities of codimension codim.
     */
    unsigned int subEntities (unsigned int cd) const
    {
      return ReferenceElements<UGCtype, dim-codim>::general(type()).size(cd-codim);
    }

    /** \brief The partition type for parallel computing
     */
    PartitionType partitionType () const
    {
#ifndef ModelP
      return InteriorEntity;
#else

      typename UG_NS<dim>::Edge *edge =
        static_cast<typename UG_NS<dim>::Edge *>(target_);

      if (UG_NS<dim>::Priority(edge)    == UG_NS<dim>::PrioHGhost
          || UG_NS<dim>::Priority(edge) == UG_NS<dim>::PrioVGhost
          || UG_NS<dim>::Priority(edge) == UG_NS<dim>::PrioVHGhost)
        return GhostEntity;
      else if (hasBorderCopy_(edge))
        return BorderEntity;
      else if (UG_NS<dim>::Priority(edge) == UG_NS<dim>::PrioMaster || UG_NS<dim>::Priority(edge) == UG_NS<dim>::PrioNone)
        return InteriorEntity;
      else
        DUNE_THROW(GridError, "Unknown priority " << UG_NS<dim>::Priority(edge));
#endif
    }

    /*! Intra-element access to entities of codimension cc > codim. Return number of entities
       with codimension cc.
     */
    //!< Default codim 1 Faces and codim == dim Vertices
    template<int cc> int count () const;

    //! geometry of this entity
    Geometry geometry () const { return Geometry( *geo_ ); }

  protected:
#ifdef ModelP
    /** \brief Return true if at least one copy of the edge on another processor is marked as 'border' */
    bool hasBorderCopy_(typename UG_NS<dim>::Edge *edge) const {
      int  *plist = UG_NS<dim>::DDD_InfoProcList(UG_NS<dim>::ParHdr(edge));
      for (int i = 0; plist[i] >= 0; i += 2)
        if (plist[i + 1] == UG_NS<dim>::PrioBorder)
          return true;

      return false;
    }
#endif

    /** \brief Set edge object to a UG edge object */
    void setToTarget(typename UG_NS<dim>::template Entity<codim>::T* target, const GridImp* gridImp) {
      target_ = target;

      // obtain the corner coordinates from UG
      UGCtype* cornerCoords[2*dim];
      UG_NS<dim>::Corner_Coordinates(target_, cornerCoords);

      // convert to the type required by MultiLinearGeometry
      std::vector<FieldVector<UGCtype, dim> > geometryCoords(2);
      for(size_t i = 0; i < 2; i++)
        for (size_t j = 0; j < dim; j++)
            geometryCoords[i][j] = cornerCoords[i][j];

      geo_ = std::make_shared<GeometryImpl>(type(), geometryCoords);

      gridImp_ = gridImp;
    }
  public:
    /** \brief Set edge object to a UG edge object using the center element and the side number
     * \note This method is only here to please the compiler, it should never actually be called!
     */
    void setToTarget(typename UG_NS<dim>::Element* target, unsigned int side) {
      DUNE_THROW(Dune::Exception, "Programming error, this method should never be called!");
    }

    typename UG_NS<dim>::template Entity<codim>::T* getTarget() const
    {
      return target_;
    }

    //! equality
    bool equals(const UGEdgeEntity& other) const
    {
      return getTarget() == other.getTarget();
    }

  protected:
    std::shared_ptr<GeometryImpl> geo_;

    typename UG_NS<dim>::template Entity<codim>::T* target_;

    /** \brief gridImp Not actually used, only the codim-0 specialization needs it.
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
  {

  public:
    /** \brief The type of UGGrid Entity seeds */
    typedef typename GridImp::Traits::template Codim<1>::EntitySeed EntitySeed;

    UGGridEntity()
    {}

    UGGridEntity(typename UG_NS<2>::template Entity<1>::T* target, const GridImp* gridImp)
      : UGEdgeEntity<2,GridImp>(target,gridImp)
    {}

    /** \brief Get the seed corresponding to this entity
     *
     * This method cannot be moved to the base class, because the 'this' pointer has the wrong type there.
     */
    EntitySeed seed () const
    {
      return EntitySeed( *this );
    }

  };

  /*! \brief Specialization for edge in 3D
   * \ingroup UGGrid
   */
  template<class GridImp>
  class UGGridEntity<2,3,GridImp>
    : public UGEdgeEntity<3,GridImp>
  {

  public:
    /** \brief The type of UGGrid Entity seeds */
    typedef typename GridImp::Traits::template Codim<2>::EntitySeed EntitySeed;

    UGGridEntity()
    {}

    UGGridEntity(typename UG_NS<3>::template Entity<2>::T* target, const GridImp* gridImp)
      : UGEdgeEntity<3,GridImp>(target,gridImp)
    {}

    /** \brief Get the seed corresponding to this entity
     *
     * This method cannot be moved to the base class, because the 'this' pointer has the wrong type there.
     */
    EntitySeed seed () const
    {
      return EntitySeed( *this );
    }

  };

  /*! \brief UGGrid face entity in 3d grids
   * \ingroup UGGrid
   */
  template<class GridImp>
  class UGGridEntity<1,3,GridImp>
  {
    enum {codim = 1};
    enum {dim = 3};
    typedef typename GridImp::ctype UGCtype;

    typedef typename GridImp::Traits::template Codim<codim>::GeometryImpl GeometryImpl;

  public:

    typedef typename GridImp::template Codim<codim>::Geometry Geometry;

    /** \brief The type of UGGrid Entity seeds */
    typedef typename GridImp::Traits::template Codim<codim>::EntitySeed EntitySeed;

    UGGridEntity()
      : target_(nullptr)
      , gridImp_(nullptr)
    {}

    UGGridEntity(typename UG_NS<dim>::template Entity<codim>::T* target, const GridImp* gridImp)
    {
      setToTarget(target,gridImp);
    }

    /** \brief Return the entity type identifier */
    GeometryType type() const
    {
      // retrieve an element that this face is a side of, and retrieve the side number
      typename UG_NS<dim>::Element* center;
      unsigned int side;
      UG_NS<dim>::GetElementAndSideFromSideVector(target_, center, side);

      switch (UG_NS<dim>::Tag(center)) {

      case UG::D3::TETRAHEDRON :
        return GeometryType(GeometryType::simplex,2);
      case UG::D3::PYRAMID :
        return (side==0)
               ? GeometryType(GeometryType::cube,2)
               : GeometryType(GeometryType::simplex,2);
      case UG::D3::PRISM :
        return (side==0 or side==4)
               ? GeometryType(GeometryType::simplex,2)
               : GeometryType(GeometryType::cube,2);
      case UG::D3::HEXAHEDRON :
        return GeometryType(GeometryType::cube,2);
      default :
        DUNE_THROW(GridError, "UGFaceEntity::type():  ERROR:  Unknown type "
                   << UG_NS<dim>::Tag(center) << " found!");

      }

    }

    //! level of the face
    int level () const {
      return UG_NS<dim>::myLevel(target_);
    }

    /** \brief The partition type for parallel computing */
    PartitionType partitionType () const
    {
#ifndef ModelP
      return InteriorEntity;
#else

      typename UG_NS<dim>::Vector *face = getTarget();

      if (UG_NS<dim>::Priority(face)    == UG_NS<dim>::PrioHGhost
          || UG_NS<dim>::Priority(face) == UG_NS<dim>::PrioVGhost
          || UG_NS<dim>::Priority(face) == UG_NS<dim>::PrioVHGhost)
        return GhostEntity;
      else if (hasBorderCopy_(face))
        return BorderEntity;
      else if (UG_NS<dim>::Priority(face) == UG_NS<dim>::PrioMaster || UG_NS<dim>::Priority(face) == UG_NS<dim>::PrioNone)
        return InteriorEntity;
      else
        DUNE_THROW(GridError, "Unknown priority " << UG_NS<dim>::Priority(face));
#endif
    }

#ifdef ModelP
    bool hasBorderCopy_(typename UG_NS<dim>::Vector *face) const {
      int  *plist = UG_NS<dim>::DDD_InfoProcList(UG_NS<dim>::ParHdr(face));
      for (int i = 0; plist[i] >= 0; i += 2)
        if (plist[i + 1] == UG_NS<dim>::PrioBorder)
          return true;

      return false;
    }
#endif

    //! geometry of this entity
    Geometry geometry () const { return Geometry( *geo_ ); }

    /** \brief Get the seed corresponding to this entity */
    EntitySeed seed () const
    {
      return EntitySeed( *this );
    }

    /** \brief Return the number of subEntities of codimension codim.
     */
    unsigned int subEntities (unsigned int cd) const
    {
      return ReferenceElements<UGCtype, dim-codim>::general(type()).size(cd-codim);
    }

    /** \brief Set to a UG side vector object
        \param target The UG side vector to point to
        \param gridImp Dummy argument, only for consistency with codim-0 entities
     */
    void setToTarget(typename UG_NS<dim>::template Entity<codim>::T* target, const GridImp* gridImp) {
      target_ = target;

      // obtain the corner coordinates from UG
      UGCtype* cornerCoords[4*dim];
      UG_NS<dim>::Corner_Coordinates(target_, cornerCoords);

      // convert to the type required by MultiLinearGeometry
      size_t numCorners = type().isTriangle() ? 3 : 4;
      std::vector<FieldVector<UGCtype, dim> > geometryCoords(numCorners);
      for(size_t i = 0; i < numCorners; i++)
        for (size_t j = 0; j < dim; j++)
            geometryCoords[i][j] = cornerCoords[i][j];

      geo_ = std::make_shared<GeometryImpl>(type(), geometryCoords);

      gridImp_ = gridImp;
    }

    typename UG_NS<dim>::template Entity<codim>::T* getTarget() const
    {
      return target_;
    }

    //! equality
    bool equals(const UGGridEntity& other) const
    {
      return getTarget() == other.getTarget();
    }

protected:
    std::shared_ptr<GeometryImpl> geo_;

    /** \brief The UG object (a side vector) that represents this face */
    typename UG_NS<dim>::template Entity<codim>::T* target_;

    /** \brief gridImp Not actually used, only the codim-0 specialization needs it.
     * But code is simpler if we just keep it everywhere.
     */
    const GridImp* gridImp_;
  };



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

    typedef typename GridImp::ctype UGCtype;

    typedef typename GridImp::Traits::template Codim< 0 >::GeometryImpl GeometryImpl;
    typedef typename GridImp::Traits::template Codim< 0 >::LocalGeometryImpl LocalGeometryImpl;

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
      : target_(nullptr)
      , gridImp_(nullptr)
    {}

    UGGridEntity(typename UG_NS<dim>::Element* target, const GridImp* gridImp)
    {
      setToTarget(target,gridImp);
    }

    //! Level of this element
    int level () const {
      return UG_NS<dim>::myLevel(target_);
    }

    /** \brief Return the entity type identifier */
    GeometryType type() const;

    /** \brief The partition type for parallel computing */
    PartitionType partitionType () const {
#ifndef ModelP
      return InteriorEntity;
#else
      if (UG_NS<dim>::EPriority(target_) == UG_NS<dim>::PrioHGhost
          || UG_NS<dim>::EPriority(target_) == UG_NS<dim>::PrioVGhost
          || UG_NS<dim>::EPriority(target_) == UG_NS<dim>::PrioVHGhost)
        return GhostEntity;
      else
        return InteriorEntity;
#endif
    }

    //! Geometry of this entity
    Geometry geometry () const { return Geometry( geo_ ); }

    /** \brief Get the seed corresponding to this entity */
    EntitySeed seed () const { return EntitySeed( *this ); }

    /** \brief Return the number of subEntities of codimension cc.
     */
    template<int cc>
    int count () const;

    /** \brief Return the number of subEntities of codimension codim.
     */
    unsigned int subEntities (unsigned int codim) const
    {
      if (dim==3) {

        switch (codim) {
        case 0 :
          return 1;
        case 1 :
          return UG_NS<dim>::Sides_Of_Elem(target_);
        case 2 :
          return UG_NS<dim>::Edges_Of_Elem(target_);
        case 3 :
          return UG_NS<dim>::Corners_Of_Elem(target_);
        }

      } else {

        switch (codim) {
        case 0 :
          return 1;
        case 1 :
          return UG_NS<dim>::Edges_Of_Elem(target_);
        case 2 :
          return UG_NS<dim>::Corners_Of_Elem(target_);
        default :
          // do nothing for codim = 3
          break;
        }

      }
      DUNE_THROW(GridError, "You can't call UGGridEntity<0,dim>::count "
             << "with dim==" << dim << " and codim==" << codim << "!");
    }


    /** \brief Provide access to sub entity i of given codimension. Entities
     *  are numbered 0 ... count<cc>()-1
     */
    template<int cc>
    typename GridImp::template Codim<cc>::Entity subEntity (int i) const;

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
    typename GridImp::template Codim<0>::Entity father () const {
      return typename GridImp::template Codim<0>::Entity(UGGridEntity(UG_NS<dim>::EFather(target_),gridImp_));
    }

    //! returns true if father entity exists
    bool hasFather () const
    {
      return UG_NS<dim>::EFather(target_) != NULL;
    }

    /*! Location of this element relative to the reference element element of the father.
     */
    LocalGeometry geometryInFather () const;

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

    typename UG_NS<dim>::template Entity<0>::T* getTarget() const
    {
      return target_;
    }

    //! equality
    bool equals(const UGGridEntity& other) const
    {
      return getTarget() == other.getTarget();
    }

    //! the current geometry
    GeometryImpl geo_;

    /** \brief The corresponding UG-internal data structure */
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
