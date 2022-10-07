// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
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
#include <dune/grid/uggrid/uggridrenumberer.hh>


namespace Dune {

  // Forward declarations
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
    friend class UGGridLeafIterator;

    template <int codim_, PartitionIteratorType PiType_, class GridImp_>
    friend class UGGridLevelIterator;

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

    /** \brief Copy constructor */
    UGGridEntity(const UGGridEntity& other)
      : target_(other.target_)
      , gridImp_(other.gridImp_)
    {
      if constexpr (codim==dim)
        geo_ = other.geo_;
      else
        geo_ = std::make_unique<GeometryImpl>(*other.geo_);
    }

    UGGridEntity(typename UG_NS<dim>::template Entity<codim>::T* target, const GridImp* gridImp)
    {
      setToTarget(target,gridImp);
    }

    typedef typename GridImp::template Codim<codim>::Geometry Geometry;

    /** \brief The type of UGGrid Entity seeds */
    typedef typename GridImp::Traits::template Codim<codim>::EntitySeed EntitySeed;

    //! level of this entity
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
      if (UG_NS<dim>::Priority(target_)    == UG_NS<dim>::PrioHGhost
          || UG_NS<dim>::Priority(target_) == UG_NS<dim>::PrioVGhost
          || UG_NS<dim>::Priority(target_) == UG_NS<dim>::PrioVHGhost)
        return GhostEntity;
      else if (hasBorderCopy())
        return BorderEntity;
      else if (UG_NS<dim>::Priority(target_) == UG_NS<dim>::PrioMaster || UG_NS<dim>::Priority(target_) == UG_NS<dim>::PrioNone)
        return InteriorEntity;
      else
        DUNE_THROW(GridError, "Unknown priority " << UG_NS<dim>::Priority(target_));
#endif
    }

    /** \brief Get the partition type of each copy of a distributed entity
     *
     * This is a non-interface method, intended mainly for debugging and testing.
     *
     * \return Each pair in the return value contains a process number and the corresponding partition type
     */
    std::vector<std::pair<int,PartitionType> > partitionTypes () const
    {
      std::vector<std::pair<int,PartitionType> > result;

#ifndef ModelP
      result.push_back(std::make_pair(0, InteriorEntity));
#elif DUNE_UGGRID_DDD_InfoProcListRange
      const auto range = UG_NS<dim>::DDD_InfoProcListRange(gridImp_->multigrid_->dddContext(), UG_NS<dim>::ParHdr(target_));
      for (auto&& [rank, priority] : range)
      {
        if (priority == UG_NS<dim>::PrioHGhost || priority == UG_NS<dim>::PrioVGhost || priority == UG_NS<dim>::PrioVHGhost)
          result.push_back(std::make_pair(rank, GhostEntity));
        else
        {
          // The entity is not ghost.  If it is (UG)PrioBorder somewhere, it is (Dune)Border.
          // Otherwise it is (Dune)Interior.
          bool hasBorderCopy = false;
          for (auto&& [rank2, priority2] : range) {
            if (priority2 == UG_NS<dim>::PrioBorder)
            {
              hasBorderCopy = true;
              break;
            }
          }

          result.push_back(std::make_pair(rank, (hasBorderCopy) ? BorderEntity : InteriorEntity));
        }
      }
#else
      int *plist = UG_NS<dim>::DDD_InfoProcList(gridImp_->multigrid_->dddContext(),
                                                UG_NS<dim>::ParHdr(target_));

      for (int i = 0; plist[i] >= 0; i += 2)
      {
        int rank = plist[i];
        auto priority = plist[i + 1];

        if (priority == UG_NS<dim>::PrioHGhost || priority == UG_NS<dim>::PrioVGhost || priority == UG_NS<dim>::PrioVHGhost)
          result.push_back(std::make_pair(rank, GhostEntity));
        else
        {
          // The entity is not ghost.  If it is (UG)PrioBorder somewhere, it is (Dune)Border.
          // Otherwise it is (Dune)Interior.
          bool hasBorderCopy = false;
          for (int i = 0; plist[i] >= 0; i += 2)
            if (plist[i + 1] == UG_NS<dim>::PrioBorder)
            {
              hasBorderCopy = true;
              break;
            }

          result.push_back(std::make_pair(rank, (hasBorderCopy) ? BorderEntity : InteriorEntity));
        }
      }
#endif
      return result;
    }

  protected:
#ifdef ModelP
    // \todo Unify with the following method
    bool hasBorderCopy() const
    {
#if DUNE_UGGRID_DDD_InfoProcListRange
      for (auto&& [rank, priority] : UG_NS<dim>::DDD_InfoProcListRange(
             gridImp_->multigrid_->dddContext(),
             UG_NS<dim>::ParHdr(target_))) {
        if (priority == UG_NS<dim>::PrioBorder)
          return true;
      }
#else
      int  *plist = UG_NS<dim>::DDD_InfoProcList(
        gridImp_->multigrid_->dddContext(),
        UG_NS<dim>::ParHdr(target_));
      for (int i = 0; plist[i] >= 0; i += 2)
        if (plist[i + 1] == UG_NS<dim>::PrioBorder)
          return true;
#endif

      return false;
    }
#endif


  public:

    //! geometry of this entity
    Geometry geometry () const
    {
      if constexpr ((dim-codim)==0)
        return Geometry( geo_ );
      else
        return Geometry( *geo_ );
    }

    /** \brief Get the seed corresponding to this entity */
    EntitySeed seed () const { return EntitySeed( *this ); }

    /** \brief Return the number of subEntities of codimension codim.
     */
    unsigned int subEntities (unsigned int cd) const
    {
      return referenceElement<UGCtype, dim-codim>(type()).size(cd-codim);
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
      gridImp_ = gridImp;
      target_ = target;

      if constexpr ((dim-codim)==1)   // Edge entity
      {
        // Obtain the corner coordinates from UG
        UGCtype* cornerCoords[2*dim];
        UG_NS<dim>::Corner_Coordinates(target_, cornerCoords);

        // convert to the type required by MultiLinearGeometry
        std::vector<FieldVector<UGCtype, dim> > geometryCoords(2);
        for (size_t i=0; i < 2; i++)
          for (size_t j=0; j < dim; j++)
            geometryCoords[i][j] = cornerCoords[i][j];

        geo_ = std::make_unique<GeometryImpl>(type(), geometryCoords);
      }
      else if constexpr ((dim-codim)==2)   // Facet entity
      {
        // obtain the corner coordinates from UG
        UGCtype* cornerCoords[4*dim];
        UG_NS<dim>::Corner_Coordinates(target_, cornerCoords);

        // convert to the type required by MultiLinearGeometry
        size_t numCorners = type().isTriangle() ? 3 : 4;
        std::vector<FieldVector<UGCtype, dim> > geometryCoords(numCorners);
        for(size_t i = 0; i < numCorners; i++)
          for (size_t j = 0; j < dim; j++)
            geometryCoords[UGGridRenumberer<dim-1>::verticesUGtoDUNE(i, type())][j] = cornerCoords[i][j];

        geo_ = std::make_unique<GeometryImpl>(type(), geometryCoords);
      }
      else
        geo_.setToTarget(target);
    }

    // The geometric realization of the entity.  It is a native UG type for vertices,
    // and a dune-geometry MultiLinearGeometry for edges and facets.
    // TODO: Maybe use MultiLinearGeometry for all these cases!
    std::conditional_t<(dim-codim)==0, GeometryImpl, std::unique_ptr<GeometryImpl> > geo_;

    /** \brief The UG object that represents this entity
     * - 0d entities: node
     * - 1d entities: edge
     * - 2d entities in 3d grids: side vector
     */
    typename UG_NS<dim>::template Entity<codim>::T* target_;

    const GridImp* gridImp_;
  };

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
    friend class UGGridLeafIterator;

    template <int codim_, PartitionIteratorType PiType_, class GridImp_>
    friend class UGGridLevelIterator;

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

    /** \brief Get the partition type of each copy of a distributed entity
     *
     * This is a non-interface method, intended mainly for debugging and testing.
     *
     * \return Each pair in the return value contains a process number and the corresponding partition type
     */
    std::vector<std::pair<int,PartitionType> > partitionTypes () const
    {
      std::vector<std::pair<int,PartitionType> > result;

#ifndef ModelP
      result.push_back(std::make_pair(0, InteriorEntity));
#elif DUNE_UGGRID_DDD_InfoProcListRange
      for (auto&& [rank, priority] : UG_NS<dim>::DDD_InfoProcListRange(
             gridImp_->multigrid_->dddContext(), &target_->ge.ddd)) {
        if (priority == UG_NS<dim>::PrioHGhost || priority == UG_NS<dim>::PrioVGhost || priority == UG_NS<dim>::PrioVHGhost)
          result.push_back(std::make_pair(rank, GhostEntity));
        else
          result.push_back(std::make_pair(rank, InteriorEntity));
      }
#else
      int *plist = UG_NS<dim>::DDD_InfoProcList(gridImp_->multigrid_->dddContext(),
                                                &target_->ge.ddd);

      for (int i = 0; plist[i] >= 0; i += 2)
      {
        int rank = plist[i];
        auto priority = plist[i + 1];

        if (priority == UG_NS<dim>::PrioHGhost || priority == UG_NS<dim>::PrioVGhost || priority == UG_NS<dim>::PrioVHGhost)
          result.push_back(std::make_pair(rank, GhostEntity));
        else
          result.push_back(std::make_pair(rank, InteriorEntity));
      }
#endif
      return result;
    }

    //! Geometry of this entity
    Geometry geometry () const { return Geometry( geo_ ); }

    /** \brief Get the seed corresponding to this entity */
    EntitySeed seed () const { return EntitySeed( *this ); }

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
      DUNE_THROW(GridError, "You can't call UGGridEntity<0,dim>::subEntities "
             << "with dim==" << dim << " and codim==" << codim << "!");
    }


    /** \brief Provide access to sub entity i of given codimension. Entities
     *  are numbered 0 ... subEntities(cc)-1
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
      return UG_NS<dim>::EFather(target_) != nullptr;
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
