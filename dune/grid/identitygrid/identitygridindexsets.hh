// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IDENTITYGRID_INDEXSETS_HH
#define DUNE_IDENTITYGRID_INDEXSETS_HH

/** \file
 * \brief The index and id sets for the IdentityGrid class
 */

#include <dune/grid/common/indexidset.hh>

#include <vector>

namespace Dune {

  /** \todo Take the index types from the host grid */
  template<class GridImp>
  class IdentityGridLevelIndexSet :
    public IndexSet<GridImp,
                    IdentityGridLevelIndexSet<GridImp>,
                    typename std::remove_const<GridImp>::type::HostGridType::LevelGridView::IndexSet::IndexType,
                    typename std::remove_const<GridImp>::type::HostGridType::LevelGridView::IndexSet::Types
                    >
  {
  public:

    typedef typename std::remove_const<GridImp>::type::HostGridType HostGrid;
    typedef typename HostGrid::LevelGridView::IndexSet::Types Types;

    constexpr static int dim = GridImp::dimension;

    //! get index of an entity
    template<int codim>
    int index (const typename GridImp::Traits::template Codim<codim>::Entity& e) const
    {
      return grid_->hostgrid_->levelIndexSet(level_).template index<codim>(grid_->template getHostEntity<codim>(e));
    }


    //! get index of subEntity of a codim 0 entity
    template<int cc>
    int subIndex (const typename GridImp::Traits::template Codim<cc>::Entity& e, int i, int codim) const
    {
      return grid_->hostgrid_->levelIndexSet(level_).subIndex(grid_->template getHostEntity<cc>(e), i, codim);
    }


    //! get number of entities of given codim, type and on this level
    std::size_t size (int codim) const {
      return grid_->hostgrid_->levelIndexSet(level_).size(codim);
    }


    //! get number of entities of given codim, type and on this level
    std::size_t size (GeometryType type) const
    {
      return grid_->hostgrid_->levelIndexSet(level_).size(type);
    }


    /** \brief Deliver all geometry types used in this grid */
    const std::vector<GeometryType>& geomTypes (int codim) const
    {
      return grid_->hostgrid_->levelIndexSet(level_).geomTypes(codim);
    }

    /** \brief Deliver all geometry types used in this grid */
    Types types (int codim) const
    {
      return grid_->hostgrid_->levelIndexSet(level_).types(codim);
    }

    /** \brief Return true if the given entity is contained in the index set */
    template<class EntityType>
    bool contains (const EntityType& e) const
    {
      return grid_->hostgrid_->levelIndexSet(level_).contains(grid_->template getHostEntity<EntityType::codimension>(e));
    }

    /** \brief Set up the index set */
    void update(const GridImp& grid, int level)
    {
      grid_ = &grid;
      level_ = level;
    }


    GridImp* grid_;

    int level_;
  };


  template<class GridImp>
  class IdentityGridLeafIndexSet :
    public IndexSet<GridImp,
                    IdentityGridLeafIndexSet<GridImp>,
                    typename std::remove_const<GridImp>::type::HostGridType::LeafGridView::IndexSet::IndexType,
                    typename std::remove_const<GridImp>::type::HostGridType::LeafGridView::IndexSet::Types
                    >
  {
    typedef typename std::remove_const<GridImp>::type::HostGridType HostGrid;

  public:

    typedef typename HostGrid::LevelGridView::IndexSet::Types Types;

    /*
     * We use the remove_const to extract the Type from the mutable class,
     * because the const class is not instantiated yet.
     */
    constexpr static int dim = std::remove_const<GridImp>::type::dimension;


    //! constructor stores reference to a grid and level
    IdentityGridLeafIndexSet (const GridImp& grid)
      : grid_(&grid)
    {}


    //! get index of an entity
    /*
        We use the RemoveConst to extract the Type from the mutable class,
        because the const class is not instantiated yet.
     */
    template<int codim>
    int index (const typename std::remove_const<GridImp>::type::template Codim<codim>::Entity& e) const
    {
      return grid_->hostgrid_->leafIndexSet().template index<codim>(grid_->template getHostEntity<codim>(e));
    }


    //! get index of subEntity of a codim 0 entity
    /*
        We use the RemoveConst to extract the Type from the mutable class,
        because the const class is not instantiated yet.
     */
    template<int cc>
    int subIndex (const typename std::remove_const<GridImp>::type::Traits::template Codim<cc>::Entity& e, int i, int codim) const
    {
      return grid_->hostgrid_->leafIndexSet().subIndex(grid_->template getHostEntity<cc>(e),i, codim);
    }


    //! get number of entities of given type
    std::size_t size (GeometryType type) const
    {
      return grid_->hostgrid_->leafIndexSet().size(type);
    }


    //! get number of entities of given codim
    std::size_t size (int codim) const
    {
      return grid_->hostgrid_->leafIndexSet().size(codim);
    }


    /** \brief Deliver all geometry types used in this grid */
    const std::vector<GeometryType>& geomTypes (int codim) const
    {
      return grid_->hostgrid_->leafIndexSet().geomTypes(codim);
    }

    /** \brief Deliver all geometry types used in this grid */
    Types types (int codim) const
    {
      return grid_->hostgrid_->leafIndexSet().types(codim);
    }

    /** \brief Return true if the given entity is contained in the index set */
    template<class EntityType>
    bool contains (const EntityType& e) const
    {
      return grid_->hostgrid_->leafIndexSet().contains(grid_->template getHostEntity<EntityType::codimension>(e));
    }



    /** \todo Currently we support only vertex and element indices */
    void update(const GridImp& grid)
    {
      grid_ = &grid;
    }


    GridImp* grid_;
  };




  template <class GridImp>
  class IdentityGridGlobalIdSet :
    public IdSet<GridImp,IdentityGridGlobalIdSet<GridImp>,
        typename std::remove_const<GridImp>::type::HostGridType::Traits::GlobalIdSet::IdType>
  {

    typedef typename std::remove_const<GridImp>::type::HostGridType HostGrid;


  public:
    //! constructor stores reference to a grid
    IdentityGridGlobalIdSet (const GridImp& g) : grid_(&g) {}

    //! define the type used for persistent indices
    typedef typename HostGrid::Traits::GlobalIdSet::IdType IdType;


    //! get id of an entity
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    template<int cd>
    IdType id (const typename std::remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      // Return id of the host entity
      return grid_->hostgrid_->globalIdSet().id(e.impl().hostEntity_);
    }


    //! get id of subEntity
    /*
        We use the remove_const to extract the Type from the mutable class,
        because the const class is not instantiated yet.
     */
    IdType subId (const typename std::remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i, int codim) const
    {
      // Return sub id of the host entity
      return grid_->hostgrid_->globalIdSet().subId(e.impl().hostEntity_,i, codim);
    }


    /** \todo Should be private */
    void update() {}


    const GridImp* grid_;
  };




  template<class GridImp>
  class IdentityGridLocalIdSet :
    public IdSet<GridImp,IdentityGridLocalIdSet<GridImp>,
        typename std::remove_const<GridImp>::type::HostGridType::Traits::LocalIdSet::IdType>
  {
  private:

    typedef typename std::remove_const<GridImp>::type::HostGridType HostGrid;


  public:
    //! define the type used for persistent local ids
    typedef typename HostGrid::Traits::LocalIdSet::IdType IdType;


    //! constructor stores reference to a grid
    IdentityGridLocalIdSet (const GridImp& g) : grid_(&g) {}


    //! get id of an entity
    /*
        We use the remove_const to extract the Type from the mutable class,
        because the const class is not instantiated yet.
     */
    template<int cd>
    IdType id (const typename std::remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      // Return id of the host entity
      return grid_->hostgrid_->localIdSet().id(e.impl().hostEntity_);
    }


    //! get id of subEntity
    /*
     * We use the remove_const to extract the Type from the mutable class,
     * because the const class is not instantiated yet.
     */
    IdType subId (const typename std::remove_const<GridImp>::type::template Codim<0>::Entity& e, int i, int codim) const
    {
      // Return sub id of the host entity
      return grid_->hostgrid_->localIdSet().subId(e.impl().hostEntity_,i,codim);
    }


    /** \todo Should be private */
    void update() {}


    const GridImp* grid_;
  };


}  // namespace Dune


#endif
