// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGGRIDLEVELITERATOR_HH
#define DUNE_UGGRIDLEVELITERATOR_HH

/** \file
 * \brief The UGGridLevelIterator class
 */

#include <dune/grid/uggrid/uggridentity.hh>

namespace Dune {

  //**********************************************************************
  //
  // --UGGridLevelIterator
  // --LevelIterator
  /** \brief Iterator over all entities of a given codimension and level of a grid.
   * \ingroup UGGrid
   */
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class UGGridLevelIterator
  {
    constexpr static int dim = GridImp::dimension;

    friend class UGGridEntity<codim,GridImp::dimension,GridImp>;
    friend class UGGridEntity<0,    GridImp::dimension,GridImp>;

    // The type of the UG entity we're pointing to
    typedef typename UG_NS<dim>::template Entity<codim>::T UGEntity;

  public:

    typedef typename GridImp::template Codim<codim>::Entity Entity;
    constexpr static int codimension = codim;

    //! Constructor
    explicit UGGridLevelIterator() : gridImp_(nullptr)
    {
      entity_.impl().setToTarget(nullptr,nullptr);
    }

    //! Constructor
    explicit UGGridLevelIterator(const GridImp& gridImp, int level) : gridImp_(&gridImp)
    {
      auto& entity = entity_.impl();

      typename UG_NS<dim>::Grid *theGrid = const_cast<typename UG_NS<dim>::Grid* >(gridImp_->multigrid_->grids[level]);
      assert(theGrid);
      if (codim==dim) {
        if (pitype==All_Partition || pitype==Ghost_Partition)
          entity.setToTarget((UGEntity*)UG_NS<dim>::PFirstNode(theGrid),gridImp_);
        else
          entity.setToTarget((UGEntity*)UG_NS<dim>::FirstNode(theGrid),gridImp_);

      }
      else if (codim==0) {
        if (pitype==All_Partition || pitype==Ghost_Partition)
          entity.setToTarget((UGEntity*)UG_NS<dim>::PFirstElement(theGrid),gridImp_);
        else
          entity.setToTarget((UGEntity*)UG_NS<dim>::FirstElement(theGrid),gridImp_);
      }
      else
        DUNE_THROW(NotImplemented, "UGGrid level iterators for codimension " << codim);

      if (entity.getTarget() && !entityOK_())
        increment();
    }

    //! prefix increment
    void increment()
    {
      auto& entity = entity_.impl();
      assert(entity.level() == UG_NS<dim>::myLevel(entity.getTarget()));
      // Increment
      do {
        entity.setToTarget(UG_NS<dim>::succ(entity.getTarget()),gridImp_);
      }
      while (entity.getTarget() && !entityOK_());
    }

    //! dereferencing
    const Entity& dereference() const {return entity_;}

    //! equality
    bool equals(const UGGridLevelIterator<codim,pitype,GridImp>& other) const {
      return entity_ == other.entity_;
    }

  private:
    /**
     * \brief Return true iff the current entity is within the right
     *        partition.
     */
    bool entityOK_()
    {
      Dune::PartitionType entityPIType = entity_.impl().partitionType();
      switch (pitype) {
      case All_Partition:
        return true;
      case Ghost_Partition:
        if (entityPIType == GhostEntity)
          return true;
        else
          return false;
      case Interior_Partition:
        if (entityPIType == InteriorEntity)
          return true;
        else
          return false;
      case InteriorBorder_Partition:
      case Overlap_Partition:
      case OverlapFront_Partition:
        if (entityPIType == BorderEntity || entityPIType == InteriorEntity)
          return true;
        else
          return false;
      default:
        DUNE_THROW(NotImplemented, "Unhandled partition iterator type " << pitype);
      }
    }

    // /////////////////////////////////////
    //   Data members
    // /////////////////////////////////////
    const GridImp* gridImp_;

    //! The makeable entity that the iterator is pointing to
    Entity entity_;
  };

}  // namespace Dune

#endif
