// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGGRIDLEAFITERATOR_HH
#define DUNE_UGGRIDLEAFITERATOR_HH

/** \file
 * \brief The UGGridLeafIterator class
 */

#include <dune/grid/uggrid/uggridentity.hh>

namespace Dune {

  /** \brief Iterator over all entities of a given codimension and level of a grid.
   * \ingroup UGGrid
   */
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class UGGridLeafIterator
  {
    constexpr static int dim = GridImp::dimension;

    // The type of the UG entity we're pointing to
    typedef typename UG_NS<dim>::template Entity<codim>::T UGEntity;

  public:

    typedef typename GridImp::template Codim<codim>::Entity Entity;
    constexpr static int codimension = codim;

    /** \brief Constructor setting up a 'begin' iterator
     */
    UGGridLeafIterator(const GridImp& grid) : grid_(&grid) {
      // Entities below this level are certainly not leaf entities
      int levelCounter = grid.leafIndexSet_.coarsestLevelWithLeafElements_;

      // If the grid is distributed, the grid on the 'coarsestLevelWithLeafElements_' may actually be empty, because
      // it is all on other processors.  Therefore we have to also look at finer levels.
      // In a sequential program, the while loops always iterate exactly once.
      auto& entity = entity_.impl();
      if (codim==dim) {
        if (pitype==All_Partition || pitype==Ghost_Partition) {
          do {
            entity.setToTarget((UGEntity*)UG_NS<dim>::PFirstNode(grid_->multigrid_->grids[levelCounter++]), grid_);
          } while (not entity.getTarget() and levelCounter <= grid.maxLevel());
        } else {
          do {
            entity.setToTarget((UGEntity*)UG_NS<dim>::FirstNode(grid_->multigrid_->grids[levelCounter++]), grid_);
          } while (not entity.getTarget() and levelCounter <= grid.maxLevel());
        }

      } else if (codim==0) {
        if (pitype==All_Partition || pitype==Ghost_Partition) {
          do {
            entity.setToTarget((UGEntity*)UG_NS<dim>::PFirstElement(grid_->multigrid_->grids[levelCounter++]), grid_);
          } while (not entity.getTarget() and levelCounter <= grid.maxLevel());
        } else {
          do {
            entity.setToTarget((UGEntity*)UG_NS<dim>::FirstElement(grid_->multigrid_->grids[levelCounter++]), grid_);
          } while (not entity.getTarget() and levelCounter <= grid.maxLevel());
        }

      } else
        DUNE_THROW(NotImplemented, "UGGrid leaf iterators for codimension " << codim);

      if (entity.getTarget() && ! (isLeaf() && isInPartition()))
        increment();
    }

    /** \brief Constructor setting up an 'end' iterator
     */
    UGGridLeafIterator() : grid_(nullptr)
    {
      entity_.impl().setToTarget(nullptr,nullptr);
    }

    //! prefix increment
    void increment() {
      // Increment until you find a leaf entity
      do {
        globalIncrement();
      }
      while (entity_.impl().getTarget() && ! (isLeaf() && isInPartition()));
    }

    //! dereferencing
    const Entity& dereference() const {return entity_;}

    //! equality
    bool equals(const UGGridLeafIterator<codim,pitype,GridImp>& other) const {
      return entity_ == other.entity_;
    }

  private:
    /**
     * \brief Return true iff the current entity is a leaf (element version)
     */
    bool isLeaf(const typename UG_NS<dim>::Element* theElement)
    {
      return UG_NS<dim>::isLeaf(theElement);
    }

    /**
     * \brief Return true iff the current entity is a leaf (node version)
     */
    bool isLeaf(const typename UG_NS<dim>::Node* theNode)
    {
      // in the case of nodes: only visit the uppermost entity if more than one is leaf
      if (theNode->son != nullptr)
        return false;

      return UG_NS<dim>::isLeaf(theNode);
    }

    /**
     * \brief Return true iff the current entity has the right PartitionType.
     */
    bool isInPartition()
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

    /**
     * \brief Return true iff the current entity is a leaf entity
     */
    bool isLeaf()
    {
      return isLeaf(entity_.impl().getTarget());
    }

    /** \brief This increment makes the iterator wander over all entities on all levels */
    void globalIncrement() {
      auto& entity = entity_.impl();

      // Store the current level, so we know it even when pointing to an invalid entity
      int level = entity.level();

      // Increment on this level
      entity.setToTarget(UG_NS<dim>::succ(entity.getTarget()), grid_);

      // If beyond the end of this level set to first of next level
      if (!entity.getTarget() && level < grid_->maxLevel()) {

        if (codim==dim) {

          if (pitype==All_Partition || pitype==Ghost_Partition)
            entity.setToTarget((UGEntity*)UG_NS<dim>::PFirstNode(grid_->multigrid_->grids[level+1]), grid_);
          else
            entity.setToTarget((UGEntity*)UG_NS<dim>::FirstNode(grid_->multigrid_->grids[level+1]), grid_);

        } else if (codim==0) {

          if (pitype==All_Partition || pitype==Ghost_Partition)
            entity.setToTarget((UGEntity*)UG_NS<dim>::PFirstElement(grid_->multigrid_->grids[level+1]), grid_);
          else
            entity.setToTarget((UGEntity*)UG_NS<dim>::FirstElement(grid_->multigrid_->grids[level+1]), grid_);

        } else
          DUNE_THROW(NotImplemented, "UGGrid leaf iterators for codimension " << codim);

      }

    }

    // /////////////////////////////////////
    //   Data members
    // /////////////////////////////////////
    const GridImp* grid_;

    //! The makeable entity that the iterator is pointing to
    Entity entity_;
  };

}  // namespace Dune

#endif
