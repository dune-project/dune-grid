// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGGRIDLEAFITERATOR_HH
#define DUNE_UGGRIDLEAFITERATOR_HH

#include <dune/grid/uggrid/uggridentitypointer.hh>

/** \file
 * \brief The UGGridLeafIterator class
 */

namespace Dune {

  /** \brief Iterator over all entities of a given codimension and level of a grid.
   * \ingroup UGGrid
   */
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class UGGridLeafIterator :
    public Dune::UGGridEntityPointer <codim,GridImp>
  {
    enum {dim = GridImp::dimension};

    // The type of the UG entity we're pointing to
    typedef typename UG_NS<dim>::template Entity<codim>::T UGEntity;

  public:

    UGGridLeafIterator(const GridImp& grid) : grid_(&grid) {
      // Entities below this level are certainly not leaf entities
      const unsigned int startingLevel = grid.leafIndexSet_.coarsestLevelWithLeafElements_;

      if (codim==dim) {
        if (pitype==All_Partition || pitype==Ghost_Partition)
          setToTarget((UGEntity*)UG_NS<dim>::PFirstNode(grid_->multigrid_->grids[startingLevel]));
        else if (pitype == Dune::Interior_Partition || pitype == Dune::InteriorBorder_Partition)
          setToTarget((UGEntity*)UG_NS<dim>::FirstNode(grid_->multigrid_->grids[startingLevel]));
        else     // overlap and overlap-front
          setToTarget((UGEntity*) 0);

      } else if (codim==0) {
        if (pitype==All_Partition || pitype==Ghost_Partition)
          setToTarget((UGEntity*)UG_NS<dim>::PFirstElement(grid_->multigrid_->grids[startingLevel]));
        else if (pitype == Dune::Interior_Partition || pitype == Dune::InteriorBorder_Partition)
          setToTarget((UGEntity*)UG_NS<dim>::FirstElement(grid_->multigrid_->grids[startingLevel]));
        else     // overlap and overlap-front
          setToTarget((UGEntity*) 0);

      } else
        DUNE_THROW(NotImplemented, "UGGrid leaf iterators for codimension " << codim);

      if (this->virtualEntity_.getTarget() && !entityOK_())
        increment();
    }

    //! Constructor
    UGGridLeafIterator()
    {
      this->virtualEntity_.setToTarget(0);
    }

    //! prefix increment
    void increment() {
      // Increment until you find a leaf entity
      do {
        globalIncrement();
      }
      while (this->virtualEntity_.getTarget() && !entityOK_());
    }

  private:
    /**
     * \brief Return true iff the current entity is a leaf and is within the right partition.
     */
    bool entityOK_()
    {
      Dune::PartitionType entityPIType = this->virtualEntity_.partitionType();
      // HACK: ghosts entities are never leafs in UG, in DUNE can be
      //       leafs
      if (entityPIType != GhostEntity &&
          !UG_NS<dim>::isLeaf(this->virtualEntity_.getTarget()))
      {
        return false;
      }

      if (pitype == All_Partition)
        return true;

      if (pitype == Ghost_Partition && entityPIType == GhostEntity)
        return true;
      else if (pitype == Interior_Partition && entityPIType == InteriorEntity)
        return true;
      else if (pitype == InteriorBorder_Partition &&
               (entityPIType == BorderEntity ||
                entityPIType == InteriorEntity))
        return true;
      return false;
    }

    /** \brief This increment makes the iterator wander over all entities on all levels */
    void globalIncrement() {

      // Store the current level, so we know it even when pointing to an invalid entity
      int level = this->level();

      // Increment on this level
      this->virtualEntity_.setToTarget(UG_NS<dim>::succ(this->virtualEntity_.getTarget()));

      // If beyond the end of this level set to first of next level
      if (!this->virtualEntity_.getTarget() && level < grid_->maxLevel()) {

        if (codim==dim) {

          if (pitype==All_Partition || pitype==Ghost_Partition)
            setToTarget((UGEntity*)UG_NS<dim>::PFirstNode(grid_->multigrid_->grids[level+1]));
          else
            setToTarget((UGEntity*)UG_NS<dim>::FirstNode(grid_->multigrid_->grids[level+1]));

        } else if (codim==0) {

          if (pitype==All_Partition || pitype==Ghost_Partition)
            setToTarget((UGEntity*)UG_NS<dim>::PFirstElement(grid_->multigrid_->grids[level+1]));
          else
            setToTarget((UGEntity*)UG_NS<dim>::FirstElement(grid_->multigrid_->grids[level+1]));

        } else
          DUNE_THROW(NotImplemented, "UGGrid leaf iterators for codimension " << codim);

      }

    }

    // /////////////////////////////////////
    //   Data members
    // /////////////////////////////////////
    const GridImp* grid_;

  };

}  // namespace Dune

#endif
