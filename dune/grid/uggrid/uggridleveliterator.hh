// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGGRIDLEVELITERATOR_HH
#define DUNE_UGGRIDLEVELITERATOR_HH

#include <dune/grid/uggrid/uggridentitypointer.hh>

/** \file
 * \brief The UGGridLevelIterator class
 */

namespace Dune {

  //**********************************************************************
  //
  // --UGGridLevelIterator
  // --LevelIterator
  /** \brief Iterator over all entities of a given codimension and level of a grid.
   * \ingroup UGGrid
   */
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class UGGridLevelIterator :
    public Dune::UGGridEntityPointer <codim,GridImp>
  {
    enum {dim = GridImp::dimension};

    friend class UGGridEntity<codim,GridImp::dimension,GridImp>;
    friend class UGGridEntity<0,    GridImp::dimension,GridImp>;

    // The type of the UG entity we're pointing to
    typedef typename UG_NS<dim>::template Entity<codim>::T UGEntity;

  public:

    typedef typename GridImp::template Codim<codim>::Entity Entity;

    //! Constructor
    explicit UGGridLevelIterator()
    {
      this->virtualEntity_.setToTarget(nullptr,nullptr);
    }

    //! Constructor
    explicit UGGridLevelIterator(const GridImp& gridImp, int level) : gridImp_(&gridImp)
    {
      typename UG_NS<dim>::Grid *theGrid = const_cast<typename UG_NS<dim>::Grid* >(gridImp_->multigrid_->grids[level]);
      assert(theGrid);
      if (codim==dim) {
        if (pitype==All_Partition || pitype==Ghost_Partition)
          this->virtualEntity_.setToTarget((UGEntity*)UG_NS<dim>::PFirstNode(theGrid),gridImp_);
        else if (pitype == Dune::Interior_Partition || pitype == Dune::InteriorBorder_Partition)
          this->virtualEntity_.setToTarget((UGEntity*)UG_NS<dim>::FirstNode(theGrid),gridImp_);
        else     // overlap and overlap-front
          this->virtualEntity_.setToTarget(0,nullptr);

      }
      else if (codim==0) {
        if (pitype==All_Partition || pitype==Ghost_Partition)
          this->virtualEntity_.setToTarget((UGEntity*)UG_NS<dim>::PFirstElement(theGrid),gridImp_);
        else if (pitype == Dune::Interior_Partition || pitype == Dune::InteriorBorder_Partition)
          this->virtualEntity_.setToTarget((UGEntity*)UG_NS<dim>::FirstElement(theGrid),gridImp_);
        else     // overlap and overlap-front
          this->virtualEntity_.setToTarget(0,nullptr);
      }
      else
        DUNE_THROW(NotImplemented, "UGGrid leaf iterators for codimension " << codim);

      if (this->virtualEntity_.getTarget() && !entityOK_())
        increment();
    }

    //! prefix increment
    void increment()
    {
      assert(this->level() == UG_NS<dim>::myLevel(this->virtualEntity_.getTarget()));
      // Increment
      do {
        this->virtualEntity_.setToTarget(UG_NS<dim>::succ(this->virtualEntity_.getTarget()),gridImp_);
      }
      while (this->virtualEntity_.getTarget() && !entityOK_());
    }

  private:
    /**
     * \brief Return true iff the current entity is within the right
     *        partition.
     */
    bool entityOK_()
    {
      if (pitype == All_Partition)
        return true;

      Dune::PartitionType entityPIType = this->virtualEntity_.partitionType();
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

    // /////////////////////////////////////
    //   Data members
    // /////////////////////////////////////
    const GridImp* gridImp_;

  };

}  // namespace Dune

#endif
