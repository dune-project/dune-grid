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
    enum {dim = GridImp::dimension};

    friend class UGGridEntity<codim,GridImp::dimension,GridImp>;
    friend class UGGridEntity<0,    GridImp::dimension,GridImp>;

    // The type of the UG entity we're pointing to
    typedef typename UG_NS<dim>::template Entity<codim>::T UGEntity;

  public:

    typedef typename GridImp::template Codim<codim>::Entity Entity;
    enum {codimension = codim};

    //! Constructor
    explicit UGGridLevelIterator()
    {
      virtualEntity_.setToTarget(nullptr,nullptr);
    }

    //! Constructor
    explicit UGGridLevelIterator(const GridImp& gridImp, int level) : gridImp_(&gridImp)
    {
      typename UG_NS<dim>::Grid *theGrid = const_cast<typename UG_NS<dim>::Grid* >(gridImp_->multigrid_->grids[level]);
      assert(theGrid);
      if (codim==dim) {
        if (pitype==All_Partition || pitype==Ghost_Partition)
          virtualEntity_.setToTarget((UGEntity*)UG_NS<dim>::PFirstNode(theGrid),gridImp_);
        else
          virtualEntity_.setToTarget((UGEntity*)UG_NS<dim>::FirstNode(theGrid),gridImp_);

      }
      else if (codim==0) {
        if (pitype==All_Partition || pitype==Ghost_Partition)
          virtualEntity_.setToTarget((UGEntity*)UG_NS<dim>::PFirstElement(theGrid),gridImp_);
        else
          virtualEntity_.setToTarget((UGEntity*)UG_NS<dim>::FirstElement(theGrid),gridImp_);
      }
      else
        DUNE_THROW(NotImplemented, "UGGrid level iterators for codimension " << codim);

      if (virtualEntity_.getTarget() && !entityOK_())
        increment();
    }

    //! prefix increment
    void increment()
    {
      assert(virtualEntity_.level() == UG_NS<dim>::myLevel(virtualEntity_.getTarget()));
      // Increment
      do {
        virtualEntity_.setToTarget(UG_NS<dim>::succ(virtualEntity_.getTarget()),gridImp_);
      }
      while (virtualEntity_.getTarget() && !entityOK_());
    }

    //! dereferencing
    const Entity& dereference() const {return virtualEntity_;}

    //! equality
    bool equals(const UGGridLevelIterator<codim,pitype,GridImp>& other) const {
      return virtualEntity_ == other.virtualEntity_;
    }

  private:
    /**
     * \brief Return true iff the current entity is within the right
     *        partition.
     */
    bool entityOK_()
    {
      Dune::PartitionType entityPIType = virtualEntity_.partitionType();
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
    UGMakeableEntity<codim,dim,GridImp> virtualEntity_;
  };

}  // namespace Dune

#endif
