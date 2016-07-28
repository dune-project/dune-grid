// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ONE_D_GRID_LEAFITERATOR_HH
#define DUNE_ONE_D_GRID_LEAFITERATOR_HH

/** \file
 * \brief The OneDGridLeafIterator class
 */

namespace Dune {

  /** \brief Iterator over all entities of a given codimension and level of a grid.
   * \ingroup OneDGrid
   */
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class OneDGridLeafIterator
  {
    enum {dim = GridImp::dimension};

    friend class OneDGridEntity<codim,dim,GridImp>;

  public:

    typedef typename GridImp::template Codim<codim>::Entity Entity;
    enum {codimension = codim};

    OneDGridLeafIterator(const GridImp& grid) : grid_(&grid) {

      /** \todo Can a make the fullRefineLevel work somehow? */
      const int fullRefineLevel = 0;

      GridImp::getRealImplementation(this->virtualEntity_).setToTarget((OneDEntityImp<1-codim>*) std::get<1-codim>(grid_->entityImps_[fullRefineLevel]).begin());

      if (!GridImp::getRealImplementation(this->virtualEntity_).target_->isLeaf())
        increment();
    }

    //! Constructor
    OneDGridLeafIterator()
    {
      GridImp::getRealImplementation(this->virtualEntity_).setToTarget(OneDGridNullIteratorFactory<1-codim>::null());
    }

    //! prefix increment
    void increment() {
      // Increment until you find a leaf entity
      do {
        globalIncrement();
      } while (GridImp::getRealImplementation(this->virtualEntity_).target_
               && !GridImp::getRealImplementation(this->virtualEntity_).target_->isLeaf());
    }

    //! dereferencing
    const Entity& dereference() const {return virtualEntity_;}

    //! equality
    bool equals(const OneDGridLeafIterator<codim,pitype,GridImp>& other) const {
      return virtualEntity_ == other.virtualEntity_;
    }

  private:

    /** \brief This increment makes the iterator wander over all entities on all levels */
    void globalIncrement() {

      // Backup current level because it may not be accessible anymore after
      // setting the pointer to the next entity.
      const int oldLevel = this->virtualEntity_.level();

      // Increment on this level
      GridImp::getRealImplementation(this->virtualEntity_).setToTarget(GridImp::getRealImplementation(this->virtualEntity_).target_->succ_);

      // If beyond the end of this level set to first of next level
      if (!GridImp::getRealImplementation(this->virtualEntity_).target_ && oldLevel < grid_->maxLevel()) {

        GridImp::getRealImplementation(this->virtualEntity_).setToTarget(const_cast<OneDEntityImp<dim-codim>*>(std::get<1-codim>(grid_->entityImps_[oldLevel+1]).begin()));

      }

    }

    // /////////////////////////////////////
    //   Data members
    // /////////////////////////////////////
    const GridImp* grid_;

    //! The entity that the iterator is pointing to
    Entity virtualEntity_;
  };

}  // namespace Dune

#endif
