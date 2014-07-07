// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IDENTITYGRIDLEVELITERATOR_HH
#define DUNE_IDENTITYGRIDLEVELITERATOR_HH

#include "identitygridentitypointer.hh"

/** \file
 * \brief The IdentityGridLevelIterator class and its specializations
 */

namespace Dune {




  //**********************************************************************
  //
  // --IdentityGridLevelIterator
  /** \brief Iterator over all entities of a given codimension and level of a grid.
   * \ingroup IdentityGrid
   */
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class IdentityGridLevelIterator :
    public Dune::IdentityGridEntityPointer <codim,GridImp>
  {
  public:

    //! Constructor
    explicit IdentityGridLevelIterator(const GridImp* identityGrid, int level)
      : IdentityGridEntityPointer<codim,GridImp>(identityGrid, identityGrid->hostgrid_->template lbegin<codim>(level)),
        hostGridLevelIterator_(identityGrid->hostgrid_->template lbegin<codim>(level))
    {
      this->virtualEntity_.setToTarget(hostGridLevelIterator_);
    }


    /** \brief Constructor which create the end iterator
        \param endDummy      Here only to distinguish it from the other constructor
        \param identityGrid  pointer to IdentityGrid instance
        \param level         grid level on which the iterator shall be created
     */
    explicit IdentityGridLevelIterator(const GridImp* identityGrid, int level, bool endDummy)
      :
        IdentityGridEntityPointer<codim,GridImp>(identityGrid, identityGrid->hostgrid_->template lend<codim>(level)),
        hostGridLevelIterator_(identityGrid->hostgrid_->template lend<codim>(level))
    {}


    //! prefix increment
    void increment() {
      ++hostGridLevelIterator_;
      this->virtualEntity_.setToTarget(hostGridLevelIterator_);
    }


  private:

    // LevelIterator to the equivalent entity in the host grid
    typedef typename GridImp::HostGridType::Traits::template Codim<codim>::LevelIterator HostGridLevelIterator;

    //! \todo Please doc me !
    HostGridLevelIterator hostGridLevelIterator_;

  };


}  // namespace Dune

#endif
