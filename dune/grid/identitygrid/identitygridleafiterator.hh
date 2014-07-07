// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IDENTITYGRIDLEAFITERATOR_HH
#define DUNE_IDENTITYGRIDLEAFITERATOR_HH

#include "identitygridentitypointer.hh"

/** \file
 * \brief The IdentityGridLeafIterator class
 */

namespace Dune {


  /** \brief Iterator over all entities of a given codimension and level of a grid.
   *  \ingroup IdentityGrid
   */
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class IdentityGridLeafIterator :
    public Dune::IdentityGridEntityPointer <codim,GridImp>
  {
  private:

    enum {dim = GridImp::dimension};


  public:

    //! \todo Please doc me !
    explicit IdentityGridLeafIterator(const GridImp* identityGrid) :
      IdentityGridEntityPointer<codim,GridImp>(identityGrid, identityGrid->hostgrid_->template leafbegin<codim>()),
      hostGridLeafIterator_(identityGrid->hostgrid_->template leafbegin<codim>()),
      hostGridLeafEndIterator_(identityGrid->hostgrid_->template leafend<codim>())
    {
      this->virtualEntity_.setToTarget(hostGridLeafIterator_);
    }


    /** \brief Constructor which create the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     *  \param identityGrid  pointer to grid instance
     */
    explicit IdentityGridLeafIterator(const GridImp* identityGrid, bool endDummy) :
      IdentityGridEntityPointer<codim,GridImp>(identityGrid, identityGrid->hostgrid_->template leafend<codim>()),
      hostGridLeafIterator_(identityGrid->hostgrid_->template leafbegin<codim>()),
      hostGridLeafEndIterator_(identityGrid->hostgrid_->template leafend<codim>())
    {}


    //! prefix increment
    void increment() {
      ++hostGridLeafIterator_;
      this->virtualEntity_.setToTarget(hostGridLeafIterator_);
    }


  private:

    // /////////////////////////////////////
    //   Data members
    // /////////////////////////////////////

    // LevelIterator to the equivalent entity in the host grid
    typedef typename GridImp::HostGridType::template Codim<codim>::LeafIterator HostGridLeafIterator;

    //! \todo Please doc me !
    HostGridLeafIterator hostGridLeafIterator_;

    //! \todo Please doc me !
    HostGridLeafIterator hostGridLeafEndIterator_;

  };


}  // namespace Dune

#endif
