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
    public Dune::IdentityGridEntityPointer<codim,GridImp,typename GridImp::HostGridType::template Codim<codim>::template Partition<pitype>::LeafIterator>
  {
  private:

    // LevelIterator to the equivalent entity in the host grid
    typedef typename GridImp::HostGridType::template Codim<codim>::template Partition<pitype>::LeafIterator HostGridLeafIterator;

    typedef Dune::IdentityGridEntityPointer<codim,GridImp,HostGridLeafIterator> Base;

  public:

    //! \todo Please doc me !
    explicit IdentityGridLeafIterator(const GridImp* identityGrid) :
      Base(identityGrid, identityGrid->hostgrid_->leafGridView().template begin<codim,pitype>())
    {}

    /** \brief Constructor which create the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     *  \param identityGrid  pointer to grid instance
     */
    explicit IdentityGridLeafIterator(const GridImp* identityGrid, bool endDummy) :
      Base(identityGrid, identityGrid->hostgrid_->leafGridView().template end<codim,pitype>())
    {}


    //! prefix increment
    void increment() {
      ++this->hostEntityPointer_;
    }

  };


}  // namespace Dune

#endif
