// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IDENTITYGRIDLEAFITERATOR_HH
#define DUNE_IDENTITYGRIDLEAFITERATOR_HH

#include <dune/grid/common/gridenums.hh>

/** \file
 * \brief The IdentityGridLeafIterator class
 */

namespace Dune {


  /** \brief Iterator over all entities of a given codimension and level of a grid.
   *  \ingroup IdentityGrid
   */
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class IdentityGridLeafIterator
  {
  private:

    // LevelIterator to the equivalent entity in the host grid
    typedef typename GridImp::HostGridType::template Codim<codim>::template Partition<pitype>::LeafIterator HostGridLeafIterator;

  public:

    constexpr static int codimension = codim;

    typedef typename GridImp::template Codim<codim>::Entity Entity;

    //! \todo Please doc me !
    explicit IdentityGridLeafIterator(const GridImp* identityGrid) :
      identityGrid_(identityGrid),
      hostLeafIterator_(identityGrid->hostgrid_->leafGridView().template begin<codim,pitype>())
    {}

    /** \brief Constructor which create the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     *  \param identityGrid  pointer to grid instance
     */
    explicit IdentityGridLeafIterator(const GridImp* identityGrid, [[maybe_unused]] bool endDummy) :
      identityGrid_(identityGrid),
      hostLeafIterator_(identityGrid->hostgrid_->leafGridView().template end<codim,pitype>())
    {}


    //! prefix increment
    void increment() {
      ++hostLeafIterator_;
    }

    //! dereferencing
    Entity dereference() const {
      return Entity{{identityGrid_,*hostLeafIterator_}};
    }

    //! equality
    bool equals(const IdentityGridLeafIterator& i) const {
      return hostLeafIterator_ == i.hostLeafIterator_;
    }

  private:
    const GridImp* identityGrid_;

    HostGridLeafIterator hostLeafIterator_;

  };


}  // namespace Dune

#endif
