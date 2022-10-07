// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IDENTITYGRIDLEVELITERATOR_HH
#define DUNE_IDENTITYGRIDLEVELITERATOR_HH

#include <dune/grid/common/gridenums.hh>

/** \file
 * \brief The IdentityGridLevelIterator class
 */

namespace Dune {

  /** \brief Iterator over all entities of a given codimension and level of a grid.
   * \ingroup IdentityGrid
   */
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class IdentityGridLevelIterator
  {

    typedef typename GridImp::HostGridType::Traits::template Codim<codim>::template Partition<pitype>::LevelIterator HostGridLevelIterator;

  public:

    constexpr static int codimension = codim;

    typedef typename GridImp::template Codim<codim>::Entity Entity;

    //! Constructor
    explicit IdentityGridLevelIterator(const GridImp* identityGrid, int level)
    : identityGrid_(identityGrid),
      hostLevelIterator_(identityGrid->hostgrid_->levelGridView(level).template begin<codim,pitype>())
    {}


    /** \brief Constructor which create the end iterator
        \param endDummy      Here only to distinguish it from the other constructor
        \param identityGrid  pointer to IdentityGrid instance
        \param level         grid level on which the iterator shall be created
     */
    explicit IdentityGridLevelIterator(const GridImp* identityGrid, int level, [[maybe_unused]] bool endDummy)
    : identityGrid_(identityGrid),
      hostLevelIterator_(identityGrid->hostgrid_->levelGridView(level).template end<codim,pitype>())
    {}


    //! prefix increment
    void increment() {
      ++hostLevelIterator_;
    }

    //! dereferencing
    Entity dereference() const {
      return Entity{{identityGrid_,*hostLevelIterator_}};
    }

    //! equality
    bool equals(const IdentityGridLevelIterator& i) const {
      return hostLevelIterator_ == i.hostLevelIterator_;
    }

  private:
    const GridImp* identityGrid_;

    HostGridLevelIterator hostLevelIterator_;
  };


}  // namespace Dune

#endif
