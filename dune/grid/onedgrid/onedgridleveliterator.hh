// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ONE_D_GRID_LEVELITERATOR_HH
#define DUNE_ONE_D_GRID_LEVELITERATOR_HH

/** \file
 * \brief The OneDGridLevelIterator class
 */

#include <dune/grid/common/gridenums.hh>

namespace Dune {



  //**********************************************************************
  //
  // --OneDGridLevelIterator
  // --LevelIterator
  /** \brief Iterator over all entities of a given codimension and level of a grid.
   * \ingroup OneDGrid
   */
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class OneDGridLevelIterator
  {
  public:
    constexpr static int dim = GridImp::dimension;
    friend class OneDGrid;
    friend class OneDGridEntity<codim,dim,GridImp>;
    friend class OneDGridEntity<0,dim,GridImp>;
    friend class OneDGridLevelGridView<GridImp>;

    typedef typename GridImp::template Codim<codim>::Entity Entity;
    constexpr static int codimension = codim;

  protected:

    /** \brief Constructor from a given iterator */
    OneDGridLevelIterator(OneDEntityImp<dim-codim>* it)
    {
      virtualEntity_.impl().setToTarget(it);
    }

  public:

    //! prefix increment
    void increment() {
      this->virtualEntity_.impl().setToTarget(this->virtualEntity_.impl().target_->succ_);
    }

    //! dereferencing
    const Entity& dereference() const {return virtualEntity_;}

    //! equality
    bool equals(const OneDGridLevelIterator<codim,pitype,GridImp>& other) const {
      return virtualEntity_ == other.virtualEntity_;
    }

  protected:

    //! The entity that the iterator is pointing to
    Entity virtualEntity_;
  };

}  // namespace Dune

#endif
