// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_ENTITY_SEED_HH
#define DUNE_GRID_COMMON_ENTITY_SEED_HH

#include <dune/grid/common/grid.hh>

/** \file
 *  \brief Interface class EntitySeed
 */

namespace Dune {

  /** \brief Store a reference to an entity with a minimal memory footprint
   *
   * The EntitySeed provides a light-weight way to store an entity.  It is supposed
   * to be implemented as memory-efficiently as possible.  To get back the actual
   * entity, you need the corresponding grid.
   * On the grid, there is the method entity(const EntitySeed&), which
   * gives you an Entity in exchange for an EntitySeed.
   */
  template<class GridImp, class EntitySeedImp>
  class EntitySeed
  {
  public:

    //! codimension of underlying entity
    constexpr static int codimension = EntitySeedImp::codimension;

    /**
     * \brief type of underlying implementation
     *
     * \warning Implementation details may change without prior notification.
     **/
    typedef EntitySeedImp Implementation;

    /** \brief Construct an empty (i.e. isValid() == false) seed */
    EntitySeed()
    {}

    /** \brief Construct from implementation class */
    EntitySeed(const EntitySeedImp& implementation)
      : implementation_(implementation)
    {}

    /** \brief check whether it is safe to create an Entity from this Seed */
    bool isValid() const
    {
      return implementation_.isValid();
    }

    /**
     * \brief access to the underlying implementation
     *
     * \warning Implementation details may change without prior notification.
     **/
    Implementation& impl() { return implementation_; }

    /**
     * \brief access to the underlying implementation
     *
     * \warning Implementation details may change without prior notification.
     **/
    const Implementation& impl() const { return implementation_; }

  private:
    /** \brief The actual implementation class */
    EntitySeedImp implementation_;
  };

} // end namespace Dune

#endif
