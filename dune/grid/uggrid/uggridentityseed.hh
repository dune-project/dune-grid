// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGGRID_ENTITY_SEED_HH
#define DUNE_UGGRID_ENTITY_SEED_HH

/** \file
 *  \brief Implementation of EntitySeed for the UGGrid grid manager
 */

namespace Dune {

  /** \brief Store a reference to an entity with a minimal memory footprint (one pointer)
   */
  template<int codim, class GridImp>
  class UGGridEntitySeed
  {
    // grid dimension
    constexpr static int dim = GridImp::dimension;
  public:

    //! codimension of underlying entity
    constexpr static int codimension = codim;

    //! default construct an invalid entity seed
    UGGridEntitySeed ()
      : target_(nullptr)
    {}

    //! construct entity seed from entity
    UGGridEntitySeed (const UGGridEntity<codim,dim,GridImp>& entity)
      : target_(entity.getTarget())
    {}

    //! check whether the EntitySeed refers to a valid Entity
    bool isValid() const
    {
      return target_ != nullptr;
    }

    /** \brief Access to the underlying UG data structure */
    typename UG_NS<dim>::template Entity<codim>::T* target() const
    {
      return target_;
    }

  private:
    /** \brief Plain old pointer to the corresponding UG data structure */
    typename UG_NS<dim>::template Entity<codim>::T* target_;
  };

} // end namespace Dune

#endif
