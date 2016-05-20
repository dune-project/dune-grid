// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ONE_D_GRID_ENTITY_SEED_HH
#define DUNE_ONE_D_GRID_ENTITY_SEED_HH

#include "onedgridentity.hh"

/** \file
 *  \brief Implementation of EntitySeed for the OneDGrid grid manager
 */

namespace Dune {

  /** \brief Store a reference to an entity with a minimal memory footprint
   */
  template<int codim, class GridImp>
  class OneDGridEntitySeed
  {
    // grid dimension
    enum { dim = GridImp::dimension };
  public:

    enum { codimension = codim };

    //! default construct an invalid entity seed
    OneDGridEntitySeed ()
      : target_(nullptr)
    {}

    //! construct entity seed from entity
    OneDGridEntitySeed (const OneDGridEntity<codim,dim,GridImp>& entity)
      : target_(entity.target_)
    {}

    //! check whether the EntitySeed refers to a valid Entity
    bool isValid() const
    {
      return target_ != nullptr;
    }

    /** \brief Access to the underlying OneDGrid data structure */
    OneDEntityImp<dim-codim>* target() const
    {
      return target_;
    }

  private:
    /** \brief Plain old pointer to the corresponding OneDGrid data structure */
    OneDEntityImp<dim-codim>* target_;
  };

} // end namespace Dune

#endif
