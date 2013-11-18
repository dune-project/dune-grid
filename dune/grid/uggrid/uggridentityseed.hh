// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGGRID_ENTITY_SEED_HH
#define DUNE_UGGRID_ENTITY_SEED_HH

/** \file
 *  \brief Implementation of EntitySeed for the UGGrid grid manager
 */

#include <dune/common/nullptr.hh>

namespace Dune {

  /** \brief Store a reference to an entity with a minimal memory footprint (one pointer)
   */
  template<int codim, class GridImp>
  class UGGridEntitySeed
  {
    // grid dimension
    enum { dim = GridImp::dimension };
  public:

    //! codimension of underlying entity
    enum { codimension = codim };

    //! default construct an invalid entity seed
    UGGridEntitySeed ()
      : target_(nullptr)
    {}

    //! construct entity seed from entity
    UGGridEntitySeed (const UGGridEntity<codim,dim,GridImp>& entity)
      : target_(entity.target_)
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
