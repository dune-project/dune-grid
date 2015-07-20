// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_ENTITY_SEED_HH
#define DUNE_GRID_ENTITY_SEED_HH

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
    enum { codimension = EntitySeedImp::codimension };

    //! Export the implementation type
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

#if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
  public:
#else
  protected:
    // give the GridDefaultImplementation class access to the impl
    friend class GridDefaultImplementation<
        GridImp::dimension, GridImp::dimensionworld,
        typename GridImp::ctype,
        typename GridImp::GridFamily> ;
#endif

    /** \brief Access to the actual implementation */
    Implementation& impl()
    {
      return implementation_;
    }
    /** \brief const Access to the actual implementation */
    const Implementation& impl() const
    {
      return implementation_;
    }

  private:
    /** \brief The actual implementation class */
    EntitySeedImp implementation_;
  };

} // end namespace Dune

#endif
