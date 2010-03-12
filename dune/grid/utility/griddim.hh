// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRIDDIM_HH
#define DUNE_GRIDDIM_HH

/**
 * @file
 * @brief Determine default grid and world dimension.
 *
 * Set integer constants during make process which can be used to construct grids
 * of given domain and coordinate dimensions. To set these constants use
 * <tt> GRIDDIM=d WORDDIM=w </tt> to the make command.
 * A default value can be set with during configure using
 * <tt>--with-grid-dim</tt> and
 * <tt>--with-world-dim</tt>. If no default is set and no
 * value is provided during make, an compile time error is generated.
 * The constants are placed in the namespace Dune::GridSelector.
 * This construction only works if the <tt> GRIDDIM_CPPFLAGS </tt> or the
 * <tt> ALL_PKG_CPPFLAGS </tt> are used.
 *
 * \note This is in general not "the" dimension of the grid or world, it only
 *       provides a default for certain cases and can be used together
 *       with the definition of a GridTtype through gridtype.hh.
 *
 */

#include <dune/common/deprecated.hh>

// first check if the constants have been set
#ifndef GRIDDIM
  #error "GRIDDIM not defined, use 'ALL_PKG_CPPFLAGS'."
#endif
#if not (GRIDDIM > 0)
  #if HEADERCHECK
    #undef GRIDDIM
    #define GRIDDIM 2
    #undef WORLDDIM
    #define WORLDDIM 2
  #else
    #error "GRIDDIM must be a positive integer. Specify GRIDDIM in make command."
  #endif
#endif

#ifndef WORLDDIM
  #define WORLDDIM GRIDDIM
#endif
#if not (WORLDDIM >= GRIDDIM)
  #error "WORLDDIM < GRIDDIM does not make sense."
#endif

namespace Dune
{
  namespace GridSelector
  {
    /**
     * @brief Dimension of grid
     *
     * This just stores \ref GRIDDIM in a C++ constant in the namespace
     * Dune::GridSelector.
     *
     * \note This is in general not "the" dimension of the grid, it only provides
     *       a default for certain cases.
     */
    const int dimgrid = GRIDDIM;
    /**
     * @brief Dimension of world
     *
     * This just stores \ref WORLDDIM in a C++ constant in the namespace
     * Dune::GridSelector.
     *
     * \note This is in general not "the" dimension of the world, it only provides
     *       a default for certain cases.
     */
    const int dimworld = WORLDDIM;
  }
}

const int DUNE_DEPRECATED dimgrid  = Dune::GridSelector::dimgrid;
const int DUNE_DEPRECATED dimworld = Dune::GridSelector::dimworld;

#endif
