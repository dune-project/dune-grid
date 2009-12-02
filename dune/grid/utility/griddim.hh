// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRIDDIM_HH
#define DUNE_GRIDDIM_HH

/**
 * @file
 * @brief Determine default grid and world dimension.
 */

#ifndef GRIDDIM
  #warning "GRIDDIM not defined, defaulting to GRIDDIM=3"
//! Dimension of grid
/**
 * Set with <tt>configure --with-grid-dim</tt>.  If unset, defaults to 3 and
 * issues a warning.
 *
 * \note This is in general not "the" dimension of the grid, it only
 *       provides a default for certain cases.
 */
  #define GRIDDIM 3
#endif
#if not (GRIDDIM > 0)
  #error "GRIDDIM must be a positive integer."
#endif

#ifndef WORLDDIM
//! Dimension of world
/**
 * Set with <tt>configure --with-world-dim</tt>.  If unset, defaults to \ref
 * GRIDDIM.
 *
 * \note This is in general not "the" dimension of the world, it only
 *       provides a default for certain cases.
 */
  #define WORLDDIM GRIDDIM
#endif
#if (WORLDDIM < GRIDDIM)
  #error "WORLDDIM < GRIDDIM does not make sense."
#endif

//! Dimension of grid
/**
 * This just stores \ref GRIDDIM in a C++ constant.
 *
 * \note This is in general not "the" dimension of the grid, it only provides
 *       a default for certain cases.
 */
const int dimgrid = GRIDDIM;
//! Dimension of world
/**
 * This just stores \ref WORLDDIM in a C++ constant.
 *
 * \note This is in general not "the" dimension of the world, it only provides
 *       a default for certain cases.
 */
const int dimworld = WORLDDIM;

#endif
