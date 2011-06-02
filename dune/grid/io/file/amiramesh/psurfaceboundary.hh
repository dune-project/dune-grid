// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef PSURFACE_BOUNDARY_HH
#define PSURFACE_BOUNDARY_HH

/** \file
 *  \brief A domain boundary implemented by the psurface library
 */

#include <psurface/PSurface.h>

namespace Dune {

  /** \brief A domain boundary implemented by the psurface library
   *
   * \warning This code is experimental.  It may change in all kinds of unexpected ways
   * without prior notice.  Use it only if you know what you are doing.
   *
   * \tparam dim The dimension of the <b>boundary</b>, not the dimension of the grid.
   */
  template <int dim>
  class PSurfaceBoundary
  {
    dune_static_assert((dim==1 or dim==2), "PSurfaceBoundaries can only have dimensions 1 or 2!");

  public:

    /** \brief Constructor from a given PSurface object */
    PSurfaceBoundary(PSurface<dim,float>* psurface)
      : psurface_(psurface)
    {}

  private:

    std::auto_ptr<PSurface<dim,float> > psurface_;

  };

}

#endif
