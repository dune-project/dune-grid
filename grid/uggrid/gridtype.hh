// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGGRID_GRIDTYPE_HH
#define DUNE_UGGRID_GRIDTYPE_HH

#include <dune/grid/utility/griddim.hh>

#if defined UGGRID
  #if HAVE_GRIDTYPE
    #error "Ambiguous definition of GRIDTYPE."
  #endif
  #if not HAVE_UG
    #error "UGGRID defined but no UG version found!"
  #endif
  #if (GRIDDIM < 2) || (GRIDDIM > 3)
    #error "UGGRID is only available for GRIDDIM=2 and GRIDDIM=3."
  #endif
  #if (GRIDDIM != WORLDDIM)
    #error "UGGRID only supports GRIDDIM=WORLDDIM."
  #endif

  #include <dune/grid/uggrid.hh>
typedef Dune :: UGGrid< dimgrid > GridType;
  #define HAVE_GRIDTYPE 1
#endif

#endif
