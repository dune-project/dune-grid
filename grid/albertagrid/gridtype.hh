// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTAGRID_GRIDTYPE_HH
#define DUNE_ALBERTAGRID_GRIDTYPE_HH

#include <dune/grid/utility/griddim.hh>

#if defined ALBERTAGRID
  #if HAVE_GRIDTYPE
    #error "Ambiguous definition of GRIDTYPE."
  #endif
  #if not HAVE_ALBERTA
    #error "ALBERTAGRID defined but no ALBERTA version found!"
  #endif
  #if (GRIDDIM < 2) || (GRIDDIM > 3)
    #error "ALBERTAGRID is only available for GRIDDIM=2 and GRIDDIM=3."
  #endif
  #if (WORLDDIM != GRIDDIM)
    #error "ALBERTAGRID currently only supports WORLDDIM=GRIDDIM."
  #endif

  #include <dune/grid/albertagrid.hh>
typedef Dune :: AlbertaGrid< dimgrid, dimworld > GridType;
  #define HAVE_GRIDTYPE 1
#endif

#endif
