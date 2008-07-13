// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_YASPGRID_GRIDTYPE_HH
#define DUNE_YASPGRID_GRIDTYPE_HH

#include <dune/grid/utility/griddim.hh>

#if defined YASPGRID
  #if HAVE_GRIDTYPE
    #error "Ambiguous definition of GRIDTYPE."
  #endif
  #if (GRIDDIM != WORLDDIM)
    #error "YASPGRID only supports GRIDDIM=WORLDDIM."
  #endif

  #include <dune/grid/yaspgrid.hh>
typedef Dune :: YaspGrid< dimgrid > GridType;
  #define HAVE_GRIDTYPE 1
#endif

#endif
