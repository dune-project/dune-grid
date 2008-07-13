// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_SGRID_GRIDTYPE_HH
#define DUNE_SGRID_GRIDTYPE_HH

#include <dune/grid/utility/griddim.hh>

#if defined SGRID
  #if HAVE_GRIDTYPE
    #error "Ambiguous definition of GRIDTYPE."
  #endif
  #if (GRIDDIM != WORLDDIM)
    #error "SGRID is only available for GRIDDIM=WORLDDIM."
  #endif

  #include <dune/grid/sgrid.hh>
typedef Dune :: SGrid< dimgrid, dimworld > GridType;
  #define HAVE_GRIDTYPE 1
#endif

#endif
