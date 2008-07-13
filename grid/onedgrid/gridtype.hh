// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ONEDGRID_GRIDTYPE_HH
#define DUNE_ONEDGRID_GRIDTYPE_HH

#include <dune/grid/utility/griddim.hh>

#if defined ONEDGRID
  #if HAVE_GRIDTYPE
    #error "Ambiguous definition of GRIDTYPE."
  #endif
  #if (GRIDDIM != 1) || (WORLDDIM != GRIDDIM)
    #error "ONEDGRID is only available for GRIDDIM=1 and WORLDDIM=GRIDDIM."
  #endif

  #include <dune/grid/onedgrid.hh>
typedef Dune :: OneDGrid GridType;
  #define HAVE_GRIDTYPE 1
#endif

#endif
