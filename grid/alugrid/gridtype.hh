// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_GRIDTYPE_HH
#define DUNE_ALUGRID_GRIDTYPE_HH

#include <dune/grid/utility/griddim.hh>

#if defined ALUGRID_CUBE
  #if HAVE_GRIDTYPE
    #error "Ambiguous definition of GRIDTYPE."
  #endif
  #if not HAVE_ALUGRID
    #error "ALUGRID_CUBE defined but no ALUGRID version found!"
  #endif
  #if (GRIDDIM != 3) || (WORLDDIM != GRIDDIM)
    #error "ALUGRID_CUBE is only available for GRIDDIM=3 and WORLDDIM=GRIDDIM."
  #endif

  #include <dune/grid/alugrid.hh>
typedef Dune :: ALUCubeGrid< dimgrid, dimworld > GridType;
  #define HAVE_GRIDTYPE 1
#endif

#if defined ALUGRID_SIMPLEX
  #if HAVE_GRIDTYPE
    #error "Ambiguous definition of GRIDTYPE."
  #endif
  #if not HAVE_ALUGRID
    #error "ALUGRID_SIMPLEX defined but no ALUGRID version found!"
  #endif
  #if (GRIDDIM < 2) || (GRIDDIM > 3)
    #error "ALUGRID_SIMPLEX is only available for GRIDDIM=2 and GRIDDIM=3."
  #endif
  #if (WORLDDIM != GRIDDIM)
    #error "ALUGRID_SIMPLEX is only available for WORLDDIM=GRIDDIM."
  #endif

  #include <dune/grid/alugrid.hh>
typedef Dune :: ALUSimplexGrid< dimgrid, dimworld > GridType;
  #define HAVE_GRIDTYPE 1
#endif

#if defined ALUGRID_CONFORM
  #if HAVE_GRIDTYPE
    #error "Ambiguous definition of GRIDTYPE."
  #endif
  #if not HAVE_ALUGRID
    #error "ALUGRID_CONFORM defined but no ALUGRID version found!"
  #endif
  #if (GRIDDIM != 2) || (WORLDDIM != GRIDDIM)
    #error "ALUGRID_CONFORM is only available for GRIDDIM=2 and WORLDDIM=GRIDDIM."
  #endif
  #include <dune/grid/alugrid.hh>
typedef Dune :: ALUConformGrid< dimgrid, dimworld > GridType;
  #define HAVE_GRIDTYPE 1
#endif

#endif
