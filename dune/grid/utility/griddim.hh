// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRIDDIM_HH
#define DUNE_GRIDDIM_HH

#ifndef GRIDDIM
  #warning "GRIDDIM not defined, defaulting to GRIDDIM=3"
  #define GRIDDIM 3
#endif
#if not (GRIDDIM > 0)
  #error "GRIDDIM must be a positive integer."
#endif

#ifndef WORLDDIM
  #define WORLDDIM GRIDDIM
#endif
#if (WORLDDIM < GRIDDIM)
  #error "WORLDDIM < GRIDDIM does not make sense."
#endif

const int dimgrid = GRIDDIM;
const int dimworld = WORLDDIM;

#endif
