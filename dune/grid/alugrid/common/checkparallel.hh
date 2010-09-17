// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALUGRID_CHECKPARALLEL_HH
#define DUNE_ALUGRID_CHECKPARALLEL_HH

#if HAVE_ALUGRID
#include <alugrid_defineparallel.h>
#endif

#if HAVE_MPI
// if this variable is defined,
// // then parallel version of ALUGrid is compiled
  #if ALU3DGRID_BUILD_FOR_PARALLEL == 0
    #warning "The ALUGrid-library wasn't compiled for parallel usage. Reconfigure\
  using the MPI compiler script or compile Dune without the MPI support!\
  Defaulting to serial ALUGrid!"
    #define ALU3DGRID_PARALLEL 0
  #else
    #define ALU3DGRID_PARALLEL 1
  #endif
#else
  #define ALU3DGRID_PARALLEL 0
#endif

#endif // #ifndef DUNE_ALUGRID_CHECKPARALLEL_HH
