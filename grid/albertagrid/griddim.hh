// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_GRIDDIM_HH
#define DUNE_ALBERTA_GRIDDIM_HH

// only use GRIDDIM when 1, 2 or 3
#if defined GRIDDIM && (GRIDDIM >= 1) && (GRIDDIM <= 3)
  #if defined GRIDDIMWORLD
    #define DIM_OF_WORLD GRIDDIMWORLD
  #else
    #define DIM_OF_WORLD GRIDDIM
  #endif
#else
  #ifndef ALBERTA_DIM
    #error "ALBERTA_DIM needed to compile AlbertaGrid."
  #endif

  #define DIM_OF_WORLD ALBERTA_DIM
#endif

#endif
