// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_GRIDDIM_HH
#define DUNE_ALBERTA_GRIDDIM_HH

// only use GRIDDIM when 1, 2 or 3
#if defined GRIDDIM && (GRIDDIM >= 1) && (GRIDDIM <= 3)
  #define DIM GRIDDIM
  #if defined GRIDDIMWORLD
    #define DIM_OF_WORLD GRIDDIMWORLD
  #else
// DIM_OF_WORLD is set to DIM by default
    #define DIM_OF_WORLD GRIDDIM
  #endif
#else
  #ifndef ALBERTA_DIM
    #error "ALBERTA_DIM needed to compile AlbertaGrid! \n"
  #endif

  #ifndef ALBERTA_WORLD_DIM
    #define ALBERTA_WORLD_DIM ALBERTA_DIM
  #endif

  #define DIM ALBERTA_DIM
  #define DIM_OF_WORLD ALBERTA_WORLD_DIM
#endif

#endif
