// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_GRIDDIM_HH
#define DUNE_ALBERTA_GRIDDIM_HH

// only use GRIDDIM when 2 or 3
#if defined GRIDDIM && (GRIDDIM >= 2) && (GRIDDIM <= 3)
  #undef ALBERTA_DIM
  #undef ALBERTA_WORLD_DIM

  #define ALBERTA_DIM GRIDDIM
  #if defined GRIDDIMWORLD
    #define ALBERTA_WORLD_DIM GRIDDIMWORLD
  #endif
#endif

#ifndef ALBERTA_DIM
  #error "ALBERTA_DIM needed to compile AlbertaGrid! \n"
#endif

#ifndef ALBERTA_WORLD_DIM
  #define ALBERTA_WORLD_DIM ALBERTA_DIM
#endif

#define DIM ALBERTA_DIM
#define DIM_OF_WORLD ALBERTA_WORLD_DIM

#endif
