// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef __ALBERTGRID_ELMEM_CC__
#define __ALBERTGRID_ELMEM_CC__

#include "misc.hh"

namespace AlbertHelp
{

  // for vertices not needed now
  // dim = 2 means, we have a vector for element and edge numbering
#if DUNE_ALBERTA_VERSION < 0x200
  enum { numOfElNumVec = DIM };
#else
  enum { numOfElNumVec = DIM_OF_WORLD };
#endif

  // IndexManagerType defined in albertgrid.hh
  // not thread save !!!
  static IndexManagerType * tmpIndexStack[numOfElNumVec];

  inline static void initIndexManager_elmem_cc(IndexManagerType (&newIm)[numOfElNumVec])
  {
    for(int i=0; i<numOfElNumVec; i++)
    {
      tmpIndexStack[i] = &newIm[i];
      assert(tmpIndexStack[i] != 0);
    }
  }

  inline static void removeIndexManager_elmem_cc(int numOfVec)
  {
    for(int i=0; i<numOfVec; i++) tmpIndexStack[i] = 0;
  }

} // end namespace AlbertHelp

#endif
