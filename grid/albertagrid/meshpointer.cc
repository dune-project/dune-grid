// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#if HAVE_ALBERTA

#include <dune/grid/albertagrid/meshpointer.hh>

namespace Dune
{

  namespace Alberta
  {

    template<>
    template<>
    int MeshPointer< 1 >::Library< dimWorld >::boundaryCount = 0;

    template<>
    template<>
    const void *MeshPointer< 1 >::Library< dimWorld >::projectionFactory = 0;


    template<>
    template<>
    int MeshPointer< 2 >::Library< dimWorld >::boundaryCount = 0;

    template<>
    template<>
    const void *MeshPointer< 2 >::Library< dimWorld >::projectionFactory = 0;


    template<>
    template<>
    int MeshPointer< 3 >::Library< dimWorld >::boundaryCount = 0;

    template<>
    template<>
    const void *MeshPointer< 3 >::Library< dimWorld >::projectionFactory = 0;

  }

}

#endif // #if HAVE_ALBERTA
