// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_EXCEPTIONS_HH
#define DUNE_GRID_COMMON_EXCEPTIONS_HH

#include <dune/common/exceptions.hh>

namespace Dune
{

  // GridError
  // ---------

  /** \brief Base class for exceptions in Dune grid modules
   */
  class GridError
    : public Exception
  {};

}

#endif // #ifndef DUNE_GRID_COMMON_EXCEPTIONS_HH
