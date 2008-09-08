// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ALBERTA_EXCEPTIONS_HH
#define DUNE_ALBERTA_EXCEPTIONS_HH

#include <dune/common/exceptions.hh>

namespace Dune
{
  // own exception classes
  class AlbertaError   : public Exception {};
  class AlbertaIOError : public IOError {};
}

#endif
