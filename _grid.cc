// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <dune/corepy/pybind11/pybind11.h>

PYBIND11_PLUGIN( _grid )
{
  pybind11::module module( "_grid" );

  return module.ptr();
}
