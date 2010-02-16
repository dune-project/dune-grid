// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"
#include "../quadraturerules.hh"

namespace Dune {

  template CubeQuadratureRule<float, 1>::CubeQuadratureRule(int);
  template CubeQuadratureRule<double, 1>::CubeQuadratureRule(int);

} // namespace
