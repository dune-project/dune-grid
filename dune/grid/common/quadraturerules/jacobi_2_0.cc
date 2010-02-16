// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#include "config.h"
#include "../quadraturerules.hh"

namespace Dune {

  template Jacobi2QuadratureRule<float, 1>::Jacobi2QuadratureRule(int);
  template Jacobi2QuadratureRule<double, 1>::Jacobi2QuadratureRule(int);

} // namespace
