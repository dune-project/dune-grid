// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <iostream>

#include "../matrix.hh"
#include "../geometrytraits.hh"

using namespace Dune;

int main()
{
  // Test whether I can compute the determinant of A A^T of a nearly singular matrix.
  // This particular matrix makes the detAAT method abort as of dune-grid revision 6631.
  FieldMatrix<double,2,2> A;
  A[0][0] =  0.099999999999999867;
  A[0][1] = -0.010000000000002118;
  A[1][0] =  0.099999999999999867;
  A[1][1] = -0.0099999999999998979;

  double detAAT = GenericGeometry::MatrixHelper<GenericGeometry::DuneCoordTraits<double> >::detAAT<2, 2>(A);
  std::cout << "detAAT = " << detAAT << std::endl;
}
