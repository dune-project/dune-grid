// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <iomanip>

#include <dune/common/gmpfield.hh>

#include "../matrix.hh"
#include "../geometrytraits.hh"

using namespace Dune;

int main()
{
  typedef double Field;
  //typedef long double Field;
  //typedef GMPField< 72 > Field;
  //typedef GMPField< 160 > Field;

  typedef GenericGeometry::MatrixHelper< GenericGeometry::DuneCoordTraits< Field > > MatrixHelper;

  // Test whether I can compute the square root of the determinant of A A^T of a nearly singular matrix.
  // This particular matrix makes the sqrtDetAAT method abort in dune-grid revision 6631.
  FieldMatrix< Field, 2, 2 > A;
  A[0][0] =  0.099999999999999867;
  A[0][1] = -0.010000000000002118;
  A[1][0] =  0.099999999999999867;
  A[1][1] = -0.0099999999999998979;

  std::cout << std::scientific << std::setprecision( 20 );

  Field sqrtDetAAT = MatrixHelper::sqrtDetAAT< 2, 2 >( A );
  std::cout << "sqrtDetAAT = " << sqrtDetAAT << std::endl;

  FieldMatrix< Field, 2, 2 > invA;
  Field detA = MatrixHelper::rightInvA< 2, 2 >( A, invA );
  std::cout << "detA = " << detA << std::endl;
  std::cout << "invA = [ " << invA[ 0 ] << ", " << invA[ 1 ] << " ]" << std::endl;

  // Lets do the same crap for a non-square matrix.
  FieldMatrix< Field, 2, 3 > B;
  B[0][0] =  0.099999999999999867;
  B[0][1] = -0.010000000000002118;
  B[0][2] = 0;
  B[1][0] =  0.099999999999999867;
  B[1][1] = -0.0099999999999998979;
  B[1][2] = 0;

  std::cout << std::scientific << std::setprecision( 20 );

  Field sqrtDetBBT = MatrixHelper::sqrtDetAAT< 2, 3 >( B );
  std::cout << "sqrtDetBBT = " << sqrtDetBBT << std::endl;

  FieldMatrix< Field, 3, 2 > invB;
  Field detB = MatrixHelper::rightInvA< 2, 3 >( B, invB );
  std::cout << "detB = " << detB << std::endl;
  std::cout << "invB = [ " << invB[ 0 ] << ", " << invB[ 1 ] << invB[ 2 ] << " ]" << std::endl;

}
