// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <dune/grid/uggrid/p2geometry.hh>

template <int mydim, int coorddim, class ctype>
Dune::FieldVector<ctype,coorddim> Dune::PiecewiseQuadraticGeometry<mydim,coorddim,ctype>
::global(const FieldVector<ctype,mydim>& local,
         const array<FieldVector<ctype,coorddim>, maxNodes>& nodes,
         const GeometryType& type) {

  Dune::FieldVector<ctype,coorddim> result(0);

  const UGShapeFunctions::LagrangeShapeFunctionSet<ctype, mydim>& sfs
    = UGShapeFunctions::LagrangeShapeFunctions<ctype, mydim>::general(type,2);

  for (int i=0; i<sfs.size(); i++)
    result.axpy(sfs[i].evaluateFunction(0,local), nodes[i]);

  return result;
}

template <int mydim, int coorddim, class ctype>
void Dune::PiecewiseQuadraticGeometry<mydim,coorddim,ctype>
::jacobianInverseTransposed(const FieldVector<ctype,mydim>& local,
                            const array<FieldVector<ctype,coorddim>, maxNodes>& nodes,
                            const GeometryType& type,
                            FieldMatrix<ctype,mydim,mydim>& matrix) {

  matrix = 0;

  const UGShapeFunctions::LagrangeShapeFunctionSet<ctype, mydim>& sfs
    = UGShapeFunctions::LagrangeShapeFunctions<ctype, mydim>::general(type,2);

  for (int i=0; i<sfs.size(); i++) {

    for (int j=0; j<mydim; j++)
      for (int k=0; k<mydim; k++)
        matrix[k][j] += sfs[i].evaluateDerivative(0,k,local) * nodes[i][j];

  }

  matrix.invert();

  //std::cout << "P2:\n" << matrix << std::endl;
}

// Explicit template instantiations
template class Dune::PiecewiseQuadraticGeometry<1,2,double>;
template class Dune::PiecewiseQuadraticGeometry<2,2,double>;
template class Dune::PiecewiseQuadraticGeometry<1,3,double>;
template class Dune::PiecewiseQuadraticGeometry<2,3,double>;
template class Dune::PiecewiseQuadraticGeometry<3,3,double>;
