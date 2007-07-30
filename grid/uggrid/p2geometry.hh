// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PIECEWISE_QUADRATIC_GEOMETRY_HH
#define DUNE_PIECEWISE_QUADRATIC_GEOMETRY_HH

#include <dune/common/geometrytype.hh>
#include <dune/common/fixedarray.hh>

#include "lagrangeshapefunctions.hh"

namespace Dune {

  template <int mydim, int coorddim, class ctype>
  class PiecewiseQuadraticGeometry {

    enum {maxNodes = Power_m_p<3,mydim>::power };

  public:

    static FieldVector<ctype,coorddim> global(const FieldVector<ctype,mydim>& local,
                                              const array<FieldVector<ctype,coorddim>, maxNodes>& nodes,
                                              const GeometryType& type);

    static void jacobianInverseTransposed(const FieldVector<ctype,mydim>& local,
                                          const array<FieldVector<ctype,coorddim>, maxNodes>& nodes,
                                          const GeometryType& type,
                                          FieldMatrix<ctype, mydim,mydim>& matrix);

  };

}


#endif
