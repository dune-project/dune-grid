// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <algorithm>

#include <dune/grid/uggrid.hh>
#include <dune/grid/uggrid/uggridgeometry.hh>

namespace Dune {

///////////////////////////////////////////////////////////////////////
//
// General implementation of UGGridGeometry mydim-> coorddim
//
///////////////////////////////////////////////////////////////////////


template< int mydim, int coorddim, class GridImp>
GeometryType UGGridGeometry<mydim,coorddim,GridImp>::type() const
{
  switch (mydim)
  {
  case 0 : return GeometryTypes::vertex;
  case 1 : return GeometryTypes::line;
  case 2 :

    switch (UG_NS<coorddim>::Tag(target_)) {
    case UG::D2::TRIANGLE :
      return GeometryTypes::triangle;
    case UG::D2::QUADRILATERAL :
      return GeometryTypes::quadrilateral;
    default :
      DUNE_THROW(GridError, "UGGridGeometry::type():  ERROR:  Unknown type "
                 << UG_NS<coorddim>::Tag(target_) << " found!");
    }

  case 3 :
    switch (UG_NS<coorddim>::Tag(target_)) {

    case UG::D3::TETRAHEDRON :
      return GeometryTypes::tetrahedron;
    case UG::D3::PYRAMID :
      return GeometryTypes::pyramid;
    case UG::D3::PRISM :
      return GeometryTypes::prism;
    case UG::D3::HEXAHEDRON :
      return GeometryTypes::hexahedron;
    default :
      DUNE_THROW(GridError, "UGGridGeometry::type():  ERROR:  Unknown type "
                 << UG_NS<coorddim>::Tag(target_) << " found!");

    }
  }

}


template<int mydim, int coorddim, class GridImp>
FieldVector<typename GridImp::ctype, coorddim> UGGridGeometry<mydim,coorddim,GridImp>::
corner(int i) const
{
  // This geometry is a vertex
  if (mydim==0) {
    FieldVector<typename GridImp::ctype, coorddim> result;
    for (size_t j=0; j<coorddim; j++)
      // The cast is here to make the code compile even when target_ is not a node
      result[j] = ((typename UG_NS<coorddim>::Node*)target_)->myvertex->iv.x[j];
    return result;
  }

  // ////////////////////////////////
  //  This geometry is an element
  // ////////////////////////////////
  assert(mydim==coorddim);

  i = UGGridRenumberer<mydim>::verticesDUNEtoUG(i,type());

  FieldVector<typename GridImp::ctype, coorddim> result;
  for (size_t j=0; j<coorddim; j++)
    // The cast is here to make the code compile even when target_ is not an element
    result[j] = UG_NS<coorddim>::Corner(((typename UG_NS<coorddim>::Element*)target_),i)->myvertex->iv.x[j];
  return result;
}

template< int mydim, int coorddim, class GridImp>
FieldVector<typename GridImp::ctype, coorddim> UGGridGeometry<mydim,coorddim,GridImp>::
global(const FieldVector<UGCtype, mydim>& local) const
{
  FieldVector<UGCtype, coorddim> globalCoord(0.0);

  // we are an actual element in UG
  // coorddim*coorddim is an upper bound for the number of vertices
  UGCtype* cornerCoords[coorddim*coorddim];
  int n = UG_NS<coorddim>::Corner_Coordinates(target_, cornerCoords);

  // Actually do the computation
  UG_NS<coorddim>::Local_To_Global(n, cornerCoords, local, globalCoord);

  return globalCoord;
}


// Maps a global coordinate within the element to a
// local coordinate in its reference element
template< int mydim, int coorddim, class GridImp>
FieldVector<typename GridImp::ctype, mydim> UGGridGeometry<mydim,coorddim, GridImp>::
local (const FieldVector<typename GridImp::ctype, coorddim>& global) const
{
  FieldVector<UGCtype, mydim> result;

  // do nothing for a vertex
  if (mydim==0)
    return result;

  // coorddim*coorddim is an upper bound for the number of vertices
  UGCtype* cornerCoords[coorddim*coorddim];
  int n = UG_NS<coorddim>::Corner_Coordinates(target_, cornerCoords);

  // Actually do the computation
  /** \todo Why is this const_cast necessary? */
  UG_NS<coorddim>::GlobalToLocal(n, const_cast<const double**>(cornerCoords), &global[0], &result[0]);

  return result;
}


template< int mydim, int coorddim, class GridImp>
typename GridImp::ctype UGGridGeometry<mydim,coorddim,GridImp>::
integrationElement (const FieldVector<typename GridImp::ctype, mydim>& local) const
{
  if (mydim==0)
    return 1;
  else
    /** \todo No need to recompute the determinant every time on a simplex */
    return std::abs(1/jacobianInverseTransposed(local).determinant());
}


template< int mydim, int coorddim, class GridImp>
FieldMatrix<typename GridImp::ctype, coorddim,mydim> UGGridGeometry<mydim,coorddim, GridImp>::
jacobianInverseTransposed (const FieldVector<typename GridImp::ctype, mydim>& local) const
{
  FieldMatrix<UGCtype,coorddim,mydim> jIT;

  // compile array of pointers to corner coordinates
  // coorddim*coorddim is an upper bound for the number of vertices
  UGCtype* cornerCoords[coorddim*coorddim];
  int n = UG_NS<coorddim>::Corner_Coordinates(target_, cornerCoords);

  // compute the transformation onto the reference element (or vice versa?)
  UG_NS<coorddim>::Transformation(n, cornerCoords, local, jIT);

  return jIT;
}
template< int mydim, int coorddim, class GridImp>
FieldMatrix<typename GridImp::ctype, mydim,coorddim> UGGridGeometry<mydim,coorddim, GridImp>::
jacobianTransposed (const FieldVector<typename GridImp::ctype, mydim>& local) const
{
  FieldMatrix<UGCtype,mydim,coorddim> jac;

  // compile array of pointers to corner coordinates
  // coorddim*coorddim is an upper bound for the number of vertices
  UGCtype* cornerCoords[coorddim*coorddim];
  int n = UG_NS<coorddim>::Corner_Coordinates(target_, cornerCoords);

  // compute the transformation onto the reference element (or vice versa?)
  UG_NS<coorddim>::JacobianTransformation(n, cornerCoords, local, jac);

  return jac;
}


/////////////////////////////////////////////////////////////////////////////////
//   Explicit template instantiations
/////////////////////////////////////////////////////////////////////////////////

template class UGGridGeometry<0,2, const UGGrid<2> >;
template class UGGridGeometry<2,2, const UGGrid<2> >;

template class UGGridGeometry<0,3, const UGGrid<3> >;
template class UGGridGeometry<3,3, const UGGrid<3> >;

} /* namespace Dune */
