// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGGRID_LOCALGEOMETRY_HH
#define DUNE_UGGRID_LOCALGEOMETRY_HH

/** \file
 * \brief The UGGridLocalGeometry class
 */

#include <dune/geometry/multilineargeometry.hh>

namespace Dune {


  /** \brief Geometry of an entity embedded in another element, not in the world space
   * \ingroup UGGrid

     \tparam mydim Dimension of the corresponding reference element
     \tparam coorddim Dimension of the coordinate space

     This class is just an adaptor that allows to use MultiLinearGeometry as the implementation
     for the UGGrid local geometries.  We cannot use MultiLinearGeometry directly, because the
     template parameters are different.
   */
  template<int mydim, int coorddim, class GridImp>
  class UGGridLocalGeometry :
    public MultiLinearGeometry<typename GridImp::ctype, mydim, coorddim>
  {
  public:

    /** \brief Constructor from a given geometry type and a vector of corner coordinates */
    UGGridLocalGeometry(const GeometryType& type, const std::vector<FieldVector<typename GridImp::ctype,coorddim> >& coordinates)
      : MultiLinearGeometry<typename GridImp::ctype, mydim, coorddim>(type, coordinates)
    {}

  };

}  // namespace Dune

#endif
