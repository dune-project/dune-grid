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

#include <dune/grid/uggrid/uggridgeometry.hh>

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
    public MultiLinearGeometry<typename GridImp::ctype, mydim, coorddim, UGGridGeometryTraits<typename GridImp::ctype>>
  {
  public:
    // inherit constructor from MultiLinearGeometry
    using MultiLinearGeometry<typename GridImp::ctype, mydim, coorddim, UGGridGeometryTraits<typename GridImp::ctype>>::MultiLinearGeometry;

    // factory of uninitialized corner storage used to construct this geometry
    static auto makeCornerStorage(std::size_t count) {
      if constexpr (mydim < 2) // storage when simplex(dim) == cube(dim)
        return std::array<FieldVector<typename GridImp::ctype, coorddim>, (1 << mydim)>{};
      else // storage when simplex(dim) != cube(dim)
        return ReservedVector<FieldVector<typename GridImp::ctype, coorddim>, (1 << mydim)>(count);
    }
  };

}  // namespace Dune

#endif
