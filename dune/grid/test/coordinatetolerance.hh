// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_COORDINATETOLERANCE_HH
#define DUNE_GRID_TEST_COORDINATETOLERANCE_HH

#include <algorithm>
#include <limits>

#include <dune/common/float_cmp.hh>

#include <dune/grid/common/rangegenerators.hh>

/** \file
 *  \brief Comparison of global coordinates in the grid checks
 *
 *  Coordinates computed in different ways (on different processes, from
 *  different elements, or after subtracting a period) differ by round-off
 *  of the size of the domain times the machine epsilon. A relative
 *  comparison (the default style of FloatCmp) fails for coordinates near
 *  zero, where it becomes an exact comparison, and a fixed absolute
 *  tolerance fails for large domains. Hence, coordinates are compared with
 *  an absolute tolerance that is relative to the size of the domain.
 */

namespace Dune {
  namespace GridCheck {

    /** \brief Absolute tolerance for comparing global coordinates of a grid view
     *
     *  Returns relTol times the diameter of the bounding box of the element
     *  corners stored on this process (or relTol if there are no elements).
     */
    template<class GridView>
    typename GridView::ctype coordinateTolerance (const GridView& gridView,
                                                  typename GridView::ctype relTol = 1e-8)
    {
      using ctype = typename GridView::ctype;
      using Coordinate = typename GridView::template Codim<0>::Geometry::GlobalCoordinate;

      Coordinate lower(std::numeric_limits<ctype>::max());
      Coordinate upper(std::numeric_limits<ctype>::lowest());
      bool empty = true;
      for (const auto& element : elements(gridView))
      {
        const auto& geometry = element.geometry();
        for (int c = 0; c < geometry.corners(); ++c)
        {
          const Coordinate x = geometry.corner(c);
          for (std::size_t i = 0; i < x.size(); ++i)
          {
            lower[i] = std::min(lower[i], x[i]);
            upper[i] = std::max(upper[i], x[i]);
          }
        }
        empty = false;
      }
      return empty ? relTol : relTol * (upper - lower).two_norm();
    }

    /** \brief Compare two coordinates (or coordinate components) with an absolute tolerance
     *
     *  \param tol  absolute tolerance, see coordinateTolerance()
     */
    template<class T>
    bool coordinatesEqual (const T& a, const T& b, typename FloatCmp::EpsilonType<T>::Type tol)
    {
      return FloatCmp::eq<T, FloatCmp::absolute>(a, b, tol);
    }

  } // end namespace GridCheck
} // end namespace Dune

#endif // DUNE_GRID_TEST_COORDINATETOLERANCE_HH
