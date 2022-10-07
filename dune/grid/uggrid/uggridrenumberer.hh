// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_UGGRID_RENUMBERER_HH
#define DUNE_UGGRID_RENUMBERER_HH

/** \file
    \brief Contains a helper class that encapsulates the vertex numbering conversions
    between UG and DUNE

    The UG numbering conventions can be found in the file 'gm/element.c' of UG.
 */

#include <dune/geometry/type.hh>

namespace Dune {

  /** \brief Generic class which doesn't actually renumber anything.
   *
   *  It is needed as the 1d implementation, used to set up intersections for 2d grids.
   */
  template <int dim>
  class UGGridRenumberer {

  public:

    /** \brief Turn a local vertex number from DUNE numbering to UG numbering
     */
    static int verticesDUNEtoUG(int i, const GeometryType& type) {
      return i;
    }

    /** \brief Turn a local vertex number from UG numbering to DUNE numbering
     */
    static int verticesUGtoDUNE(int i, const GeometryType& type) {
      return i;
    }

  };

  /** \brief DUNE and UG use different local numberings for the subentities of elements.
      This class does the conversions for 2d-grids.

      \todo Is there an efficient and elegant way to remove one of the redundant
      facesUGtoDUNE methods?
   */
  template <>
  class UGGridRenumberer<2> {

  public:

    /** \brief Turn a local edge number from DUNE numbering to UG numbering */
    static int verticesDUNEtoUG(int i, const GeometryType& type) {

      if (type.isCube()) {
        // vertices of a quadrilateral
        const int renumbering[4] = {0, 1, 3, 2};
        return renumbering[i];
      }

      return i;
    }

    /** \brief Turn a local vertex number from UG numbering to DUNE numbering */
    static int verticesUGtoDUNE(int i, const GeometryType& type) {

      if (type.isCube()) {
        // vertices of a quadrilateral
        const int renumbering[4] = {0, 1, 3, 2};
        return renumbering[i];
      }

      return i;
    }

    /** \brief Turn a local face number from DUNE numbering to UG numbering */
    static int facesDUNEtoUG(int i, const GeometryType& type) {

      if (type.isCube()) {

        // faces of a quadrilateral
        const int renumbering[4] = {3, 1, 0, 2};
        return renumbering[i];

      }

      if (type.isSimplex()) {

        // faces of a triangle
        const int renumbering[3] = {0, 2, 1};
        return renumbering[i];
      }

      return i;
    }

    /** \brief Turn a local edge number from DUNE numbering to UG numbering */
    static int edgesDUNEtoUG(int i, const GeometryType& type) {
      return facesDUNEtoUG(i, type);
    }

    /** \brief Turn a local face number from UG numbering to DUNE numbering */
    static int facesUGtoDUNE(int i, const GeometryType& type) {

      if (type.isCube()) {

        // faces of a quadrilateral
        const int renumbering[4] = {2, 1, 3, 0};
        return renumbering[i];

      } else if (type.isSimplex()) {

        // faces of a triangle
        const int renumbering[3] = {0, 2, 1};
        return renumbering[i];
      }

      return i;
    }

    /** \brief Turn a local edge number from UG numbering to DUNE numbering */
    static int edgesUGtoDUNE(int i, const GeometryType& type) {
      return facesUGtoDUNE(i, type);
    }

    /** \brief Turn a local face number from UG numbering to DUNE numbering
     * \param tag The UG way to specify element types.  See the file ug/gm/gm.h for the possible values
     */
    static int facesUGtoDUNE(int i, unsigned int tag) {

      if (tag == UG::D2::QUADRILATERAL) {

        // faces of a quadrilateral
        const int renumbering[4] = {2, 1, 3, 0};
        return renumbering[i];

      } else if (tag == UG::D2::TRIANGLE) {

        // faces of a triangle
        const int renumbering[3] = {0, 2, 1};
        return renumbering[i];
      }

      return i;
    }

    /** \brief Turn a local edge number from UG numbering to DUNE numbering
     * \param tag The UG way to specify element types.  See the file ug/gm/gm.h for the possible values
     */
    static int edgesUGtoDUNE(int i, unsigned int tag) {
      return facesUGtoDUNE(i, tag);
    }

  };

  /** \brief DUNE and UG use different local numberings for the subentities of elements.
      This class does the conversions for 3d-grids.

      \todo Is there an efficient and elegant way to remove one of the redundant
      facesUGtoDUNE methods?
   */
  template <>
  class UGGridRenumberer<3> {

  public:

    /** \brief Turn a local edge number from new DUNE numbering to UG numbering */
    static int verticesDUNEtoUG(int i, const GeometryType& type) {

      if (type.isCube()) {
        const int renumbering[8] = {0, 1, 3, 2, 4, 5, 7, 6};
        return renumbering[i];
      }

      if (type.isPyramid()) {
        const int renumbering[5] = {0, 1, 3, 2, 4};
        return renumbering[i];
      }

      return i;
    }

    /** \brief Turn a local vertex number from UG numbering to DUNE numbering */
    static int verticesUGtoDUNE(int i, const GeometryType& type) {

      if (type.isCube()) {
        const int renumbering[8] = {0, 1, 3, 2, 4, 5, 7, 6};
        return renumbering[i];
      }

      if (type.isPyramid()) {
        const int renumbering[5] = {0, 1, 3, 2, 4};
        return renumbering[i];
      }

      return i;
    }

    /** \brief Turn a local edge number from DUNE numbering to UG numbering */
    static int edgesDUNEtoUG(int i, const GeometryType& type) {

      if (type.isCube()) {

        // edges of a hexahedron
        const int renumbering[12] = {4, 5, 7, 6, 3, 1, 0, 2, 11, 9, 8, 10};
        return renumbering[i];

      }

      if (type.isPrism()) {

        // edges of a prism
        const int renumbering[9] = {3, 4, 5, 0, 2, 1, 6, 8, 7};
        return renumbering[i];

      }

      if (type.isPyramid()) {

        // edges of a pyramid
        const int renumbering[8] = {3, 1, 0, 2, 4, 5, 7, 6};
        return renumbering[i];

      }

      if (type.isSimplex()) {

        // edges of a tetrahedon
        const int renumbering[6] = {0, 2, 1, 3, 4, 5};
        return renumbering[i];
      }

      return i;
    }

    /** \brief Turn a local edge number from UG numbering to DUNE numbering */
    static int edgesUGtoDUNE(int i, const GeometryType& type) {

      if (type.isCube()) {

        // edges of a hexahedron
        const int renumbering[12] = {6, 5, 7, 4, 0, 1, 3, 2, 10, 9, 11, 8};
        return renumbering[i];

      }

      if (type.isPrism()) {

        // edges of a prism
        const int renumbering[9] = {3, 5, 4, 0, 1, 2, 6, 8, 7};
        return renumbering[i];

      }

      if (type.isPyramid()) {

        // edges of a pyramid
        const int renumbering[8] = {2, 1, 3, 0, 4, 5, 7, 6};
        return renumbering[i];

      }

      if (type.isSimplex()) {

        // edges of a tetrahedon
        const int renumbering[6] = {0, 2, 1, 3, 4, 5};
        return renumbering[i];
      }

      return i;
    }

    /** \brief Turn a local edge number from UG numbering to DUNE numbering
     * \param tag The UG way to specify element types.  See the file ug/gm/gm.h for the possible values
     */
    static int edgesUGtoDUNE(int i, unsigned int tag) {

      if (tag == UG::D3::HEXAHEDRON) {

        // edges of a hexahedron
        const int renumbering[12] = {6, 5, 7, 4, 0, 1, 3, 2, 10, 9, 11, 8};
        return renumbering[i];

      } if (tag == UG::D3::PRISM) {

        // edges of a prism
        const int renumbering[9] = {3, 5, 4, 0, 1, 2, 6, 8, 7};
        return renumbering[i];

      } if (tag == UG::D3::PYRAMID) {

        // edges of a pyramid
        const int renumbering[8] = {2, 1, 3, 0, 4, 5, 7, 6};
        return renumbering[i];

      } else if (tag == UG::D3::TETRAHEDRON) {

        // edges of a tetrahedon
        const int renumbering[6] = {0, 2, 1, 3, 4, 5};
        return renumbering[i];
      }

      return i;
    }

    /** \brief Turn a local face number from DUNE numbering to UG numbering */
    static int facesDUNEtoUG(int i, const GeometryType& type) {

      if (type.isCube()) {

        // faces of a hexahedron
        const int renumbering[6] = {4, 2, 1, 3, 0, 5};
        return renumbering[i];

      }

      if (type.isPrism()) {

        // faces of a prism
        const int renumbering[5] = {1, 3, 2, 0, 4};
        return renumbering[i];

      }

      if (type.isPyramid()) {

        // faces of a pyramid
        const int renumbering[5] = {0, 4, 2, 1, 3};
        return renumbering[i];

      }

      if (type.isSimplex()) {

        // faces of a tetrahedon
        const int renumbering[4] = {0, 3, 2, 1};
        return renumbering[i];
      }

      return i;
    }

    /** \brief Turn a local face number from UG numbering to DUNE numbering */
    static int facesUGtoDUNE(int i, const GeometryType& type) {

      if (type.isCube()) {

        // faces of a hexahedron
        const int renumbering[6] = {4, 2, 1, 3, 0, 5};
        return renumbering[i];

      }

      if (type.isPrism()) {

        // faces of a hexahedron
        const int renumbering[5] = {3, 0, 2, 1, 4};
        return renumbering[i];

      }

      if (type.isPyramid()) {

        // faces of a hexahedron
        const int renumbering[5] = {0, 3, 2, 4, 1};
        return renumbering[i];

      }

      if (type.isSimplex()) {

        // faces of a tetrahedon
        const int renumbering[4] = {0, 3, 2, 1};
        return renumbering[i];
      }

      return i;
    }

    /** \brief Turn a local face number from UG numbering to DUNE numbering
     * \param tag The UG way to specify element types.  See the file ug/gm/gm.h for the possible values
     */
    static int facesUGtoDUNE(int i, unsigned int tag) {

      if (tag == UG::D3::HEXAHEDRON) {

        // faces of a hexahedron
        const int renumbering[6] = {4, 2, 1, 3, 0, 5};
        return renumbering[i];

      } if (tag == UG::D3::PRISM) {

        // faces of a hexahedron
        const int renumbering[5] = {3, 0, 2, 1, 4};
        return renumbering[i];

      } if (tag == UG::D3::PYRAMID) {

        // faces of a hexahedron
        const int renumbering[5] = {0, 3, 2, 4, 1};
        return renumbering[i];

      } else if (tag == UG::D3::TETRAHEDRON) {

        // faces of a tetrahedon
        const int renumbering[4] = {0, 3, 2, 1};
        return renumbering[i];
      }

      return i;
    }

  };
}

#endif
