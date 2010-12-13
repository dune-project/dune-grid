// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGGRID_RENUMBERER_HH
#define DUNE_UGGRID_RENUMBERER_HH

/** \file
    \brief Contains a helper class that encapsulates the vertex numbering conversions
    between UG and DUNE

    The UG numbering conventions can be found in the file 'gm/element.c' of UG.
 */

#include <dune/common/geometrytype.hh>

#include <dune/grid/genericgeometry/conversion.hh>

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

    /** \brief Turn a local edge number from DUNE numbering to UG numbering */
    static int edgesDUNEtoUG(int i, const GeometryType& type) {

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
    static int facesDUNEtoUG(int i, const GeometryType& type) {
      // faces are edges in 2d
      return edgesDUNEtoUG(i,type);
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

    /** \brief Turn a local face number from UG numbering to DUNE numbering */
    static int facesUGtoDUNE(int i, int nSides) {

      if (nSides == 4) {

        // faces of a quadrilateral
        const int renumbering[4] = {2, 1, 3, 0};
        return renumbering[i];

      } else if (nSides == 3) {

        // faces of a triangle
        const int renumbering[3] = {0, 2, 1};
        return renumbering[i];
      }

      return i;
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
      typedef Dune::GenericGeometry::MapNumberingProvider<3> Numbering;
      const int j = Numbering::generic2dune( type.id(), i, 2 );

      if (type.isCube()) {
        DUNE_THROW(NotImplemented, "edgesDUNEtoUG for hexahedra");
        // edges of a hexahedron
        const int renumbering[6] = {4, 2, 1, 3, 0, 5};
        return renumbering[j];

      }

      return j;
    }

    /** \brief Turn a local edge number from DUNE numbering to UG numbering */
    static int facesDUNEtoUG(int i, const GeometryType& type) {
      typedef Dune::GenericGeometry::MapNumberingProvider<3> Numbering;
      const int j = Numbering::generic2dune( type.id(), i, 1 );

      if (type.isCube()) {

        // faces of a hexahedron
        const int renumbering[6] = {4, 2, 1, 3, 0, 5};
        return renumbering[j];

      }

      if (type.isSimplex()) {

        // faces of a tetrahedon
        const int renumbering[4] = {1, 2, 3, 0};
        return renumbering[j];
      }

      return j;
    }

    /** \brief Turn a local face number from UG numbering to DUNE numbering */
    static int facesUGtoDUNE(int i, const GeometryType& type) {

      int j = i;

      if (type.isCube()) {

        // faces of a hexahedron
        const int renumbering[6] = {4, 2, 1, 3, 0, 5};
        j = renumbering[i];

      } else if (type.isSimplex()) {

        // faces of a tetrahedon
        const int renumbering[4] = {3, 0, 1, 2};
        j =  renumbering[i];
      }

      typedef Dune::GenericGeometry::MapNumberingProvider<3> Numbering;
      return Numbering::dune2generic( type.id(), j, 1 );
    }

    /** \brief Turn a local face number from UG numbering to DUNE numbering */
    static int facesUGtoDUNE(int i, int nSides) {

      if (nSides==6) {

        // faces of a hexahedron
        const int renumbering[6] = {4, 2, 1, 3, 0, 5};
        return renumbering[i];

      } else if (nSides==4) {

        // faces of a tetrahedon
        const int renumbering[4] = {3, 0, 1, 2};
        return renumbering[i];
      }

      return i;
    }

  };
}

#endif
