// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_UGGRIDGEOMETRY_HH
#define DUNE_UGGRIDGEOMETRY_HH

/** \file
 * \brief The UGGridGeometry class and its specializations
 */

#include <dune/common/fmatrix.hh>

#include <dune/geometry/multilineargeometry.hh>

namespace Dune {


  /** \brief Defines the geometry part of a mesh entity.
   * \ingroup UGGrid

     \tparam mydim Dimension of the corresponding reference element
     \tparam coorddim Each corner is a point with coorddim coordinates.

     This version is actually used only for mydim==coorddim. The mydim<coorddim cases
     are in specializations below.
   */
  template<int mydim, int coorddim, class GridImp>
  class UGGridGeometry :
    public GeometryDefaultImplementation <mydim, coorddim, GridImp, UGGridGeometry>
  {
    typedef typename GridImp::ctype UGCtype;

    template <int codim_, int dim_, class GridImp_>
    friend class UGGridEntity;

  public:

    /** \brief Default constructor
     */
    UGGridGeometry()
    {}

    /** \brief Return the element type identifier
     *
     * UGGrid supports triangles and quadrilaterals in 2D, and
     * tetrahedra, pyramids, prisms, and hexahedra in 3D.
     */
    GeometryType type () const;

    //! returns true if type is simplex, false otherwise
    bool affine() const { return type().isSimplex(); }

    //! return the number of corners of this element.
    int corners () const {
      return UG_NS<coorddim>::Corners_Of_Elem(target_);
    }

    //! access to coordinates of corners. Index is the number of the corner
    FieldVector<UGCtype, coorddim> corner (int i) const;

    /** \brief Maps a local coordinate within reference element to
     * global coordinate in element  */
    FieldVector<UGCtype, coorddim> global (const FieldVector<UGCtype, mydim>& local) const;

    /** \brief Maps a global coordinate within the element to a
     * local coordinate in its reference element */
    FieldVector<UGCtype, mydim> local (const FieldVector<UGCtype, coorddim>& global) const;

    /**
       Integration over a general element is done by integrating over the reference element
       and using the transformation from the reference element to the global element as follows:
       \f[\int\limits_{\Omega_e} f(x) dx = \int\limits_{\Omega_{ref}} f(g(l)) A(l) dl \f] where
       \f$g\f$ is the local to global mapping and \f$A(l)\f$ is the integration element.

       For a general map \f$g(l)\f$ involves partial derivatives of the map (surface element of
       the first kind if \f$d=2,w=3\f$, determinant of the Jacobian of the transformation for
       \f$d=w\f$, \f$\|dg/dl\|\f$ for \f$d=1\f$).

       For linear elements, the derivatives of the map with respect to local coordinates
       do not depend on the local coordinates and are the same over the whole element.

       For a structured mesh where all edges are parallel to the coordinate axes, the
       computation is the length, area or volume of the element is very simple to compute.

       Each grid module implements the integration element with optimal efficiency. This
       will directly translate in substantial savings in the computation of finite element
       stiffness matrices.
     */
    UGCtype integrationElement (const FieldVector<UGCtype, mydim>& local) const;

    UGCtype volume() const {

      if (mydim==0)
        return 1;

      // coorddim*coorddim is an upper bound for the number of vertices
      UGCtype* cornerCoords[coorddim*coorddim];
      int n = UG_NS<coorddim>::Corner_Coordinates(target_, cornerCoords);
      return UG_NS<coorddim>::Area_Of_Element(n,
                                              const_cast<const double**>(cornerCoords));
    }

    //! The inverse transpose of the Jacobian matrix of the mapping from the reference element to this element
    FieldMatrix<UGCtype, coorddim,mydim> jacobianInverseTransposed (const FieldVector<UGCtype, mydim>& local) const;
    //! The transpose of the Jacobian matrix of the mapping from the reference element to this element
    FieldMatrix<UGCtype, mydim,coorddim> jacobianTransposed (const FieldVector<UGCtype, mydim>& local) const;


  private:

    /** \brief Init the element with a given UG element */
    void setToTarget(typename UG_NS<coorddim>::template Entity<coorddim-mydim>::T* target)
    {
      target_ = target;
    }

    // in element mode this points to the element we map to
    // in coord_mode this is the element whose reference element is mapped into the father's one
    typename UG_NS<coorddim>::template Entity<coorddim-mydim>::T* target_;

  };



  /****************************************************************/
  /*                                                              */
  /*       Specialization for faces in 3d                         */
  /*                                                              */
  /****************************************************************/

  template<class GridImp>
  class UGGridGeometry<2, 3, GridImp> :
    public MultiLinearGeometry<typename GridImp::ctype, 2, 3>
  {
  public:

    /** \brief Constructor from a given geometry type and a vector of corner coordinates */
    UGGridGeometry(const GeometryType& type, const std::vector<FieldVector<typename GridImp::ctype,3> >& coordinates)
      : MultiLinearGeometry<typename GridImp::ctype, 2, 3>(type, coordinates)
    {}

  };


  /****************************************************************/
  /*                                                              */
  /*       Specialization for edges in 3d                         */
  /*                                                              */
  /****************************************************************/

  template<class GridImp>
  class UGGridGeometry<1, 3, GridImp> :
    public MultiLinearGeometry<typename GridImp::ctype, 1, 3>
  {
  public:

    /** \brief Constructor from a given geometry type and a vector of corner coordinates */
    UGGridGeometry(const GeometryType& type, const std::vector<FieldVector<typename GridImp::ctype,3> >& coordinates)
      : MultiLinearGeometry<typename GridImp::ctype, 1, 3>(type, coordinates)
    {}

  };


  /****************************************************************/
  /*                                                              */
  /*       Specialization for faces in 2d                         */
  /*                                                              */
  /****************************************************************/

  template<class GridImp>
  class UGGridGeometry <1, 2, GridImp> :
    public MultiLinearGeometry<typename GridImp::ctype,1,2>
  {
  public:

    /** \brief Constructor from a given geometry type and a vector of corner coordinates */
    UGGridGeometry(const GeometryType& type, const std::vector<FieldVector<typename GridImp::ctype,2> >& coordinates)
      : MultiLinearGeometry<typename GridImp::ctype, 1, 2>(type, coordinates)
    {}

  };

}  // namespace Dune

#endif
