// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IDENTITYGRIDGEOMETRY_HH
#define DUNE_IDENTITYGRIDGEOMETRY_HH

/** \file
 * \brief The IdentityGridGeometry class and its specializations
 */

#include <dune/common/fmatrix.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/geometry.hh>

namespace Dune {

  template<int mydim, int coorddim, class GridImp>
  class IdentityGridGeometry :
    public GeometryDefaultImplementation <mydim, coorddim, GridImp, IdentityGridGeometry>
  {
  private:

    typedef typename GridImp::ctype ctype;


  public:

    // The codimension of this entitypointer wrt the host grid
    constexpr static int CodimInHostGrid = GridImp::HostGridType::dimension - mydim;
    constexpr static int DimensionWorld = GridImp::HostGridType::dimensionworld;

    // select appropriate hostgrid geometry via typeswitch
    typedef typename GridImp::HostGridType::Traits::template Codim<CodimInHostGrid>::Geometry HostGridGeometryType;
    typedef typename GridImp::HostGridType::Traits::template Codim<CodimInHostGrid>::Geometry HostGridLocalGeometryType;

    typedef typename std::conditional<coorddim==DimensionWorld, HostGridGeometryType, HostGridLocalGeometryType>::type HostGridGeometry;

    //! type of jacobian transposed
    typedef typename HostGridGeometryType::JacobianInverseTransposed JacobianInverseTransposed;
    typedef typename HostGridGeometryType::JacobianTransposed JacobianTransposed;


    /** constructor from host geometry
     */
    IdentityGridGeometry(const HostGridGeometry& hostGeometry)
      : hostGeometry_(hostGeometry)
    {}


    /** \brief Return the element type identifier
     */
    GeometryType type () const {
      return hostGeometry_.type();
    }

    // return wether we have an affine mapping
    bool affine() const {
      return hostGeometry_.affine();
    }

    //! return the number of corners of this element. Corners are numbered 0...n-1
    int corners () const {
      return hostGeometry_.corners();
    }


    //! access to coordinates of corners. Index is the number of the corner
    const FieldVector<ctype, coorddim> corner (int i) const {
      return hostGeometry_.corner(i);
    }


    /** \brief Maps a local coordinate within reference element to
     * global coordinate in element  */
    FieldVector<ctype, coorddim> global (const FieldVector<ctype, mydim>& local) const {
      return hostGeometry_.global(local);
    }

    /** \brief Return the transposed of the Jacobian
     */
    JacobianTransposed
    jacobianTransposed ( const FieldVector<ctype, mydim>& local ) const {
      return hostGeometry_.jacobianTransposed(local);
    }

    /** \brief Maps a global coordinate within the element to a
     * local coordinate in its reference element */
    FieldVector<ctype, mydim> local (const FieldVector<ctype, coorddim>& global) const {
      return hostGeometry_.local(global);
    }


    //! Returns true if the point is in the current element
    bool checkInside(const FieldVector<ctype, mydim> &local) const {
      return hostGeometry_.checkInside(local);
    }


    /**
     */
    ctype integrationElement (const FieldVector<ctype, mydim>& local) const {
      return hostGeometry_.integrationElement(local);
    }


    //! The Jacobian matrix of the mapping from the reference element to this element
    JacobianInverseTransposed jacobianInverseTransposed (const FieldVector<ctype, mydim>& local) const {
      return hostGeometry_.jacobianInverseTransposed(local);
    }


    HostGridGeometry hostGeometry_;

  };

}  // namespace Dune

#endif
