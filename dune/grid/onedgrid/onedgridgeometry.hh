// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ONE_D_GEOMETRY_HH
#define DUNE_ONE_D_GEOMETRY_HH

#include <dune/geometry/axisalignedcubegeometry.hh>

#include <dune/grid/common/geometry.hh>

#include <dune/grid/onedgrid/onedgridentity.hh>

/** \file
 * \brief The OneDGridElement class and its specializations
 */

namespace Dune {

  // forward declaration
  template <int codim, int dim, class GridImp>
  class OneDGridEntity;

  /** \brief Unspecialized class.  Not used for anything */
  template<int mydim, int coorddim, class GridImp>
  class OneDGridGeometry;

  //**********************************************************************
  //
  // --OneDGridGeometry
  /** \brief Defines the geometry part of a vertex.
   * \ingroup OneDGrid
   */
  template<class GridImp>
  class OneDGridGeometry <0, 1, GridImp> :
    public AxisAlignedCubeGeometry<typename GridImp::ctype,0,1>
  {
  public:
    explicit OneDGridGeometry(const FieldVector<typename GridImp::ctype,1>& p)
      : AxisAlignedCubeGeometry<typename GridImp::ctype,0,1>(p)
    {}
  };

  //**********************************************************************
  //
  // --OneDGridGeometry
  /** \brief Defines the geometry part of a mesh entity.
   * \ingroup OneDGrid
   */
  template<int mydim, int coorddim, class GridImp>
  class OneDGridGeometry :
    public AxisAlignedCubeGeometry<typename GridImp::ctype, mydim, coorddim>
  {
  public:
    OneDGridGeometry(const FieldVector<double,1>& left, const FieldVector<double,1>& right)
      : AxisAlignedCubeGeometry<typename GridImp::ctype, mydim, coorddim>(left,right)
    {}
  };

}  // namespace Dune

#endif
