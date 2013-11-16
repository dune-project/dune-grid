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

     This class is just an adaptor that allows to use CachedMultiLinearGeometry as the implementation
     for the UGGrid local geometries.  We cannot use CachedMultiLinearGeometry directly, because the
     template parameters are different.
   */
  template<int mydim, int coorddim, class GridImp>
  class UGGridLocalGeometry :
    public CachedMultiLinearGeometry<typename GridImp::ctype, mydim, coorddim>
  {
  public:

    /** \brief Constructor from a given geometry type and a vector of corner coordinates */
    UGGridLocalGeometry(const GeometryType& type, const std::vector<FieldVector<typename GridImp::ctype,coorddim> >& coordinates)
      : CachedMultiLinearGeometry<typename GridImp::ctype, mydim, coorddim>(type, coordinates)
    {}

  };

  namespace FacadeOptions
  {

    /** \brief Switch on the new implementation for the Geometry interface class
     * \deprecated Eventually the new implementation will be hardwired,
     *             and this switch may disappear without prior warning!
     */
    template< int mydim, int cdim, class GridImp>
    struct StoreGeometryReference<mydim,cdim,GridImp,UGGridLocalGeometry>
    {
      static const bool v = false;
    };

  }


}  // namespace Dune

#endif
